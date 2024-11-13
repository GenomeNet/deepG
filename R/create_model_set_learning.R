#' @title Create LSTM/CNN network for combining multiple sequences 
#' 
#' @description Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers.
#' Input is a 4D tensor, where axis correspond to:
#' \enumerate{
#'    \item batch size
#'    \item number of samples in one batch
#'    \item length of one sample
#'    \item size of vocabulary
#' }
#' After LSTM/CNN part all representations get aggregated by summation.
#' Can be used to make single prediction for combination of multiple input sequences. Architecture
#' is equivalent to \link{create_model_lstm_cnn_multi_input} but instead of multiple input layers with 3D input, 
#' input here in one 4D tensor.  
#'     
#' @inheritParams create_model_lstm_cnn
#' @param samples_per_target Number of samples to combine for one target.
#' @param aggregation_method At least one of the options `"sum", "mean", "max"`.
#' @param gap_time_dist Pooling or flatten method after last time distribution wrapper. Same options as for `flatten_method` argument
#' in \link{create_model_transformer} function.
#' @param lstm_time_dist Vector containing number of units per LSTM cell. Applied after time distribution part. 
#' @param transformer_args List of arguments for transformer blocks; see \link{layer_transformer_block_wrapper}.
#' Additionally, list can contain `pool_flatten` argument to apply global pooling or flattening after last transformer block (same options
#' as `flatten_method` argument in \link{create_model_transformer} function).
#' @examplesIf reticulate::py_module_available("tensorflow")
#' 
#' # Examples needs keras attached to run 
#' maxlen <- 50
#' \donttest{
#' library(keras)
#' create_model_lstm_cnn_time_dist(
#'   maxlen = maxlen,
#'   vocabulary_size = 4,
#'   samples_per_target = 7,
#'   kernel_size = c(10, 10),
#'   filters = c(64, 128),
#'   pool_size = c(2, 2),
#'   layer_lstm = c(32),
#'   aggregation_method = c("max"),
#'   layer_dense = c(64, 2),
#'   learning_rate = 0.001)
#'  }
#'  
#' @returns A keras model with time distribution wrapper applied to LSTM and CNN layers.   
#' @export
create_model_lstm_cnn_time_dist <- function(
    maxlen = 50,
    dropout_lstm = 0,
    recurrent_dropout_lstm = 0,
    layer_lstm = NULL,
    layer_dense = c(4),
    solver = "adam",
    learning_rate = 0.001,
    vocabulary_size = 4,
    bidirectional = FALSE,
    stateful = FALSE,
    batch_size = NULL,
    compile = TRUE,
    kernel_size = NULL,
    filters = NULL,
    strides = NULL,
    pool_size = NULL,
    padding = "same",
    dilation_rate = NULL,
    gap_time_dist = NULL,
    use_bias = TRUE,
    zero_mask = FALSE,
    label_smoothing = 0,
    label_noise_matrix = NULL,
    last_layer_activation = "softmax",
    loss_fn = "categorical_crossentropy",
    auc_metric = FALSE,
    f1_metric = FALSE,
    samples_per_target,
    batch_norm_momentum = 0.99,
    verbose = TRUE,
    model_seed = NULL,
    aggregation_method = NULL, 
    transformer_args = NULL,
    lstm_time_dist = NULL,
    mixed_precision = FALSE,
    bal_acc = FALSE,
    mirrored_strategy = NULL) {
  
  if (mixed_precision) tensorflow::tf$keras$mixed_precision$set_global_policy("mixed_float16")
  
  if (is.null(mirrored_strategy)) mirrored_strategy <- ifelse(count_gpu() > 1, TRUE, FALSE)
  if (mirrored_strategy) {
    mirrored_strategy <- tensorflow::tf$distribute$MirroredStrategy()
    with(mirrored_strategy$scope(), { 
      argg <- as.list(environment())
      argg$mirrored_strategy <- FALSE
      model <- do.call(create_model_lstm_cnn_time_dist, argg)
    })
    return(model)
  }
  
  #stopifnot(aggregation_method %in% c("sum", "lstm", "lstm_sum"))
  
  layer_dense <- as.integer(layer_dense)
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  num_output_layers <- 1
  num_input_layers <- 1
  
  num_targets <- layer_dense[length(layer_dense)]
  layers.lstm <- length(layer_lstm)
  use.cnn <- ifelse(!is.null(kernel_size), TRUE, FALSE)
  
  if (!is.null(layer_lstm)) {
    stopifnot(length(layer_lstm) == 1 | (length(layer_lstm) ==  layers.lstm))
  }
  
  if (layers.lstm == 0 & !use.cnn) {
    stop("Model does not use LSTM or CNN layers.")
  }
  
  if (is.null(strides)) strides <- rep(1L, length(filters))
  if (is.null(dilation_rate) & use.cnn) dilation_rate <- rep(1L, length(filters))
  
  if (use.cnn) {
    same_length <- (length(kernel_size) == length(filters)) &
      (length(filters) == length(strides)) &
      (length(strides) == length(dilation_rate))
    if (!same_length) {
      stop("kernel_size, filters, dilation_rate and strides must have the same length")
    }
  }
  
  stopifnot(maxlen > 0)
  stopifnot(dropout_lstm <= 1 & dropout_lstm >= 0)
  stopifnot(recurrent_dropout_lstm <= 1 & recurrent_dropout_lstm >= 0)
  
  if (length(layer_lstm) == 1) {
    layer_lstm <- rep(layer_lstm, layers.lstm)
  }
  
  if (stateful) {
    input_tensor <- keras::layer_input(batch_shape = c(batch_size, maxlen, vocabulary_size))
  } else {
    input_tensor <- keras::layer_input(shape = c(samples_per_target, maxlen, vocabulary_size))
  }
  
  if (use.cnn) {
    for (i in 1:length(filters)) {
      if (i == 1) {
        output_tensor <- input_tensor %>%
          keras::time_distributed(keras::layer_conv_1d(
            kernel_size = kernel_size[i],
            padding = padding,
            activation = "relu",
            filters = filters[i],
            strides = strides[i],
            dilation_rate = dilation_rate[i],
            input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
            use_bias = use_bias
          ))
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_max_pooling_1d(pool_size = pool_size[i]))
        }
        output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_batch_normalization(momentum = batch_norm_momentum))
      } else {
        
        output_tensor <- output_tensor %>%
          keras::time_distributed(keras::layer_conv_1d(
            kernel_size = kernel_size[i],
            padding = padding,
            activation = "relu",
            filters = filters[i],
            strides = strides[i],
            dilation_rate = dilation_rate[i],
            input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
            use_bias = use_bias
          ))
        output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_batch_normalization(momentum = batch_norm_momentum))
        
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_max_pooling_1d(pool_size = pool_size[i]))
        }
        #output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
      }
    }
  } else {
    if (zero_mask) {
      output_tensor <- input_tensor %>% keras::time_distributed(keras::layer_masking())
    } else {
      output_tensor <- input_tensor
    }
  }
  # lstm layers
  if (layers.lstm > 0) {
    if (layers.lstm > 1) {
      if (bidirectional) {
        for (i in 1:(layers.lstm - 1)) {
          output_tensor <- output_tensor %>%
            keras::time_distributed(keras::bidirectional(
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              keras::layer_lstm(
                units = layer_lstm[i],
                return_sequences = TRUE,
                dropout = dropout_lstm,
                recurrent_dropout = recurrent_dropout_lstm,
                stateful = stateful,
                recurrent_activation = "sigmoid"
              )
            ))
        }
      } else {
        for (i in 1:(layers.lstm - 1)) {
          output_tensor <- output_tensor %>%
            keras::time_distributed(keras::layer_lstm(
              units = layer_lstm[i],
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              return_sequences = TRUE,
              dropout = dropout_lstm,
              recurrent_dropout = recurrent_dropout_lstm,
              stateful = stateful,
              recurrent_activation = "sigmoid"
            ))
        }
      }
    }
    # last LSTM layer
    if (bidirectional) {
      output_tensor <- output_tensor %>%
        keras::time_distributed(keras::bidirectional(
          input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout_lstm, recurrent_dropout = recurrent_dropout_lstm,
                            stateful = stateful, recurrent_activation = "sigmoid")
        ))
    } else {
      output_tensor <- output_tensor %>%
        keras::time_distributed(keras::layer_lstm(units = layer_lstm[length(layer_lstm)],
                                                  input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
                                                  dropout = dropout_lstm, recurrent_dropout = recurrent_dropout_lstm, stateful = stateful,
                                                  recurrent_activation = "sigmoid"))
    }
  }
  
  if (!is.null(gap_time_dist)) {
    if (layers.lstm != 0) {
      stop("Global average pooling not compatible with using LSTM layer")
    }
    output_tensor <- pooling_flatten_time_dist(gap_time_dist, output_tensor)
  } else {
    if (layers.lstm == 0) {
      output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_flatten())
    }
  }
  
  num_aggr_layers <- 0
  aggr_layer_list <- list()
  
  if (!is.null(transformer_args)) {
    num_aggr_layers <- num_aggr_layers + 1
    for (i in seq_along(transformer_args$num_heads)) {
      attn_block <- layer_transformer_block_wrapper(
        num_heads = as.integer(transformer_args$num_heads[i]),
        head_size = as.integer(transformer_args$head_size[i]),
        dropout_rate = transformer_args$dropout[i],
        ff_dim = as.integer(transformer_args$ff_dim[i]),
        embed_dim = as.integer(output_tensor$shape[[3]]),
        vocabulary_size = as.integer(output_tensor$shape[[2]]),
        load_r6 = FALSE)
      if (i == 1) {
        output_tensor_attn <- output_tensor %>% attn_block
      } else {
        output_tensor_attn <- output_tensor_attn %>% attn_block
      }
    }
    output_tensor_attn <- pooling_flatten(transformer_args$pool_flatten, output_tensor_attn)
    aggr_layer_list[[num_aggr_layers]] <- output_tensor_attn
  }
  
  if (!is.null(aggregation_method)) {
    num_aggr_layers <- num_aggr_layers + 1
    layer_aggregate_td <- layer_aggregate_time_dist_wrapper(method = aggregation_method)
    output_tensor_aggregation_sum <- output_tensor %>% layer_aggregate_td
    aggr_layer_list[[num_aggr_layers]] <- output_tensor_aggregation_sum
  }
  
  if (!is.null(lstm_time_dist)) {
    num_aggr_layers <- num_aggr_layers + 1
    return_sequences <- TRUE
    for (i in 1:length(lstm_time_dist)) {
      if (i == length(lstm_time_dist)) {
        return_sequences <- FALSE 
      } 
      if (i == 1) {
        output_tensor_lstm <- output_tensor %>% keras::layer_lstm(units=lstm_time_dist[i], return_sequences = return_sequences)
      } else {
        output_tensor_lstm <- output_tensor_lstm %>% keras::layer_lstm(units=lstm_time_dist[i], return_sequences = return_sequences)
      }
    }
    aggr_layer_list[[num_aggr_layers]] <- output_tensor_lstm
  }
  
  if (num_aggr_layers == 0) {
    stop("You need to choose an aggregation method, either with aggregation_method, transformer_args or lstm_time_dist.")
  }
  
  if (num_aggr_layers == 1) {
    output_tensor <- aggr_layer_list[[1]]
  }
  
  if (num_aggr_layers > 1) {
    output_tensor <- keras::layer_concatenate(aggr_layer_list) 
  }
  
  if (length(layer_dense) > 1) {
    for (i in 1:(length(layer_dense) - 1)) {
      output_tensor <- output_tensor %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
    }
  }
  
  if (num_output_layers == 1) {
    output_tensor <- output_tensor %>%
      keras::layer_dense(units = num_targets, activation = last_layer_activation, dtype = "float32")
  } else {
    output_list <- list()
    for (i in 1:num_output_layers) {
      layer_name <- paste0("output_", i, "_", num_output_layers)
      output_list[[i]] <- output_tensor %>%
        keras::layer_dense(units = num_targets, activation = last_layer_activation, name = layer_name, dtype = "float32")
    }
  }
  
  if (num_output_layers == 1) {
    model <- keras::keras_model(inputs = input_tensor, outputs = output_tensor)
  } else {
    model <- keras::keras_model(inputs = input_tensor, outputs = output_list)
  }
  
  if (compile) {
    model <- compile_model(model = model, label_smoothing = label_smoothing, layer_dense = layer_dense,
                           solver = solver, learning_rate = learning_rate, loss_fn = loss_fn, 
                           num_output_layers = num_output_layers, label_noise_matrix = label_noise_matrix,
                           bal_acc = bal_acc, f1_metric = f1_metric, auc_metric = auc_metric)
  }
  
  argg <- as.list(environment())
  model <- add_hparam_list(model, argg)
  
  if (verbose) model$summary()
  return(model)
}

#' @title Create LSTM/CNN network that can process multiple samples for one target
#'
#' @description Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers with multiple 
#' input layers. After LSTM/CNN part all representations get aggregated by summation. 
#' Can be used to make single prediction for combination of multiple input sequences. 
#' Implements approach as described [here](https://arxiv.org/abs/1703.06114)
#'
#' @inheritParams create_model_lstm_cnn
#' @inheritParams create_model_lstm_cnn_time_dist
#' @param samples_per_target Number of samples to combine for one target.
#' @param dropout_dense Vector of dropout rates between dense layers. No dropout if `NULL`.
#' @param gap_inputs Global pooling method to apply. Same options as for `flatten_method` argument
#' in \link{create_model_transformer} function.
#' @examplesIf reticulate::py_module_available("tensorflow")
#' 
#' # Examples needs keras attached to run 
#' maxlen <- 50
#' \donttest{
#' library(keras)
#' create_model_lstm_cnn_multi_input(
#'   maxlen = maxlen,
#'   vocabulary_size = 4,
#'   samples_per_target = 7,
#'   kernel_size = c(10, 10),
#'   filters = c(64, 128),
#'   pool_size = c(2, 2),
#'   layer_lstm = c(32),
#'   layer_dense = c(64, 2),
#'   aggregation_method = c("max"),
#'   learning_rate = 0.001)
#'  }  
#'  
#' @returns A keras model with multiple input layers. Input goes through shared LSTM/CNN layers.   
#' @export
create_model_lstm_cnn_multi_input <- function(
    maxlen = 50,
    dropout_lstm = 0,
    recurrent_dropout_lstm = 0,
    layer_lstm = NULL,
    layer_dense = c(4),
    dropout_dense = NULL,
    solver = "adam",
    learning_rate = 0.001,
    vocabulary_size = 4,
    bidirectional = FALSE,
    batch_size = NULL,
    compile = TRUE,
    kernel_size = NULL,
    filters = NULL,
    strides = NULL,
    pool_size = NULL,
    padding = "same",
    dilation_rate = NULL,
    gap_inputs = NULL,
    use_bias = TRUE,
    zero_mask = FALSE,
    label_smoothing = 0,
    label_noise_matrix = NULL,
    last_layer_activation = "softmax",
    loss_fn = "categorical_crossentropy",
    auc_metric = FALSE,
    f1_metric = FALSE,
    bal_acc = FALSE,
    samples_per_target,
    batch_norm_momentum = 0.99,
    aggregation_method = c('sum'), 
    verbose = TRUE,
    model_seed = NULL,
    mixed_precision = FALSE,
    mirrored_strategy = NULL) {
  
  if (mixed_precision) tensorflow::tf$keras$mixed_precision$set_global_policy("mixed_float16")
  
  if (is.null(mirrored_strategy)) mirrored_strategy <- ifelse(count_gpu() > 1, TRUE, FALSE)
  if (mirrored_strategy) {
    mirrored_strategy <- tensorflow::tf$distribute$MirroredStrategy()
    with(mirrored_strategy$scope(), { 
      argg <- as.list(environment())
      argg$mirrored_strategy <- FALSE
      model <- do.call(create_model_lstm_cnn_multi_input, argg)
    })
    return(model)
  }
  
  layer_dense <- as.integer(layer_dense)
  if (!is.null(dropout_dense)) stopifnot(length(dropout_dense) == length(layer_dense))
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  num_targets <- layer_dense[length(layer_dense)]
  layers.lstm <- length(layer_lstm)
  use.cnn <- ifelse(!is.null(kernel_size), TRUE, FALSE)
  num_output_layers = 1
  
  if (!is.null(layer_lstm)) {
    stopifnot(length(layer_lstm) == 1 | (length(layer_lstm) ==  layers.lstm))
  }
  
  if (layers.lstm == 0 & !use.cnn) {
    stop("Model does not use LSTM or CNN layers.")
  }
  
  if (is.null(strides)) strides <- rep(1L, length(filters))
  if (is.null(dilation_rate) & use.cnn) dilation_rate <- rep(1L, length(filters))
  
  if (use.cnn) {
    same_length <- (length(kernel_size) == length(filters)) &
      (length(filters) == length(strides)) &
      (length(strides) == length(dilation_rate))
    if (!same_length) {
      stop("kernel_size, filters, dilation_rate and strides must have the same length")
    }
  }
  
  stopifnot(maxlen > 0)
  stopifnot(dropout_lstm <= 1 & dropout_lstm >= 0)
  stopifnot(recurrent_dropout_lstm <= 1 & recurrent_dropout_lstm >= 0)
  
  if (length(layer_lstm) == 1) {
    layer_lstm <- rep(layer_lstm, layers.lstm)
  }
  
  input_tensor <- keras::layer_input(shape = c(maxlen, vocabulary_size))
  
  if (use.cnn) {
    for (i in 1:length(filters)) {
      if (i == 1) {
        output_tensor <- input_tensor %>%
          keras::layer_conv_1d(
            kernel_size = kernel_size[i],
            padding = padding,
            activation = "relu",
            filters = filters[i],
            strides = strides[i],
            dilation_rate = dilation_rate[i],
            input_shape = c(maxlen, vocabulary_size),
            use_bias = use_bias
          )
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
      } else {
        
        output_tensor <- output_tensor %>%
          keras::layer_conv_1d(
            kernel_size = kernel_size[i],
            padding = padding,
            activation = "relu",
            filters = filters[i],
            strides = strides[i],
            dilation_rate = dilation_rate[i],
            input_shape = c(maxlen, vocabulary_size),
            use_bias = use_bias
          )
        output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
        
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        
      }
    }
  } else {
    if (zero_mask) {
      output_tensor <- input_tensor %>% keras::layer_masking()
    } else {
      output_tensor <- input_tensor
    }
  }
  # lstm layers
  if (layers.lstm > 0) {
    if (layers.lstm > 1) {
      if (bidirectional) {
        for (i in 1:(layers.lstm - 1)) {
          output_tensor <- output_tensor %>%
            keras::bidirectional(
              input_shape = c(maxlen, vocabulary_size),
              keras::layer_lstm(
                units = layer_lstm[i],
                return_sequences = TRUE,
                dropout = dropout_lstm,
                recurrent_dropout = recurrent_dropout_lstm,
                recurrent_activation = "sigmoid"
              )
            )
        }
      } else {
        for (i in 1:(layers.lstm - 1)) {
          output_tensor <- output_tensor %>%
            keras::layer_lstm(
              units = layer_lstm[i],
              input_shape = c(maxlen, vocabulary_size),
              return_sequences = TRUE,
              dropout = dropout_lstm,
              recurrent_dropout = recurrent_dropout_lstm,
              recurrent_activation = "sigmoid"
            )
        }
      }
    }
    # last LSTM layer
    if (bidirectional) {
      output_tensor <- output_tensor %>%
        keras::bidirectional(
          input_shape = c(maxlen, vocabulary_size),
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout_lstm, recurrent_dropout = recurrent_dropout_lstm,
                            recurrent_activation = "sigmoid")
        )
    } else {
      output_tensor <- output_tensor %>%
        keras::layer_lstm(units = layer_lstm[length(layer_lstm)],
                          input_shape = c(maxlen, vocabulary_size),
                          dropout = dropout_lstm, recurrent_dropout = recurrent_dropout_lstm,
                          recurrent_activation = "sigmoid")
    }
  }
  
  if (!is.null(gap_inputs)) {
    if (layers.lstm != 0) {
      stop("Global average pooling not compatible with using LSTM layer")
    }
    output_tensor <- output_tensor %>% pooling_flatten(global_pooling = gap_inputs)
  } else {
    if (layers.lstm == 0) {
      output_tensor <- output_tensor %>% keras::layer_flatten()
    }
  }
  
  feature_ext_model <- keras::keras_model(inputs = input_tensor, outputs = output_tensor)
  
  input_list <- list()
  representation_list <- list()
  for (i in 1:samples_per_target) {
    input_list[[i]] <- keras::layer_input(shape = c(maxlen, vocabulary_size), name = paste0("input_", i))
    representation_list[[i]] <- feature_ext_model(input_list[[i]])
  }
  
  if (!is.null(aggregation_method)) {
    layer_aggregate_td <- layer_aggregate_time_dist_wrapper(method = aggregation_method, multi_in = TRUE)
    y <- representation_list %>% layer_aggregate_td
  }
  
  if (length(layer_dense) > 1) {
    for (i in 1:(length(layer_dense) - 1)) {
      if (!is.null(dropout_dense)) y <- y %>% keras::layer_dropout(dropout_dense[i])
      y <- y %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
    }
  }
  
  y <- y %>% keras::layer_dense(units = num_targets, activation = last_layer_activation, dtype = "float32")
  model <- keras::keras_model(inputs = input_list, outputs = y)
  
  if (compile) {
    model <- compile_model(model = model, label_smoothing = label_smoothing, layer_dense = layer_dense,
                           solver = solver, learning_rate = learning_rate, loss_fn = loss_fn, 
                           num_output_layers = num_output_layers, label_noise_matrix = label_noise_matrix,
                           bal_acc = bal_acc, f1_metric = f1_metric, auc_metric = auc_metric)
  }
  
  argg <- c(as.list(environment()))
  model <- add_hparam_list(model, argg)
  
  if (verbose) model$summary()
  return(model)
}
