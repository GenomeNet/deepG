#' @title Create LSTM/CNN network
#'
#' @description Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers.
#' Last layer is a dense layer.
#' 
#' @param maxlen Length of predictor sequence.
#' @param dropout_lstm Fraction of the units to drop for inputs.
#' @param recurrent_dropout_lstm Fraction of the units to drop for recurrent state.
#' @param layer_lstm Number of cells per network layer. Can be a scalar or vector.
#' @param layer_dense Vector specifying number of neurons per dense layer after last LSTM or CNN layer (if no LSTM used).
#' @param dropout_dense Dropout rates between dense layers. No dropout if `ǸULL`.
#' @param solver Optimization method, options are `"adam", "adagrad", "rmsprop"` or `"sgd"`.
#' @param learning_rate Learning rate for optimizer.
#' @param bidirectional Use bidirectional wrapper for lstm layers.
#' @param vocabulary_size Number of unique character in vocabulary.
#' @param stateful Boolean. Whether to use stateful LSTM layer.
#' @param batch_size Number of samples that are used for one network update. Only used if \code{stateful = TRUE}.
#' @param compile Whether to compile the model.
#' @param kernel_size Size of 1d convolutional layers. For multiple layers, assign a vector. (e.g, `rep(3,2)` for two layers and kernel size 3)
#' @param filters Number of filters. For multiple layers, assign a vector.
#' @param strides Stride values. For multiple layers, assign a vector.
#' @param pool_size Integer, size of the max pooling windows. For multiple layers, assign a vector.
#' @param padding Padding of CNN layers, e.g. `"same", "valid"` or `"causal"`.
#' @param dilation_rate Integer, the dilation rate to use for dilated convolution.
#' @param gap Whether to apply global average pooling after last CNN layer.
#' @param use_bias Boolean. Usage of bias for CNN layers.
#' @param residual_block Boolean. If true, the residual connections are used in CNN. It is not used in the first convolutional layer.
#' @param residual_block_length Integer. Determines how many convolutional layers (or triplets when `size_reduction_1D_conv` is `TRUE`) exist
#  between the legs of a residual connection. e.g. if the `length kernel_size/filters` is 7 and `residual_block_length` is 2, there are 1+(7-1)*2 convolutional
#  layers in the model when `size_reduction_1Dconv` is FALSE and 1+(7-1)*2*3 convolutional layers when `size_reduction_1Dconv` is TRUE.
#' @param size_reduction_1Dconv Boolean. When `TRUE`, the number of filters in the convolutional layers is reduced to 1/4 of the number of filters of
#  the original layer by a convolution layer with kernel size 1, and number of filters are increased back to the original value by a convolution layer
#  with kernel size 1 after the convolution with original kernel size with reduced number of filters.
#' @param label_input Integer or `NULL`. If not `NULL`, adds additional input layer of \code{label_input} size.
#' @param zero_mask Boolean, whether to apply zero masking before LSTM layer. Only used if model does not use any CNN layers.
#' @param label_smoothing Float in \[0, 1\]. If 0, no smoothing is applied. If > 0, loss between the predicted
#' labels and a smoothed version of the true labels, where the smoothing squeezes the labels towards 0.5.
#' The closer the argument is to 1 the more the labels get smoothed.
#' @param label_noise_matrix Matrix of label noises. Every row stands for one class and columns for percentage of labels in that class.
#' If first label contains 5 percent wrong labels and second label no noise, then
#' 
#' \code{label_noise_matrix <- matrix(c(0.95, 0.05, 0, 1), nrow = 2, byrow = TRUE )}
#' @param last_layer_activation Either `"sigmoid"` or `"softmax"`.
#' @param loss_fn Either `"categorical_crossentropy"` or `"binary_crossentropy"`. If `label_noise_matrix` given, will use custom `"noisy_loss"`.
#' @param num_output_layers Number of output layers.
#' @param auc_metric Whether to add AUC metric.
#' @param f1_metric Whether to add F1 metric.
#' @param bal_acc Whether to add balanced accuracy.
#' @param verbose Boolean.
#' @param batch_norm_momentum Momentum for the moving mean and the moving variance.
#' @param model_seed Set seed for model parameters in tensorflow if not `NULL`.
#' @param mixed_precision Whether to use mixed precision (https://www.tensorflow.org/guide/mixed_precision).
#' @param mirrored_strategy Whether to use distributed mirrored strategy. If NULL, will use distributed mirrored strategy only if >1 GPU available.   
#' @examples 
#' create_model_lstm_cnn(
#'   maxlen = 500,
#'   vocabulary_size = 4,
#'   kernel_size = c(8, 8, 8),
#'   filters = c(16, 32, 64),
#'   pool_size = c(3, 3, 3),
#'   layer_lstm = c(32, 64),
#'   layer_dense = c(128, 4),
#'   learning_rate = 0.001)
#' @export
create_model_lstm_cnn <- function(
    maxlen = 50,
    dropout_lstm = 0,
    recurrent_dropout_lstm = 0,
    layer_lstm = NULL,
    layer_dense = c(4),
    dropout_dense = NULL,
    kernel_size = NULL,
    filters = NULL,
    strides = NULL,
    pool_size = NULL,
    solver = "adam",
    learning_rate = 0.001,
    vocabulary_size = 4,
    bidirectional = FALSE,
    stateful = FALSE,
    batch_size = NULL,
    compile = TRUE,
    padding = "same",
    dilation_rate = NULL,
    gap = FALSE,
    use_bias = TRUE,
    residual_block = FALSE,
    residual_block_length = 1,
    size_reduction_1Dconv = FALSE,
    label_input = NULL,
    zero_mask = FALSE,
    label_smoothing = 0,
    label_noise_matrix = NULL,
    last_layer_activation = "softmax",
    loss_fn = "categorical_crossentropy",
    num_output_layers = 1,
    auc_metric = FALSE,
    f1_metric = FALSE,
    bal_acc = TRUE,
    verbose = TRUE,
    batch_norm_momentum = 0.99,
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
      model <- do.call(create_model_lstm_cnn, argg)
    })
    return(model)
  }
  
  layer_dense <- as.integer(layer_dense)
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
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
    if (residual_block & (padding != "same")) {
      stop("Padding option must be same when residual block is used.")
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
    input_tensor <- keras::layer_input(shape = c(maxlen, vocabulary_size))
  }
  
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
            input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
            use_bias = use_bias
          )
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
      } else {
        if (residual_block){
          if ((strides[i] > 1) | (pool_size[i] > 1)) {
            residual_layer <- output_tensor %>% keras::layer_average_pooling_1d(pool_size=strides[i]*pool_size[i])
          } else {
            residual_layer <- output_tensor
          }
          if (filters[i-1] != filters[i]){
            residual_layer <- residual_layer %>%
              keras::layer_conv_1d(
                kernel_size = 1,
                padding = padding,
                activation = "relu",
                filters = filters[i],
                strides = 1,
                dilation_rate = dilation_rate[i],
                input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
                use_bias = use_bias
              )
            residual_layer <- residual_layer %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
          }
          if (residual_block_length > 1){
            for (j in 1:(residual_block_length-1)){
              if (size_reduction_1Dconv){
                output_tensor <- output_tensor %>%
                  keras::layer_conv_1d(
                    kernel_size = 1,
                    padding = padding,
                    activation = "relu",
                    filters = filters[i]/4,
                    strides = 1,
                    dilation_rate = dilation_rate[i],
                    input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
                    use_bias = use_bias
                  )
                output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
                
                output_tensor <- output_tensor %>%
                  keras::layer_conv_1d(
                    kernel_size = kernel_size[i],
                    padding = padding,
                    activation = "relu",
                    filters = filters[i]/4,
                    strides = 1,
                    dilation_rate = dilation_rate[i],
                    input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
                    use_bias = use_bias
                  )
                output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
                
                output_tensor <- output_tensor %>%
                  keras::layer_conv_1d(
                    kernel_size = 1,
                    padding = padding,
                    activation = "relu",
                    filters = filters[i],
                    strides = 1,
                    dilation_rate = dilation_rate[i],
                    input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
                    use_bias = use_bias
                  )
                output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
                
              } else {
                output_tensor <- output_tensor %>%
                  keras::layer_conv_1d(
                    kernel_size = kernel_size[i],
                    padding = padding,
                    activation = "relu",
                    filters = filters[i],
                    strides = 1,
                    dilation_rate = dilation_rate[i],
                    input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
                    use_bias = use_bias
                  )
                output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
              }
            }
          }
        }
        if (size_reduction_1Dconv){
          output_tensor <- output_tensor %>%
            keras::layer_conv_1d(
              kernel_size = 1,
              padding = padding,
              activation = "relu",
              filters = filters[i]/4,
              strides = 1,
              dilation_rate = dilation_rate[i],
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              use_bias = use_bias
            )
          output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
          
          output_tensor <- output_tensor %>%
            keras::layer_conv_1d(
              kernel_size = kernel_size[i],
              padding = padding,
              activation = "relu",
              filters = filters[i]/4,
              strides = strides[i],
              dilation_rate = dilation_rate[i],
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              use_bias = use_bias
            )
          output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
          
          output_tensor <- output_tensor %>%
            keras::layer_conv_1d(
              kernel_size = 1,
              padding = padding,
              activation = "relu",
              filters = filters[i],
              strides = 1,
              dilation_rate = dilation_rate[i],
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              use_bias = use_bias
            )
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              use_bias = use_bias
            )
          output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
        }
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        #output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
        if (residual_block){
          output_tensor <- keras::layer_add(list(output_tensor, residual_layer))
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              keras::layer_lstm(
                units = layer_lstm[i],
                return_sequences = TRUE,
                dropout = dropout_lstm,
                recurrent_dropout = recurrent_dropout_lstm,
                stateful = stateful,
                recurrent_activation = "sigmoid"
              )
            )
        }
      } else {
        for (i in 1:(layers.lstm - 1)) {
          output_tensor <- output_tensor %>%
            keras::layer_lstm(
              layer_lstm[i],
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              return_sequences = TRUE,
              dropout = dropout_lstm,
              recurrent_dropout = recurrent_dropout_lstm,
              stateful = stateful,
              recurrent_activation = "sigmoid"
            )
        }
      }
    }
    # last LSTM layer
    if (bidirectional) {
      output_tensor <- output_tensor %>%
        keras::bidirectional(
          input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout_lstm, recurrent_dropout = recurrent_dropout_lstm,
                            stateful = stateful, recurrent_activation = "sigmoid")
        )
    } else {
      output_tensor <- output_tensor %>%
        keras::layer_lstm(units = layer_lstm[length(layer_lstm)],
                          input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
                          dropout = dropout_lstm, recurrent_dropout = recurrent_dropout_lstm, stateful = stateful,
                          recurrent_activation = "sigmoid")
    }
  }
  
  if (gap) {
    if (layers.lstm != 0) {
      stop("Global average pooling not compatible with using LSTM layer")
    }
    output_tensor <- output_tensor %>% keras::layer_global_average_pooling_1d()
  } else {
    if (layers.lstm == 0) {
      output_tensor <- output_tensor %>% keras::layer_flatten()
    }
  }
  
  if (!is.null(label_input)) {
    input_label_list <- list()
    for (i in 1:length(label_input)) {
      if (!stateful) {
        eval(parse(text = paste0("label_input_layer_", as.character(i), "<- keras::layer_input(c(label_input[i]))")))
      } else {
        eval(parse(text = paste0("label_input_layer_", as.character(i), "<- keras::layer_input(batch_shape = c(batch_size, label_input[i]))")))
      }
      input_label_list[[i]] <- eval(parse(text = paste0("label_input_layer_", as.character(i))))
    }
    output_tensor <- keras::layer_concatenate(c(
      input_label_list, output_tensor
    )
    )
  }
  
  if (length(layer_dense) > 1) {
    for (i in 1:(length(layer_dense) - 1)) {
      if (!is.null(dropout_dense)) output_tensor <- output_tensor %>% keras::layer_dropout(dropout_dense[i])
      output_tensor <- output_tensor %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
    }
  }
  
  if (num_output_layers == 1) {
    if (!is.null(dropout_dense)) output_tensor <- output_tensor %>% keras::layer_dropout(dropout_dense[length(dropout_dense)])
    output_tensor <- output_tensor %>%
      keras::layer_dense(units = num_targets, activation = last_layer_activation, dtype = "float32")
  } else {
    output_list <- list()
    for (i in 1:num_output_layers) {
      layer_name <- paste0("output_", i, "_", num_output_layers)
      if (!is.null(dropout_dense)) {
        output_list[[i]] <- output_tensor %>% keras::layer_dropout(dropout_dense[length(dropout_dense)])
        output_list[[i]] <- output_list[[i]] %>%
          keras::layer_dense(units = num_targets, activation = last_layer_activation, name = layer_name, dtype = "float32")
      } else {
        output_list[[i]] <- output_tensor %>%
          keras::layer_dense(units = num_targets, activation = last_layer_activation, name = layer_name, dtype = "float32")
      }
    }
  }
  
  if (!is.null(label_input)) {
    label_inputs <- list()
    for (i in 1:length(label_input)) {
      eval(parse(text = paste0("label_inputs$label_input_layer_", as.character(i), "<- label_input_layer_", as.character(i))))
    }
    if (num_output_layers == 1) {
      model <- keras::keras_model(inputs = list(label_inputs, input_tensor), outputs = output_tensor)
    } else {
      model <- keras::keras_model(inputs = list(label_inputs, input_tensor), outputs = output_list)
    }
  } else {
    if (num_output_layers == 1) {
      model <- keras::keras_model(inputs = input_tensor, outputs = output_tensor)
    } else {
      model <- keras::keras_model(inputs = input_tensor, outputs = output_list)
    }
  }
  
  optimizer <- set_optimizer(solver, learning_rate) 
  
  #add metrics
  if (loss_fn == "binary_crossentropy") {
    model_metrics <- c(tf$keras$metrics$BinaryAccuracy(name = "acc"))
  } else {
    model_metrics <- c("acc")
  } 
  
  cm_dir <- NULL
  if (num_output_layers == 1) {
    cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
    while (dir.exists(cm_dir)) {
      cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
    }
    dir.create(cm_dir)
    model$cm_dir <- cm_dir
    
    if (loss_fn == "categorical_crossentropy") {
      
      if (bal_acc) {
        macro_average_cb <- balanced_acc_wrapper(num_targets, cm_dir)
        model_metrics <- c(macro_average_cb, "acc")
      }
      
      if (f1_metric) {
        f1 <- f1_wrapper(num_targets)
        model_metrics <- c(model_metrics, f1)
      }
    }
    
    if (auc_metric) {
      auc <- auc_wrapper(model_output_size = layer_dense[length(layer_dense)],
                         loss = loss_fn)
      model_metrics <- c(model_metrics, auc)
    }
  }
  
  if (label_smoothing > 0 & !is.null(label_noise_matrix)) {
    stop("Can not use label smoothing and label noise at the same time. Either set label_smoothing = 0 or label_noise_matrix = NULL")
  }
  
  if (label_smoothing > 0) {
    if (loss_fn == "categorical_crossentropy") {
      smooth_loss <- tensorflow::tf$losses$CategoricalCrossentropy(label_smoothing = label_smoothing, name = "smooth_loss")
    }
    if (loss_fn == "binary_crossentropy") {
      smooth_loss <- tensorflow::tf$losses$BinaryCrossentropy(label_smoothing = label_smoothing, name = "smooth_loss")
    }
    model %>% keras::compile(loss = smooth_loss,
                             optimizer = optimizer, metrics = model_metrics)
  } else if (!is.null(label_noise_matrix)) {
    row_sums <- rowSums(label_noise_matrix)
    if (!all(row_sums == 1)) {
      warning("Sum of noise matrix rows don't add up to 1")
    }
    noisy_loss <- noisy_loss_wrapper(solve(label_noise_matrix))
    model %>% keras::compile(loss =  noisy_loss,
                             optimizer = optimizer, metrics = model_metrics)
  } else {
    model %>% keras::compile(loss = loss_fn,
                             optimizer = optimizer, metrics = model_metrics)
  }
  
  argg <- c(as.list(environment()))
  model <- add_hparam_list(model, argg)
  model$cm_dir <- cm_dir
  
  if (verbose) summary(model)
  return_model <- model
}

#' Create wavenet model
#'
#' Creates network architecture as described [here](https://arxiv.org/abs/1609.03499). Implementation 
#' uses code from [here](https://github.com/r-tensorflow/wavenet).
#' 
#' @inheritParams wavenet::wavenet
#' @inheritParams create_model_lstm_cnn
#' @examples 
#'\dontrun{
#' model <- create_model_wavenet(residual_blocks = 2^rep(1:4, 2), maxlen = 1000)
#' }
#' @export
create_model_wavenet <- function(filters = 16, kernel_size = 2, residual_blocks, maxlen,
                                 input_tensor = NULL, initial_kernel_size = 32, initial_filters = 32,
                                 output_channels = 4, output_activation = "softmax", solver = "adam",
                                 learning_rate = 0.001, compile = TRUE, verbose = TRUE, model_seed = NULL,
                                 mixed_precision = FALSE, mirrored_strategy = NULL) {
  
  if (mixed_precision) tensorflow::tf$keras$mixed_precision$set_global_policy("mixed_float16")
  
  if (is.null(mirrored_strategy)) mirrored_strategy <- ifelse(count_gpu() > 1, TRUE, FALSE)
  if (mirrored_strategy) {
    mirrored_strategy <- tensorflow::tf$distribute$MirroredStrategy()
    with(mirrored_strategy$scope(), { 
      argg <- as.list(environment())
      argg$mirrored_strategy <- FALSE
      model <- do.call(create_model_wavenet, argg)
    })
    return(model)
  }
  
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  
  model <- wavenet::wavenet(filters = filters, kernel_size = kernel_size, residual_blocks = residual_blocks,
                            input_shape = list(maxlen, output_channels), input_tensor = input_tensor, initial_kernel_size = initial_kernel_size,
                            initial_filters = initial_filters, output_channels = output_channels, output_activation = "softmax")
  
  optimizer <- set_optimizer(solver, learning_rate)
  
  if (compile) {
    model %>% keras::compile(loss = "categorical_crossentropy",
                             optimizer = optimizer, metrics = c("acc"))
  }
  
  argg <- c(as.list(environment()))
  model <- add_hparam_list(model, argg)
  
  reticulate::py_set_attr(x = model, name = "hparam", value = model$hparam)
  if (verbose) summary(model)
  model
}

#' @title Create LSTM/CNN network to predict middle part of a sequence
#'
#' @description
#' Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers.
#' Function creates two sub networks consisting each of (optional) CNN layers followed by an arbitrary number of LSTM layers. Afterwards the last LSTM layers
#' get concatenated and followed by one or more dense layers. Last layer is a dense layer.
#' Network tries to predict target in the middle of a sequence. If input is AACCTAAGG, input tensors should correspond to x1 = AACC, x2 = GGAA and y = T.
#' @inheritParams create_model_lstm_cnn
#' @examples
#' create_model_lstm_cnn_target_middle(
#'   maxlen = 500,
#'   vocabulary_size = 4,
#'   kernel_size = c(8, 8, 8),
#'   filters = c(16, 32, 64),
#'   pool_size = c(3, 3, 3),
#'   layer_lstm = c(32, 64),
#'   layer_dense = c(128, 4),
#'   learning_rate = 0.001)
#' @export
create_model_lstm_cnn_target_middle <- function(
    maxlen = 50,
    dropout_lstm = 0,
    recurrent_dropout_lstm = 0,
    layer_lstm = 128,
    solver = "adam",
    learning_rate = 0.001,
    vocabulary_size = 4,
    bidirectional = FALSE,
    stateful = FALSE,
    batch_size = NULL,
    padding = "same",
    compile = TRUE,
    layer_dense = NULL,
    kernel_size = NULL,
    filters = NULL,
    pool_size = NULL,
    strides = NULL,
    label_input = NULL,
    zero_mask = FALSE,
    label_smoothing = 0,
    label_noise_matrix = NULL,
    last_layer_activation = "softmax",
    loss_fn = "categorical_crossentropy",
    num_output_layers = 1,
    f1_metric = FALSE,
    verbose = TRUE,
    batch_norm_momentum = 0.99,
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
      model <- do.call(create_model_lstm_cnn_target_middle, argg)
    })
    return(model)
  }
  
  layer_dense <- as.integer(layer_dense)
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  use.cnn <- ifelse(!is.null(kernel_size), TRUE, FALSE)
  num_targets <- layer_dense[length(layer_dense)]
  layers.lstm <- length(layer_lstm)
  
  stopifnot(length(layer_lstm) == 1 | (length(layer_lstm) ==  layers.lstm))
  stopifnot(maxlen > 0)
  stopifnot(dropout_lstm <= 1 & dropout_lstm >= 0)
  stopifnot(recurrent_dropout_lstm <= 1 & recurrent_dropout_lstm >= 0)
  
  if (!is.null(layer_lstm)) {
    stopifnot(length(layer_lstm) == 1 | (length(layer_lstm) ==  layers.lstm))
  }
  
  if (is.null(strides)) {
    strides <- rep(1L, length(filters))
  }
  
  if (use.cnn) {
    same_length <- (length(kernel_size) == length(filters)) & (length(filters) == length(strides))
    if (!same_length) {
      stop("kernel_size, filters and strides must have the same length")
    }
  }
  
  # length of split sequences
  maxlen_1 <- ceiling(maxlen/2)
  maxlen_2 <- floor(maxlen/2)
  if (stateful) {
    input_tensor_1 <- keras::layer_input(batch_shape = c(batch_size, maxlen_1, vocabulary_size))
  } else {
    input_tensor_1 <- keras::layer_input(shape = c(maxlen_1, vocabulary_size))
  }
  
  if (use.cnn) {
    for (i in 1:length(filters)) {
      if (i == 1) {
        output_tensor_1 <- input_tensor_1 %>%
          keras::layer_conv_1d(
            kernel_size = kernel_size[i],
            padding = padding,
            activation = "relu",
            filters = filters[i],
            strides = strides[i],
            input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL)
          )
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor_1 <- output_tensor_1 %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        output_tensor_1 <- output_tensor_1 %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
      } else {
        output_tensor_1 <- output_tensor_1 %>%
          keras::layer_conv_1d(
            kernel_size = kernel_size[i],
            padding = padding,
            activation = "relu",
            strides = strides[i],
            filters = filters[i]
          )
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor_1 <- output_tensor_1 %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        output_tensor_1 <- output_tensor_1 %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
      }
    }
  } else {
    if (zero_mask) {
      output_tensor_1 <- input_tensor_1 %>% keras::layer_masking()
    } else {
      output_tensor_1 <- input_tensor_1
    }
  }
  
  # lstm layers
  if (!is.null(layers.lstm) && layers.lstm > 0) {
    if (layers.lstm > 1) {
      if (bidirectional) {
        for (i in 1:(layers.lstm - 1)) {
          output_tensor_1 <- output_tensor_1 %>%
            keras::bidirectional(
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              keras::layer_lstm(
                units = layer_lstm[i],
                return_sequences = TRUE,
                dropout = dropout_lstm,
                recurrent_dropout = recurrent_dropout_lstm,
                stateful = stateful,
                recurrent_activation = "sigmoid"
              )
            )
        }
      } else {
        for (i in 1:(layers.lstm - 1)) {
          output_tensor_1 <- output_tensor_1 %>%
            keras::layer_lstm(
              units = layer_lstm[i],
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              return_sequences = TRUE,
              dropout = dropout_lstm,
              recurrent_dropout = recurrent_dropout_lstm,
              stateful = stateful,
              recurrent_activation = "sigmoid"
            )
        }
      }
    }
    
    # last LSTM layer
    if (bidirectional) {
      output_tensor_1 <- output_tensor_1 %>%
        keras::bidirectional(
          input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout_lstm, recurrent_dropout = recurrent_dropout_lstm,
                            stateful = stateful, recurrent_activation = "sigmoid")
        )
    } else {
      output_tensor_1 <- output_tensor_1 %>%
        keras::layer_lstm(units = layer_lstm[length(layer_lstm)],
                          input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
                          dropout = dropout_lstm, recurrent_dropout = recurrent_dropout_lstm, stateful = stateful,
                          recurrent_activation = "sigmoid")
    }
  }
  
  if (stateful) {
    input_tensor_2 <- keras::layer_input(batch_shape = c(batch_size, maxlen_2, vocabulary_size))
  } else {
    input_tensor_2 <- keras::layer_input(shape = c(maxlen_2, vocabulary_size))
  }
  
  if (use.cnn) {
    for (i in 1:length(filters)) {
      if (i == 1) {
        output_tensor_2 <- input_tensor_2 %>%
          keras::layer_conv_1d(
            kernel_size = kernel_size[i],
            padding = padding,
            activation = "relu",
            filters = filters[i],
            strides = strides[i],
            input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL)
          )
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor_2 <- output_tensor_2 %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        output_tensor_2 <- output_tensor_2 %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
      } else {
        output_tensor_2 <- output_tensor_2 %>%
          keras::layer_conv_1d(
            kernel_size = kernel_size[i],
            padding = padding,
            activation = "relu",
            strides = strides[i],
            filters = filters[i]
          )
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor_2 <- output_tensor_2 %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        output_tensor_2 <- output_tensor_2 %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
      }
    }
  } else {
    if (zero_mask) {
      output_tensor_2 <- input_tensor_2 %>% keras::layer_masking()
    } else {
      output_tensor_2 <- input_tensor_2
    }
  }
  
  
  # lstm layers
  if (!is.null(layers.lstm) && layers.lstm > 0) {
    if (layers.lstm > 1) {
      if (bidirectional) {
        for (i in 1:(layers.lstm - 1)) {
          output_tensor_2 <- output_tensor_2 %>%
            keras::bidirectional(
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              keras::layer_lstm(
                units = layer_lstm[i],
                return_sequences = TRUE,
                dropout = dropout_lstm,
                recurrent_dropout = recurrent_dropout_lstm,
                stateful = stateful,
                recurrent_activation = "sigmoid"
              )
            )
        }
      } else {
        for (i in 1:(layers.lstm - 1)) {
          output_tensor_2 <- output_tensor_2 %>%
            keras::layer_lstm(
              units = layer_lstm[i],
              input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
              return_sequences = TRUE,
              dropout = dropout_lstm,
              recurrent_dropout = recurrent_dropout_lstm,
              stateful = stateful,
              recurrent_activation = "sigmoid"
            )
        }
      }
    }
    
    # last LSTM layer
    if (bidirectional) {
      output_tensor_2 <- output_tensor_2 %>%
        keras::bidirectional(
          input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout_lstm, recurrent_dropout = recurrent_dropout_lstm,
                            stateful = stateful, recurrent_activation = "sigmoid")
        )
    } else {
      output_tensor_2 <- output_tensor_2 %>%
        keras::layer_lstm(units = layer_lstm[length(layer_lstm)],
                          input_shape = switch(stateful + 1, c(maxlen, vocabulary_size), NULL),
                          dropout = dropout_lstm, recurrent_dropout = recurrent_dropout_lstm, stateful = stateful,
                          recurrent_activation = "sigmoid")
    }
  }
  
  output_tensor <- keras::layer_concatenate(list(output_tensor_1, output_tensor_2))
  
  if (layers.lstm == 0) {
    output_tensor <- output_tensor %>% keras::layer_flatten()
  }
  
  if (!is.null(label_input)) {
    input_label_list <- list()
    for (i in 1:length(label_input)) {
      if (!stateful) {
        eval(parse(text = paste0("label_input_layer_", as.character(i), "<- keras::layer_input(c(label_input[i]))")))
      } else {
        eval(parse(text = paste0("label_input_layer_", as.character(i), "<- keras::layer_input(batch_shape = c(batch_size, label_input[i]))")))
      }
      input_label_list[[i]] <- eval(parse(text = paste0("label_input_layer_", as.character(i))))
    }
    output_tensor <- keras::layer_concatenate(c(
      input_label_list, output_tensor
    )
    )
  }
  
  if (length(layer_dense) > 1) {
    for (i in 1:(length(layer_dense) - 1)) {
      output_tensor <- output_tensor %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
    }
  }
  
  if (num_output_layers == 1) {
    output_tensor <- output_tensor %>%
      keras::layer_dense(units = num_targets, activation = last_layer_activation, dtype = "float32")
  }  else {
    output_list <- list()
    for (i in 1:num_output_layers) {
      layer_name <- paste0("output_", i, "_", num_output_layers)
      output_list[[i]] <- output_tensor %>%
        keras::layer_dense(units = num_targets, activation = last_layer_activation, name = layer_name, dtype = "float32")
    }
  }
  
  # print model layout to screen
  if (!is.null(label_input)) {
    label_inputs <- list()
    for (i in 1:length(label_input)) {
      eval(parse(text = paste0("label_inputs$label_input_layer_", as.character(i), "<- label_input_layer_", as.character(i))))
    }
    model <- keras::keras_model(inputs = c(label_inputs, input_tensor_1, input_tensor_2), outputs = output_tensor)
  } else {
    model <- keras::keras_model(inputs = list(input_tensor_1, input_tensor_2), outputs = output_tensor)
  }
  
  # choose optimization method
  optimizer <- set_optimizer(solver, learning_rate) 
  
  #add metrics
  cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  while (dir.exists(cm_dir)) {
    cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  }
  dir.create(cm_dir)
  model$cm_dir <- cm_dir
  
  #add metrics
  if (loss_fn == "binary_crossentropy") {
    model_metrics <- c(tf$keras$metrics$BinaryAccuracy(name = "acc"))
  } else {
    model_metrics <- c("acc")
  } 
  
  if (loss_fn == "categorical_crossentropy") {
    
    macro_average_cb <- balanced_acc_wrapper(num_targets, cm_dir)
    model_metrics <- c(macro_average_cb, "acc")
    
    if (f1_metric) {
      f1 <- f1_wrapper(num_targets)
      model_metrics <- c(model_metrics, f1)
    }
  }
  
  if (label_smoothing > 0 & !is.null(label_noise_matrix)) {
    stop("Can not use label smoothing and label noise at the same time. Either set label_smoothing = 0 or label_noise_matrix = NULL")
  }
  if (label_smoothing > 0) {
    if (loss_fn == "categorical_crossentropy") {
      smooth_loss <- tensorflow::tf$losses$CategoricalCrossentropy(label_smoothing = label_smoothing, name = "smooth_loss")
    }
    if (loss_fn == "binary_crossentropy") {
      smooth_loss <- tensorflow::tf$losses$BinaryCrossentropy(label_smoothing = label_smoothing, name = "smooth_loss")
    }
    model %>% keras::compile(loss = smooth_loss,
                             optimizer = optimizer, metrics = model_metrics)
  } else if (!is.null(label_noise_matrix)) {
    row_sums <- rowSums(label_noise_matrix)
    if (!all(row_sums == 1)) {
      warning("Sum of noise matrix rows don't add up to 1")
    }
    noisy_loss <- noisy_loss_wrapper(solve(label_noise_matrix))
    model %>% keras::compile(loss =  noisy_loss,
                             optimizer = optimizer, metrics = model_metrics)
  } else {
    model %>% keras::compile(loss = "categorical_crossentropy",
                             optimizer = optimizer, metrics = model_metrics)
  }
  
  argg <- c(as.list(environment()))
  model <- add_hparam_list(model, argg)
  reticulate::py_set_attr(x = model, name = "hparam", value = model$hparam)
  
  if (verbose) summary(model)
  return_model <- model
}

#' Extract hyperparameters from model
#'
#' @param model A keras model.
#' @keywords internal
get_hyper_param <- function(model) {
  layers.lstm <- 0
  use.cudnn <- FALSE
  bidirectional <- FALSE
  use.codon.cnn <- FALSE
  learning_rate <- keras::k_eval(model$optimizer$lr)
  solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
  
  layerList <- keras::get_config(model)["layers"]
  for (i in 1:length(layerList)) {
    layer_class_name <- layerList[[i]]$class_name
    
    if (layer_class_name == "Conv1D") {
      use.codon.cnn <- TRUE
    }
    
    if (layer_class_name == "MaxPooling1D") {
    }
    
    if (layer_class_name == "BatchNormalization") {
    }
    
    if (layer_class_name == "CuDNNLSTM") {
      layers.lstm <- layers.lstm + 1
      use.cudnn <- TRUE
      lstm_layer_size <- layerList[[i]]$config$units
      recurrent_dropout_lstm <- 0
      dropout_lstm <- 0
    }
    
    if (layer_class_name == "LSTM") {
      layers.lstm <- layers.lstm + 1
      lstm_layer_size <- layerList[[i]]$config$units
      recurrent_dropout_lstm <- layerList[[i]]$config$recurrent_dropout
      dropout_lstm <- layerList[[i]]$config$dropout
    }
    # TODO: wrong output since bidirectional is layer wrapper (?)
    if (layer_class_name == "Bidirectional") {
      bidirectional <- TRUE
      if (layerList[[i]]$config$layer$class_name == "LSTM") {
        use.cudnn <- FALSE
        layers.lstm <- layers.lstm + 1
        lstm_layer_size <- layerList[[i]]$config$layer$config$units
        recurrent_dropout_lstm <- layerList[[i]]$config$layer$config$recurrent_dropout
        dropout_lstm <- layerList[[i]]$config$layer$config$dropout
      } else {
        use.cudnn <- FALSE
        layers.lstm <- layers.lstm + 1
        lstm_layer_size <- layerList[[i]]$config$layer$config$units
      }
    }
    
    if (layer_class_name == "Dense") {
    }
    
    if (layer_class_name == "Activation") {
    }
  }
  
  list(dropout = dropout_lstm,
       recurrent_dropout = recurrent_dropout_lstm,
       lstm_layer_size =  lstm_layer_size,
       solver = solver,
       use.cudnn = use.cudnn,
       layers.lstm = layers.lstm,
       learning_rate = learning_rate,
       use.codon.cnn = use.codon.cnn,
       bidirectional = bidirectional
  )
}


#' Remove layers from model and add dense layers
#' 
#' Function takes a model as input and removes all layers after a certain layer, specified in \code{layer_name} argument.
#' Optional to add dense layers on top of pruned model. Model can have multiple output layers with separate loss/activation functions.
#' You can freeze all the weights of the pruned model by setting \code{freeze_base_model = TRUE}.
#'
#' @inheritParams create_model_lstm_cnn
#' @param layer_name Name of last layer to use from old model.
#' @param model A keras model. 
#' @param dense_layers List of vectors specifying number of units for each dense layer. If this is a list of length > 1, model
#' has multiple output layers.
#' @param shared_dense_layers Vector with number of units for dense layer. These layers will be connected on top of layer in 
#' argument `layer_name`. Can be used to have shared dense layers, before model has multiple output layers. Don't use if model has just one output layer 
#' (use only `dense_layers`).   
#' @param last_activation List of activations for last entry for each list entry from \code{dense_layers}. Either `"softmax"`, `"sigmoid"` or `"linear"`.
#' @param output_names List of names for each output layer.
#' @param losses List of loss function for each output.
#' @param verbose Boolean.
#' @param dropout List of vectors with dropout rates for each new dense layer.
#' @param dropout_shared Vectors of dropout rates for dense layer from `shared_dense_layers`.
#' @param freeze_base_model Whether to freeze all weights before new dense layers.
#' @param compile Boolean, whether to compile the new model.
#' @param learning_rate Learning rate if `compile = TRUE`, default learning rate of the old model
#' @param global_pooling "max_ch_first" for global max pooling with channel first ([keras docs](https://keras.io/api/layers/pooling_layers/global_average_pooling1d/)),
#' "max_ch_last" for global max pooling with channel last, "average_ch_first" for global average pooling with channel first, 
#' "average_ch_last" for global average pooling with channel last or `NULL` for no global pooling. 
#' "both_ch_first" or "both_ch_last" to combine average and max pooling.
#' @examples
#' model_1 <- create_model_lstm_cnn(layer_lstm = c(64, 64),
#'                                  maxlen = 50,
#'                                  layer_dense = c(32, 4), 
#'                                  verbose = FALSE)
#' # get name of second to last layer 
#' num_layers <- length(model_1$get_config()$layers)
#' layer_name <- model_1$get_config()$layers[[num_layers-1]]$name
#' # add dense layer with multi outputs and separate loss/activations functions
#' model_2 <- remove_add_layers(model = model_1,
#'                              layer_name = layer_name,
#'                              dense_layers = list(c(32, 16, 1), c(8, 1), c(12, 5)),
#'                              losses = list("binary_crossentropy", "mae", "categorical_crossentropy"),
#'                              last_activation = list("sigmoid", "linear", "softmax"),
#'                              freeze_base_model = TRUE,
#'                              output_names = list("out_1_binary_classsification", 
#'                                                  "out_2_regression", 
#'                                                  "out_3_classification")
#') 
#' @export
remove_add_layers <- function(model = NULL,
                              layer_name = NULL,
                              dense_layers = NULL,
                              shared_dense_layers = NULL,
                              last_activation = list("softmax"),
                              output_names = NULL,
                              losses = NULL,
                              verbose = TRUE,
                              dropout = NULL,
                              dropout_shared = NULL,
                              freeze_base_model = FALSE,
                              compile = FALSE,
                              learning_rate = 0.001,
                              solver = "adam",
                              flatten = FALSE,
                              global_pooling = NULL,
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
      model <- do.call(remove_add_layers, argg)
    })
    return(model)
  }
  
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  if (!is.null(layer_name)) check_layer_name(model, layer_name)
  if (!is.null(shared_dense_layers) & is.null(dense_layers)) {
    stop("You need to specify output layers in dense_layers argument")
  }
  if (!is.null(shared_dense_layers) & length(dense_layers) == 1) {
    stop("If your model has just one output layer, use only dense_layers argument (and set shared_dense_layers = NULL).")
  }
  if (!is.null(global_pooling)) {
    stopifnot(global_pooling %in% c("max_ch_first", "max_ch_last", "average_ch_first", "average_ch_last", "both_ch_first", "both_ch_last"))
  }
  
  if (!is.list(dense_layers)) {
    dense_layers <- list(dense_layers)
  }
  
  if (!is.null(dropout) && !is.list(dropout)) {
    dropout <- list(dropout)
  }
  
  if (is.null(losses)) {
    losses <- list()
    for (i in 1:length(last_activation)) {
      if (last_activation[[i]] == "softmax") loss <- "categorical_crossentropy"
      if (last_activation[[i]] == "sigmoid") loss <- "binary_crossentropy"
      if (last_activation[[i]] == "linear") loss <- "mse"
      losses[[i]] <- loss
    }
  }
  
  if (is.null(output_names)) {
    output_names <- vector("list", length = length(dense_layers))
  }
  
  if (length(dense_layers) != length(last_activation)) {
    stop("Length of dense_layers and last_activation must be the same")
  }
  
  if (length(dense_layers) != length(output_names)) {
    stop("Length of dense_layers and output_names must be the same")
  }
  
  if (!is.null(dropout)) {
    for (i in 1:length(dense_layers)) {
      stopifnot(length(dropout[[i]]) == length(dense_layers[[i]]))
    }
  }
  
  if (verbose) {
    print("Original model: ")
    print(model)
  }
  
  is_sequential <- any(stringr::str_detect(class(model), "sequential"))
  
  if (!is.null(layer_name)) {
    model_new <- tensorflow::tf$keras$Model(model$input, model$get_layer(layer_name)$output)
    
    if (freeze_base_model) {
      keras::freeze_weights(model_new)
    }
    
    if (!is.null(global_pooling)) {
      if (global_pooling == "max_ch_first") {
        out <- model_new$output %>% keras::layer_global_max_pooling_1d(data_format="channels_first")
      } else if (global_pooling == "max_ch_last") {
        out <- model_new$output %>% keras::layer_global_max_pooling_1d(data_format="channels_last")
      } else if (global_pooling ==  "average_ch_first") {
        out <- model_new$output %>% keras::layer_global_average_pooling_1d(data_format="channels_first")
      } else if (global_pooling ==  "average_ch_last"){ 
        out <- model_new$output %>% keras::layer_global_average_pooling_1d(data_format="channels_last")
      } else if (global_pooling ==  "both_ch_last"){ 
        out1 <- model_new$output %>% keras::layer_global_average_pooling_1d(data_format="channels_last")
        out2 <- model_new$output %>% keras::layer_global_max_pooling_1d(data_format="channels_last")
        out <- keras::layer_concatenate(list(out1, out2))
      } else {
        out1 <- model_new$output %>% keras::layer_global_average_pooling_1d(data_format="channels_first")
        out2 <- model_new$output %>% keras::layer_global_max_pooling_1d(data_format="channels_first")
        out <- keras::layer_concatenate(list(out1, out2))
      }       
      model_new <- tensorflow::tf$keras$Model(model_new$input, out)
    }
    
    if (flatten) {
      out <- model_new$output %>% keras::layer_flatten()
      model_new <- tensorflow::tf$keras$Model(model_new$input, out)
    }
    
    if (!is.null(shared_dense_layers)) {
      out <- model_new$output 
      for (i in 1:length(shared_dense_layers)) {
        if (!is.null(dropout_shared)) {
          out <- out %>% keras::layer_dropout(dropout_shared[i])
        }
        out <- out %>% keras::layer_dense(shared_dense_layers[i], activation = "relu")
      }
      model_new <- tensorflow::tf$keras$Model(model_new$input, out)
    }
    
    output_list <- list()
    name_dense_index <- 1
    
    if (!is.null(dense_layers[[1]])) {
      for (output_num in 1:length(dense_layers)) {
        for (i in 1:length(dense_layers[[output_num]])) {
          if (i == length(dense_layers[[output_num]])) {
            activation <- last_activation[[output_num]]
            dtype <- "float32"
          } else {
            activation <- "relu"
            dtype <- NULL
          }
          
          if (i == length(dense_layers[[output_num]])) {
            layer_name <- output_names[[output_num]]
          } else {
            layer_name <- paste0("dense_new_", name_dense_index)
            name_dense_index <- name_dense_index + 1
          }
          
          if (is.null(dropout)) {
            if (i == 1) {
              output_list[[output_num]] <- model_new$output %>%
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name, dtype = dtype)
            } else {
              output_list[[output_num]] <- output_list[[output_num]] %>%
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name, dtype = dtype)
            }
          } else {
            if (i == 1) {
              output_list[[output_num]] <- model_new$output %>%
                keras::layer_dropout(rate = dropout[[output_num]][i]) %>%
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name, dtype = dtype)
            } else {
              output_list[[output_num]] <- output_list[[output_num]] %>%
                keras::layer_dropout(rate = dropout[[output_num]][i]) %>%
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name, dtype = dtype)
            }
          }
        }
      }
      model_new <- tensorflow::tf$keras$Model(model_new$input, output_list)
    }
  }
  
  if (verbose) {
    print("New model: ")
    print(model_new)
    for (i in 1:length(model_new$layers)) {
      cat(model_new$layers[[i]]$name , "trainable:" , model_new$layers[[i]]$trainable, "\n")
    }
  }
  
  if (compile) {
    if (is.null(learning_rate)) {
      learning_rate <- keras::k_eval(model$optimizer$lr)
    } else {
      learning_rate <- learning_rate
    }
    
    if (is.null(solver)) {
      solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
    } 
    
    optimizer <- set_optimizer(solver, learning_rate) 
    
    metric_list <- list()
    for (i in 1:length(losses)) {
      metric_list[[i]] <- ifelse(losses[[i]] == "binary_crossentropy", "binary_accuracy", "acc")
    }
    
    model_new %>% keras::compile(loss = losses,
                                 optimizer = optimizer,
                                 metrics = metric_list)
  }
  
  model_new
}


#' Merge two models
#' 
#' Combine two models at certain layers and add dense layer(s) afterwards.
#'
#' @param models List of two models.
#' @param layer_names Vector of length 2 with names of layers to merge.
#' @param freeze_base_model Boolean vector of length 2. Whether to freeze weights of individual models.
#' @inheritParams create_model_lstm_cnn
#' @examples
#' model_1 <- create_model_lstm_cnn(layer_lstm = c(64, 64), maxlen = 50, layer_dense = c(32, 4),
#'                                  verbose = FALSE)
#' model_2 <- create_model_lstm_cnn(layer_lstm = c(32), maxlen = 40, 
#'                                  layer_dense = c(8, 2), verbose = FALSE)
#' # get names of second to last layers
#' num_layers_1 <- length(model_1$get_config()$layers)
#' layer_name_1 <- model_1$get_config()$layers[[num_layers_1 - 1]]$name
#' num_layers_2 <- length(model_2$get_config()$layers)
#' layer_name_2 <- model_2$get_config()$layers[[num_layers_2 - 1]]$name
#' # merge models
#' model <- merge_models(models = list(model_1, model_2),
#'                       layer_names = c(layer_name_1, layer_name_2),
#'                       layer_dense = c(6, 2), 
#'                       freeze_base_model = c(FALSE, FALSE)) 
#' @export
merge_models <- function(models, layer_names, layer_dense, solver = "adam",
                         learning_rate = 0.0001,
                         freeze_base_model = c(FALSE, FALSE),
                         model_seed = NULL) {
  
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  
  model_1 <- remove_add_layers(model = models[[1]],
                               layer_name = layer_names[1],
                               dense_layers = NULL,
                               verbose = FALSE,
                               dropout = NULL,
                               freeze_base_model = freeze_base_model[1],
                               compile = FALSE,
                               learning_rate = NULL)
  
  model_2 <- remove_add_layers(model = models[[2]],
                               layer_name = layer_names[2],
                               dense_layers = NULL,
                               verbose = FALSE,
                               dropout = NULL,
                               freeze_base_model = freeze_base_model[2],
                               compile = FALSE,
                               learning_rate = NULL)
  
  # choose optimization method
  optimizer <- set_optimizer(solver, learning_rate) 
  
  output_tensor <- keras::layer_concatenate(c(model_1$output, model_2$output))
  num_targets <- layer_dense[length(layer_dense)]
  
  if (length(layer_dense) > 1) {
    for (i in 1:(length(layer_dense) - 1)) {
      output_tensor <- output_tensor %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
    }
  }
  
  output_tensor <- output_tensor %>%
    keras::layer_dense(units = num_targets, activation = "softmax")
  
  model <- keras::keras_model(inputs = list(model_1$input, model_2$input), outputs = output_tensor)
  model %>% keras::compile(loss = "categorical_crossentropy",
                           optimizer = optimizer, metrics = c("acc"))
  model
}

#' Check if layer is in model
#' 
#' @keywords internal
check_layer_name <- function(model, layer_name) {
  num_layers <- length(model$layers)
  layer_names <- vector("character")
  for (i in 1:num_layers) {
    layer_names[i] <- model$layers[[i]]$name
  }
  if (!(layer_name %in% layer_names)) {
    message <- paste0("Model has no layer named ", "'", layer_name, "'")
    stop(message)
  }
}

#' Add layer 
#' 
#' Add output of time distribution representations.
#' 
#' @inheritParams create_model_lstm_cnn_time_dist
#' @param load_r6 Whether to load the R6 layer class.
#' @export
layer_add_time_dist_wrapper <- function(load_r6 = FALSE) {
  
  layer_add_time_dist <- keras::new_layer_class(
    "layer_add_time_dist",
    
    initialize = function(...) {
      super$initialize(...)
    },
    
    call = function(inputs) {
      out <- tensorflow::tf$math$reduce_sum(inputs, axis=1L)
      out
    },
    
    get_config = function() {
      config <- super$get_config()
      config
    }
  )
  
  if (load_r6) {
    return(layer_add_time_dist)
  } else {
    return(layer_add_time_dist())
  }
  
}


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
#' is equivalent to [create_model_lstm_cnn_multi_input()] but instead of multiple input layers with 3D input, 
#' input here in one 4D tensor.  
#'     
#' @inheritParams create_model_lstm_cnn
#' @param samples_per_target Number of samples to combine for one target.
#' @param aggregation_method "sum" or "lstm". How to process output of last time distributed layer.
#' If "sum", add representations. If "lstm", use output as input for LSTM layer(s) (as specified in lstm_time_dist).
#' @param lstm_time_dist Vector containing number of units per LSTM cell. Only used if `aggregation_method = "lstm`.   
#' @examples 
#' create_model_lstm_cnn_time_dist(
#'   maxlen = 50,
#'   vocabulary_size = 4,
#'   samples_per_target = 7,
#'   kernel_size = c(10, 10),
#'   filters = c(64, 128),
#'   pool_size = c(2, 2),
#'   layer_lstm = c(32),
#'   layer_dense = c(64, 2),
#'   learning_rate = 0.001)
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
    gap = FALSE,
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
    aggregation_method = "sum", 
    lstm_time_dist = c(128),
    mixed_precision = FALSE,
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
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout_lstm, recurrent_dropout = recurrent_dropout,
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
  
  if (gap) {
    if (layers.lstm != 0) {
      stop("Global average pooling not compatible with using LSTM layer")
    }
    output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_global_average_pooling_1d())
  } else {
    if (layers.lstm == 0) {
      output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_flatten())
    }
  }
  
  if (aggregation_method == "sum") {
    layer_add_td <- layer_add_time_dist_wrapper()
    output_tensor <- output_tensor %>% layer_add_td
  } else {
    return_sequences <- TRUE
    for (i in 1:length(lstm_time_dist)) {
      if (i == length(lstm_time_dist)) return_sequences <- FALSE
      output_tensor <- output_tensor %>% keras::layer_lstm(units=lstm_time_dist[i], return_sequences = return_sequences)
    }
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
  
  # choose optimization method
  optimizer <- set_optimizer(solver, learning_rate) 
  
  #add metrics
  cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  while (dir.exists(cm_dir)) {
    cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  }
  dir.create(cm_dir)
  model$cm_dir <- cm_dir
  
  #add metrics
  if (loss_fn == "binary_crossentropy") {
    model_metrics <- c(tf$keras$metrics$BinaryAccuracy(name = "acc"))
  } else {
    model_metrics <- c("acc")
  } 
  
  if (loss_fn == "categorical_crossentropy") {
    
    macro_average_cb <- balanced_acc_wrapper(num_targets, cm_dir)
    model_metrics <- c(macro_average_cb, "acc")
    
    if (f1_metric) {
      f1 <- f1_wrapper(num_targets = 2, loss = loss_fn)
      model_metrics <- c(model_metrics, f1)
    }
  }
  
  if (auc_metric) {
    auc <- auc_wrapper(model_output_size = layer_dense[length(layer_dense)],
                       loss = loss_fn)
    model_metrics <- c(model_metrics, auc)
  }
  
  model %>% keras::compile(loss = loss_fn,
                           optimizer = optimizer, metrics = model_metrics)
  
  argg <- as.list(environment())
  model <- add_hparam_list(model, argg)
  model$cm_dir <- cm_dir
  
  if (verbose) summary(model)
  return_model <- model
}


#' @title Create LSTM/CNN network that can process multiple samples for one target
#'
#' @description Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers with multiple 
#' input layers. After LSTM/CNN part all representations get aggregated by summation. 
#' Can be used to make single prediction for combination of multiple input sequences. 
#' Implements approach as described [here](https://arxiv.org/abs/1703.06114)
#'
#' @inheritParams create_model_lstm_cnn
#' @param samples_per_target Number of samples to combine for one target.
#' @param dropout_dense Vector of dropout rates between dense layers. No dropout if `NULL`.
#' @examples
#' create_model_lstm_cnn_multi_input(
#'   maxlen = 50,
#'   vocabulary_size = 4,
#'   samples_per_target = 7,
#'   kernel_size = c(10, 10),
#'   filters = c(64, 128),
#'   pool_size = c(2, 2),
#'   layer_lstm = c(32),
#'   layer_dense = c(64, 2),
#'   learning_rate = 0.001)
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
    gap = FALSE,
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
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout, recurrent_dropout = recurrent_dropout_lstm,
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
  
  if (gap) {
    if (layers.lstm != 0) {
      stop("Global average pooling not compatible with using LSTM layer")
    }
    output_tensor <- output_tensor %>% keras::layer_global_average_pooling_1d()
  } else {
    if (layers.lstm == 0) {
      output_tensor <- output_tensor %>% keras::layer_flatten()
    }
  }
  
  #output_tensor <- output_tensor %>% keras::layer_flatten()
  
  feature_ext_model <- keras::keras_model(inputs = input_tensor, outputs = output_tensor)
  
  input_list <- list()
  representation_list <- list()
  for (i in 1:samples_per_target) {
    input_list[[i]] <- keras::layer_input(shape = c(maxlen, vocabulary_size), name = paste0("input_", i))
    representation_list[[i]] <- feature_ext_model(input_list[[i]])
  }
  
  y <- keras::layer_add(representation_list)
  #y <- keras::layer_average(representation_list)
  
  if (length(layer_dense) > 1) {
    for (i in 1:(length(layer_dense) - 1)) {
      if (!is.null(dropout_dense)) y <- y %>% keras::layer_dropout(dropout_dense[i])
      y <- y %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
    }
  }
  
  y <- y %>% keras::layer_dense(units = num_targets, activation = last_layer_activation, dtype = "float32")
  model <- keras::keras_model(inputs = input_list, outputs = y)
  
  if (compile) {
    # choose optimization method
    optimizer <- set_optimizer(solver, learning_rate) 
    
    #add metrics
    cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
    while (dir.exists(cm_dir)) {
      cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
    }
    dir.create(cm_dir)
    model$cm_dir <- cm_dir
    
    #add metrics
    if (loss_fn == "binary_crossentropy") {
      model_metrics <- c(tf$keras$metrics$BinaryAccuracy(name = "acc"))
    } else {
      model_metrics <- c("acc")
    } 
    
    if (loss_fn == "categorical_crossentropy") {
      
      macro_average_cb <- balanced_acc_wrapper(num_targets, cm_dir)
      model_metrics <- c(macro_average_cb, "acc")
      
      if (f1_metric) {
        f1 <- f1_wrapper(num_targets)
        model_metrics <- c(model_metrics, f1)
      }
    }
    
    if (auc_metric) {
      auc <- auc_wrapper(model_output_size = layer_dense[length(layer_dense)],
                         loss = loss_fn)
      model_metrics <- c(model_metrics, auc)
    }
    
    model %>% keras::compile(loss = loss_fn,
                             optimizer = optimizer, metrics = model_metrics)
    
    model$cm_dir <- cm_dir
    
  }
  
  argg <- c(as.list(environment()))
  model <- add_hparam_list(model, argg)
  
  if (verbose) summary(model)
  return_model <- model
}

#' Replace input layer
#'
#' Replace first layer of model with new input layer of different shape. Only works for sequential models that
#' use CNN and LSTM layers.
#'
#' @param model A keras model.
#' @param input_shape The new input shape vector (without batch size).
#' @examples 
#' model_1 <-  create_model_lstm_cnn(
#'   maxlen = 50,
#'   kernel_size = c(10, 10),
#'   filters = c(64, 128),
#'   pool_size = c(2, 2),
#'   layer_lstm = c(32),
#'   verbose = FALSE,
#'   layer_dense = c(64, 2))
#' model <- reshape_input(model_1, input_shape = c(120, 4))
#' model
#' @export
reshape_input <- function(model, input_shape) {
  
  in_layer <- layer_input(shape = input_shape)
  for (i in 2:length(model$layers)) {
    layer_name <- model$layers[[i]]$name
    if (i == 2) {
      out_layer <- in_layer %>% model$get_layer(layer_name)()
    } else {
      out_layer <- out_layer %>% model$get_layer(layer_name)()
    }
  }
  new_model <- tensorflow::tf$keras$Model(in_layer, out_layer)
  return_model <- new_model
}


#' Get solver and learning_rate from model.
#'
#' @keywords internal
get_optimizer <- function(model) {
  solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
  learning_rate <- keras::k_eval(model$optimizer$lr)
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(learning_rate)
  }
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(learning_rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(learning_rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(learning_rate)
  }
  return(optimizer)
}

#' @title Create GenomeNet Model with Given Architecture Parameters
#'
#' @param maxlen (integer `numeric(1)`)\cr
#'   Input sequence length.
#' @param learning.rate (`numeric(1)`)\cr
#'   Used by the `keras` optimizer that is specified by `optimizer`.
#' @param number_of_cnn_layers (integer `numeric(1)`)\cr
#'   Target number of CNN-layers to use in total. If `number_of_cnn_layers` is
#'   greater than `conv_block_count`, then the effective number of CNN layers
#'   is set to the closest integer that is divisible by `conv_block_count`.
#' @param `conv_block_count` (integer `numeric(1)`)\cr
#'   Number of convolutional blocks, into which the CNN layers are divided.
#'   If this is greater than `number_of_cnn_layers`, then it is set to
#'   `number_of_cnn_layers` (the convolutional block size will then be 1).\cr
#'   Convolutional blocks are used when `model_type` is `"gap"` (the output of
#'   the last `conv_block_count * (1 - skip_block_fraction)` blocks is
#'   fed to global average pooling and then concatenated), and also when
#'   `residual_block` is `TRUE` (the number of filters is held constant within
#'   blocks). If neither of these is the case, `conv_block_count` has little
#'   effect besides the fact that `number_of_cnn_layers` is set to the closest
#'   integer divisible by `conv_block_count`.
#' @param `kernel_size_0` (`numeric(1)`)\cr
#'   Target CNN kernel size of the first CNN-layer. Although CNN kernel size is
#'   always an integer, this value can be non-integer, potentially affecting
#'   the kernel-sizes of intermediate layers (which are geometrically
#'   interpolated between `kernel_size_0` and `kernel_size_end`).
#' @param `kernel_size_end` (`numeric(1)`)\cr
#'   Target CNN kernel size of the last CNN-layer; ignored if only one
#'   CNN-layer is used (i.e. if `number_of_cnn_layers` is 1). Although CNN
#'   kernel size is always an integer, this value can be non-integer,
#'   potentially affecting the kernel-sizes of intermediate layers (which are
#'   geometrically interpolated between `kernel_size_0` and `kernel_size_end`).
#' @param filters_0 (`numeric(1)`)\cr
#'   Target filter number of the first CNN-layer. Although CNN filter number is
#'   always an integer, this value can be non-integer, potentially affecting
#'   the filter-numbers of intermediate layers (which are geometrically
#'   interpolated between `filters_0` and `filters_end`).\cr
#'   Note that filters are constant within convolutional blocks when
#'   `residual_block` is `TRUE`.
#' @param filters_end (`numeric(1)`)\cr
#'   Target filter number of the last CNN-layer; ignored if only one CNN-layer
#'   is used (i.e. if `number_of_cnn_layers` is 1). Although CNN filter number
#'   is always an integer, this value can be non-integer, potentially affecting
#'   the filter-numbers of intermediatdilation_rates layers (which are geometrically
#'   interpolated between `kernel_size_0` and `kernel_size_end`).\cr
#'   Note that filters are constant within convolutional blocks when
#'   `residual_block` is `TRUE`.
#' @param dilation_end (`numeric(1)`)\cr
#'   Dilation of the last CNN-layer *within each block*. Dilation rates within
#'   each convolutional block grows exponentially from 1 (no dilation) for the
#'   first CNN-layer to each block, to this value. Set to 1 (default) to
#'   disable dilation.
#' @param max_pool_end (`numeric(1)`)\cr
#'   Target total effective pooling of CNN part of the network. "Effective
#'   pooling" here is the product of the pooling rates of all previous
#'   CNN-layers. A network with three CNN-layers, all of which are followed
#'   by pooling layers of size 2, therefore has effective pooling of 8, with
#'   the effective pooling at intermediate positions being 1 (beginning), 2,
#'   and 4. Effective pooling after each layer is set to the power of 2 that is,
#'   on a logarithmic scale, closest to
#'   `max_pool_end ^ (<CNN layer number> / <total number of CNN layers>)`.
#'   Therefore, even though the total effective pooling size of the whole
#'   CNN part of the network will always be a power of 2, having different,
#'   possibly non-integer values of `max_pool_end`, will still lead to
#'   different networks.
#' @param dense_layer_num (integer `numeric(1)`)\cr
#'   number of dense layers at the end of the network, not counting the output
#'   layer.
#' @param dense_layer_units (integer `numeric(1)`)\cr
#'   Number of units in each dense layer, except for the output layer.
#' @param dropout (`numeric(1)`)\cr
#'   Dropout rate of dense layers, except for the output layer.
#' @param batch_norm_momentum (`numeric(1)`)\cr
#'   `momentum`-parameter of `layer_batch_normalization` layers used in the
#'   convolutional part of the network.
#' @param leaky_relu_alpha (`numeric(1)`)\cr
#'   `alpha`-parameter of the `layer_activation_leaky_relu` activation layers
#'   used in the convolutional part of the network.
#' @param dense_activation (`character(1)`)\cr
#'   Which activation function to use for dense layers. Should be one of
#'   `"relu"`, `"sigmoid"`, or `"tanh"`.
#' @param skip_block_fraction (`numeric(1)`)\cr
#'   What fraction of the first convolutional blocks to skip.
#'   Only used when `model_type` is `"gap"`.
#' @param residual_block (`logical(1)`)\cr
#'   Whether to use residual layers in the convolutional part of the network.
#' @param reverse_encoding (`logical(1)`)\cr
#'   Whether the network should have a second input for reverse-complement
#'   sequences.
#' @param optimizer (`character(1)`)\cr
#'   Which optimizer to use. One of `"adam"`, `"adagrad"`, `"rmsprop"`, or `"sgd"`.
#' @param model_type (`character(1)`)\cr
#'   Whether to use the global average pooling (`"gap"`) or recurrent
#'   (`"recurrent"`) model type.
#' @param recurrent_type (`character(1)`)\cr
#'   Which recurrent network type to use. One of `"lstm"` or `"gru"`.
#'   Only used when `model_type` is `"recurrent"`.
#' @param recurrent_layers (integer `numeric(1)`)\cr
#'   Number of recurrent layers.
#'   Only used when `model_type` is `"recurrent"`.
#' @param recurrent_bidirectional (`logical(1)`)\cr
#'   Whether to use bidirectional recurrent layers.
#'   Only used when `model_type` is `"recurrent"`.
#' @param recurrent_units (integer `numeric(1)`)\cr
#'   Number of units in each recurrent layer.
#'   Only used when `model_type` is `"recurrent"`.
#' @param vocabulary.size (integer `numeric(1)`)\cr
#'   Vocabulary size of (one-hot encoded) input strings. This determines the
#'   input tensor shape, together with `maxlen`.
#' @param last_layer_activation Either `"sigmoid"` or `"softmax"`.
#' @param loss_fn Either `"categorical_crossentropy"` or `"binary_crossentropy"`. If `label_noise_matrix` given, will use custom `"noisy_loss"`.
#' @param num_targets (integer `numeric(1)`)\cr
#'   Number of output units to create.
#' @return A keras model.
#' @export
create_model_genomenet <- function(
    maxlen = 300,
    learning_rate = 0.001,
    number_of_cnn_layers = 1,
    conv_block_count = 1,
    kernel_size_0 = 16,
    kernel_size_end = 16,
    filters_0 = 256,
    filters_end = 512,
    dilation_end = 1,
    max_pool_end = 1,
    dense_layer_num = 1,
    dense_layer_units = 100,
    dropout_lstm = 0,
    batch_norm_momentum = 0.8,
    leaky_relu_alpha = 0,
    dense_activation = "relu",
    skip_block_fraction = 0,
    residual_block = FALSE,
    reverse_encoding = FALSE,
    optimizer = "adam",
    model_type = "gap",
    recurrent_type = "lstm",
    recurrent_layers = 1,
    recurrent_bidirectional = FALSE,
    recurrent_units = 100,
    vocabulary_size = 4,
    last_layer_activation = "softmax",
    loss_fn = "categorical_crossentropy",
    num_targets = 2,
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
      model <- do.call(create_model_genomenet, argg)
    })
    return(model)
  }
  
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  stopifnot(maxlen > 0 & maxlen %% 1 == 0)
  stopifnot(learning_rate > 0)
  stopifnot(number_of_cnn_layers >= 1 & number_of_cnn_layers %% 1 == 0)
  stopifnot(conv_block_count >= 1 & conv_block_count %% 1 == 0)
  stopifnot(kernel_size_0 > 0)
  stopifnot(kernel_size_end > 0)
  stopifnot(filters_0 > 0)
  stopifnot(filters_end > 0)
  stopifnot(dilation_end >= 1)
  stopifnot(max_pool_end >= 1)
  stopifnot(dense_layer_num >= 0 & dense_layer_num %% 1 == 0)
  stopifnot(dense_layer_units >= 0 & dense_layer_units %% 1 == 0)
  stopifnot(0 <= dropout_lstm & dropout_lstm <= 1)
  stopifnot(0 <= batch_norm_momentum & batch_norm_momentum <= 1)
  
  stopifnot(0 <= leaky_relu_alpha& leaky_relu_alpha <= 1)
  dense_activation <- match.arg(dense_activation, c("relu", "sigmoid", "tanh"))
  stopifnot(0 <= skip_block_fraction & skip_block_fraction <= 1)
  
  model_type = match.arg(model_type, c("gap", "recurrent"))
  
  stopifnot(isTRUE(residual_block) || isFALSE(residual_block))
  stopifnot(isTRUE(residual_block) || isFALSE(residual_block))
  
  optimizer <- match.arg(optimizer, c("adam", "adagrad", "rmsprop", "sgd"))
  recurrent_type <- match.arg(recurrent_type, c("lstm", "gru"))
  stopifnot(recurrent_layers >= 1 & recurrent_layers %% 1 == 0)
  stopifnot(isTRUE(recurrent_bidirectional) || isFALSE(recurrent_bidirectional))
  stopifnot(recurrent_units >= 1 && recurrent_units %% 1 == 0)
  
  stopifnot(vocabulary_size >= 2 & vocabulary_size %% 1 == 0)
  stopifnot(num_targets >= 2 & num_targets %% 1 == 0)
  
  if (number_of_cnn_layers < conv_block_count){
    conv_block_size <- 1
    conv_block_count <- number_of_cnn_layers
  } else {
    conv_block_size <- round(number_of_cnn_layers / conv_block_count)
    number_of_cnn_layers <- conv_block_size * conv_block_count
  }
  
  if (residual_block) {
    filters_exponent <- rep(seq(from = 0, to = 1, length.out = conv_block_count), each = conv_block_size)
  } else {
    filters_exponent <- seq(from = 0, to = 1, length.out = number_of_cnn_layers)
  }
  
  filters <- ceiling(filters_0 * (filters_end / filters_0) ^ filters_exponent)
  
  kernel_size <- ceiling(kernel_size_0 * (kernel_size_end / kernel_size_0) ^ seq(from = 0, to = 1, length.out = number_of_cnn_layers))
  
  dilation_rates <- round(dilation_end ^ seq(0, 1, length.out = conv_block_size))
  dilation_rates <- rep(dilation_rates, conv_block_count)
  
  max_pool_divider <- round(log2(max_pool_end) * seq(0, 1, length.out = number_of_cnn_layers + 1))
  max_pool_array <- 2 ^ diff(max_pool_divider)
  
  input_tensor <- keras::layer_input(shape = c(maxlen, vocabulary_size))
  output_tensor <- input_tensor
  
  output_collection <- list()
  
  for (i in seq_len(number_of_cnn_layers)) {
    layer <- keras::layer_conv_1d(kernel_size = kernel_size[i],
                                  padding = "same",
                                  activation = "linear",
                                  filters = filters[i],
                                  dilation_rate = dilation_rates[i])
    
    output_tensor <- output_tensor %>% layer
    
    if (model_type == "gap" && i %% conv_block_size == 0) {
      output_collection[[length(output_collection) + 1]] <- keras::layer_global_average_pooling_1d(output_tensor)
    }
    
    if (max_pool_array[i] > 1) {
      layer <- keras::layer_max_pooling_1d(pool_size = max_pool_array[i], padding = "same")
      output_tensor <- output_tensor %>% layer
    }
    
    if (residual_block) {
      if (i > 1) {
        if (max_pool_array[i] > 1) {
          layer <- keras::layer_average_pooling_1d(pool_size = max_pool_array[i], padding = "same")
          residual_layer <- residual_layer %>% layer
        }
        if (filters[i - 1] != filters[i]) {
          layer <- keras::layer_conv_1d(kernel_size = 1,
                                        padding = "same",
                                        activation = "linear",
                                        filters = filters[i]
          )
          residual_layer <- residual_layer %>% layer
        }
        
        output_tensor <- keras::layer_add(list(output_tensor, residual_layer))
      }
      
      residual_layer <- output_tensor
    }
    
    layer <- keras::layer_batch_normalization(momentum = batch_norm_momentum)
    output_tensor <- output_tensor %>% layer
    
    layer <- keras::layer_activation_leaky_relu(alpha = leaky_relu_alpha)
    
    output_tensor <- output_tensor %>% layer
  }
  
  if (model_type == "gap") {
    # skip 'skip_block_fraction' of outputs we collected --> use the last (1 - skip_block_fraction) part of them
    use_blocks <- ceiling((1 - skip_block_fraction) * length(output_collection))
    use_blocks <- max(use_blocks, 1)
    
    output_collection <- tail(output_collection, use_blocks)
    
    # concatenate outputs from blocks (that we are using)
    if (length(output_collection) > 1) {
      output_tensor <- keras::layer_concatenate(output_collection)
    } else {
      output_tensor <- output_collection[[1]]
    }
  } else {
    # recurrent model
    
    recurrent_layer_constructor = switch(recurrent_type,
                                         lstm = keras::layer_lstm,
                                         gru = keras::layer_gru
    )
    
    for (i in seq_len(recurrent_layers)) {
      if (recurrent_bidirectional) {
        layer <- keras::bidirectional(
          layer = recurrent_layer_constructor(units = recurrent_units, return_sequences = (i != recurrent_layers)))
      } else {
        layer <- recurrent_layer_constructor(units = recurrent_units, return_sequences = (i != recurrent_layers))
      }
      output_tensor <- output_tensor %>% layer
    }
  }
  
  for (i in seq_len(dense_layer_num)) {
    layer <- keras::layer_dropout(rate = dropout_lstm)
    output_tensor <- output_tensor %>% layer
    layer <- keras::layer_dense(units = dense_layer_units, activation = dense_activation)
    output_tensor <- output_tensor %>% layer
  }
  
  output_tensor <- output_tensor %>%
    keras::layer_dense(units = num_targets, activation = last_layer_activation, dtype = "float32")
  
  # define "model" as the mapping from input_tensor to output_tensor
  model <- keras::keras_model(inputs = input_tensor, outputs = output_tensor)
  
  if (reverse_encoding) {
    input_tensor_reversed <- keras::layer_input(shape = c(maxlen, vocabulary_size))
    
    # define "output_tensor_reversed" as what comes out of input_tensor_reversed when model() is applied to it
    output_tensor_reversed <- model(input_tensor_reversed)
    
    # define a new model: model from above (with input_tensor, output_tensor), and
    model <- keras::keras_model(
      inputs = c(input_tensor, input_tensor_reversed),
      outputs = keras::layer_average(c(output_tensor, output_tensor_reversed))
    )
  }
  
  # assign optimization method
  keras_optimizer <- set_optimizer(optimizer, learning_rate) 
  
  #add metrics
  if (loss_fn == "binary_crossentropy") {
    model_metrics <- c(tf$keras$metrics$BinaryAccuracy(name = "acc"))
  } else {
    model_metrics <- c("acc")
  } 
  
  cm_dir <-
    file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  while (dir.exists(cm_dir)) {
    cm_dir <-
      file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  }
  dir.create(cm_dir)
  model$cm_dir <- cm_dir
  
  if (loss_fn == "categorical_crossentropy") {
    if (bal_acc) {
      macro_average_cb <- balanced_acc_wrapper(num_targets, cm_dir)
      model_metrics <- c(macro_average_cb, "acc")
    }
    
    if (f1_metric) {
      f1 <- f1_wrapper(num_targets)
      model_metrics <- c(model_metrics, f1)
    }
  }
  
  if (auc_metric) {
    auc <-
      auc_wrapper(model_output_size = layer_dense[length(layer_dense)],
                  loss = loss_fn)
    model_metrics <- c(model_metrics, auc)
  }
  
  model %>% keras::compile(loss = loss_fn, optimizer = keras_optimizer, metrics = "acc")
  
  argg <- c(as.list(environment()))
  model <- add_hparam_list(model, argg)
  
  model
}


# #' Load pretrained Genomenet model
# #'
# #' Classification model with labels "bacteria", "virus-no-phage","virus-phage".
# #' TODO: add link to paper
# #'
# #' @inheritParams create_model_lstm_cnn
# #' @param maxlen Model input size. Either 150 or 10000.
# #' @param learning_rate Learning rate for optimizer. If compile is TRUE and learning_rate is NULL,
# #' will use learning rate from previous training.
# #' @export
# load_model_self_genomenet <- function(maxlen, compile = FALSE, optimizer = "adam",
#                                       learning_rate = NULL) {
#   
#   stopifnot(any(maxlen == c(150,10000)))
#   
#   if (maxlen == 150) {
#     load(model_self_genomenet_maxlen_150)
#     model <- keras::unserialize_model(model_self_genomenet_maxlen_150, compile = FALSE)
#   }
#   
#   if (maxlen == 10000) {
#     load("data/self_genomenet_model_maxlen_10k.rda")
#     #data(model_self_genomenet_maxlen_10k)
#     model <- keras::unserialize_model(model_self_genomenet_maxlen_10k, compile = FALSE)
#   }
#   
#   if (is.null(learning_rate)) {
#     if (maxlen == 150) learning_rate <- 0.00039517784549691
#     if (maxlen == 10000) learning_rate <- 8.77530464905713e-05
#   }
#   
#   if (compile) {
#     keras_optimizer <- set_optimizer(optimizer, learning_rate)
#     model %>% keras::compile(loss = "categorical_crossentropy", optimizer = keras_optimizer, metrics = "acc")
#   }
#   
#   return_model <- model
# }

set_optimizer <- function(solver = "adam", learning_rate = 0.01) {
  
  stopifnot(solver %in% c("adam", "adagrad", "rmsprop", "sgd"))
  
  named_lr <- "lr" %in% names(formals(keras::optimizer_adam))
  if (named_lr) {
    arg_list <- list(lr = learning_rate)
  } else {
    arg_list <- list(learning_rate = learning_rate)
  }
  
  if (solver == "adam")
    keras_optimizer <- do.call(keras::optimizer_adam, arg_list)
  if (solver == "adagrad")
    keras_optimizer <- do.call(keras::optimizer_adagrad, arg_list)
  if (solver == "rmsprop")
    keras_optimizer <- do.call(keras::optimizer_rmsprop, arg_list)
  if (solver == "sgd")
    keras_optimizer <- do.call(keras::optimizer_sgd, arg_list)
  
  return(keras_optimizer)
  
}

#' Get activation functions of output layers
#' 
#' Get activation functions of output layers.
#' 
#' @param model A keras model.
#' @examples 
#' model <-  create_model_lstm_cnn(
#'   maxlen = 50,
#'   layer_lstm = 8,
#'   layer_dense = c(64, 2),
#'   verbose = FALSE)
#' get_output_activations(model)
#' @export
get_output_activations <- function(model) {
  
  out_names <- model$output_names
  act_vec <- vector("character", length(out_names))
  count <- 1
  
  for (layer_name in out_names) {
    act_name <- model$get_layer(layer_name)$get_config()$activation
    if (is.null(act_name)) act_name <- "linear"
    act_vec[count] <- act_name
    count <- count + 1
  }
  return(act_vec)
}

# temporary fix for metric bugs
manage_metrics <- function(model, compile = FALSE) {
  
  dummy_gen <- generator_dummy(model, batch_size = 1)
  z <- dummy_gen()
  x <- z[[1]]
  y <- z[[2]]
  
  if (length(model$metrics) == 0) {
    suppressMessages(
      eval <- model$evaluate(x, y, verbose = 0L)
    )
  }
  
  if (compile) {
    metric_names <- vector("character", length(model$metrics))
    for (i in 1:length(model$metrics)) {
      metric_names[i] <-  model$metrics[[i]]$name
    }
    
    duplicated_index <- duplicated(metric_names)
    loss_index <- stringr::str_detect(metric_names, "loss")
    index <- duplicated_index | loss_index
    
    # remove double metrics
    model <- model %>% keras::compile(loss = model$loss,
                                      optimizer = model$optimizer,
                                      metrics = model$metrics[!index])
    suppressMessages(
      eval <- model$evaluate(x, y, verbose = 0L)
    )
  }
  
  return(model)
  
}

#' Load checkpoint 
#' 
#' Load checkpoint from directory. Chooses best checkpoint based on some condition. Condition
#' can be best accuracy, best loss, last epoch number or a specified epoch number.
#' 
#' @inheritParams create_model_lstm_cnn
#' @param cp_path A directory containing checkpoints or a single checkpoint file. 
#' If a directory, choose checkpoint based on `cp_filter` or `ep_index`.
#' @param cp_filter Condition to choose checkpoint if `cp_path` is a directory.
#' Either "acc" for bast validation accuracy, "loss" for best validation loss or "last_ep" for last epoch.
#' @param ep_index Load checkpoint from specific epoch number. If not `ǸULL`, has priority over `cp_filter`.
#' @param compile Whether to load compiled model.
#' @param re_compile Whether to compile model with parameters from `learning_rate`,
#' `solver` and `loss`.  
#' @param add_custom_object Named list of custom objects.
#' @param verbose Whether to print chosen checkpoint path.
#' @param loss Loss function. Only used if model gets compiled.
#' @export
load_cp <- function(cp_path, cp_filter = "last_ep", ep_index = NULL, compile = FALSE,
                    learning_rate = 0.01, solver = "adam", re_compile = FALSE,
                    loss = "categorical_crossentropy",
                    add_custom_object = NULL,
                    verbose = TRUE, mirrored_strategy = NULL) {
  
  if (is.null(mirrored_strategy)) mirrored_strategy <- ifelse(count_gpu() > 1, TRUE, FALSE)
  if (mirrored_strategy) {
    mirrored_strategy <- tensorflow::tf$distribute$MirroredStrategy()
    with(mirrored_strategy$scope(), { 
      argg <- as.list(environment())
      argg$mirrored_strategy <- FALSE
      model <- do.call(load_cp, argg)
    })
    return(model)
  }
  
  # custom objects to load transformer architecture
  custom_objects <- list(
    "layer_pos_embedding" = layer_pos_embedding_wrapper(load_r6 = TRUE),
    "layer_pos_sinusoid" = layer_pos_sinusoid_wrapper(load_r6 = TRUE),
    "layer_transformer_block" = layer_transformer_block_wrapper(load_r6 = TRUE),
    "layer_euc_dist" = layer_euc_dist_wrapper(load_r6 = TRUE),
    "layer_cosine_sim" = layer_cosine_sim_wrapper(load_r6 = TRUE),
    "layer_add_time_dist" = layer_add_time_dist_wrapper(load_r6 = TRUE)
  )
  
  if (!is.null(add_custom_object)) {
    for (i in 1:length(add_custom_object)) {
      custom_objects[[names(add_custom_object)[i]]] <- add_custom_object[[i]]
    }
  }
  
  if (!is.null(cp_filter)) {
    stopifnot(cp_filter %in% c("acc", "loss", "last_ep"))
    if (!is.null(ep_index)) {
      cp_filter <- NULL
    }
  } 
  
  is_directory <- dir.exists(cp_path)
  if (is_directory) {
    cps <- list.files(cp_path, full.names = TRUE)
    files_basename <- basename(cps)
    stopifnot(xor(is.null(cp_filter), is.null(ep_index)))
  } else {
    stopifnot(file.exists(cp_path))
    cp <- cp_path
    model <- keras::load_model_hdf5(cp, compile = compile, custom_objects = custom_objects)
  } 
  
  if (is_directory & !is.null(cp_filter)) {
    
    if (cp_filter == "acc") {
      if (!all(stringr::str_detect(files_basename, "acc"))) {
        stop("No accuracy information in checkpoint names ('acc' string), use other metric.")
      }
      acc_scores <- files_basename %>% stringr::str_extract("acc\\d++\\.\\d++") %>% 
        stringr::str_remove("acc") %>% as.numeric()
      # use later epoch for ties
      index <- which.max(rank(acc_scores, ties.method = "last"))
    }
    
    if (cp_filter == "loss") {
      if (!all(stringr::str_detect(files_basename, "loss"))) {
        stop("No loss information in checkpoint names ('loss' string), use other metric.")
      }
      loss_scores <- files_basename %>% stringr::str_extract("loss\\d++\\.\\d++") %>% 
        stringr::str_remove("loss") %>% as.numeric()
      index <- which.min(rank(loss_scores, ties.method = "last"))
    }
    
    if (cp_filter == "last_ep") {
      ep_scores <- files_basename %>% stringr::str_extract("Ep\\.\\d++") %>% 
        stringr::str_remove("Ep\\.") %>% as.numeric()
      index <- which.max(ep_scores)
    }
    
  }
  
  if (is_directory & !is.null(ep_index)) {
    ep_scores <- files_basename %>% stringr::str_extract("Ep\\.\\d++") %>% 
      stringr::str_remove("Ep\\.") %>% as.numeric()
    index <- which(ep_scores == ep_index)
  }
  
  if (is_directory) {
    cp <- cps[index]
    model <- keras::load_model_hdf5(cp, compile = compile, custom_objects = custom_objects)
  }
  
  if (verbose) {
    cat("Using checkpoint", cp, "\n")
  }
  
  if (re_compile) {
    optimizer <- set_optimizer(solver, learning_rate)
    model %>% keras::compile(loss = loss,
                             optimizer = optimizer,
                             metrics = metrics)
  }
  
  return(model)
  
}

get_cp  <- function(cp_path, cp_filter = "last_ep", ep_index = NULL) {
  
  if (!is.null(cp_filter)) {
    stopifnot(cp_filter %in% c("acc", "loss", "last_ep"))
    if (!is.null(ep_index)) {
      cp_filter <- NULL
    }
  } 
  
  is_directory <- dir.exists(cp_path)
  if (is_directory) {
    cps <- list.files(cp_path, full.names = TRUE)
    files_basename <- basename(cps)
    stopifnot(xor(is.null(cp_filter), is.null(ep_index)))
  } else {
    stopifnot(file.exists(cp_path))
    cp <- cp_path
  } 
  
  if (is_directory & !is.null(cp_filter)) {
    
    if (cp_filter == "acc") {
      if (!all(stringr::str_detect(files_basename, "acc"))) {
        stop("No accuracy information in checkpoint names ('acc' string), use other metric.")
      }
      acc_scores <- files_basename %>% stringr::str_extract("acc\\d++\\.\\d++") %>% 
        stringr::str_remove("acc") %>% as.numeric()
      # use later epoch for ties
      index <- which.max(rank(acc_scores, ties.method = "last"))
    }
    
    if (cp_filter == "loss") {
      if (!all(stringr::str_detect(files_basename, "loss"))) {
        stop("No loss information in checkpoint names ('loss' string), use other metric.")
      }
      loss_scores <- files_basename %>% stringr::str_extract("loss\\d++\\.\\d++") %>% 
        stringr::str_remove("loss") %>% as.numeric()
      index <- which.min(rank(loss_scores, ties.method = "last"))
    }
    
    if (cp_filter == "last_ep") {
      ep_scores <- files_basename %>% stringr::str_extract("Ep\\.\\d++") %>% 
        stringr::str_remove("Ep\\.") %>% as.numeric()
      index <- which.max(ep_scores)
    }
    
  }
  
  if (is_directory & !is.null(ep_index)) {
    ep_scores <- files_basename %>% stringr::str_extract("Ep\\.\\d++") %>% 
      stringr::str_remove("Ep\\.") %>% as.numeric()
    index <- which(ep_scores == ep_index)
  }
  
  if (is_directory) {
    cp <- cps[index]
  }
  
  return(cp)
  
}

#' Layer for positional embedding
#' 
#' Positional encoding layer with learned embedding.
#' 
#' @inheritParams create_model_transformer
#' @param load_r6 Whether to load the R6 layer class.
#' @export
layer_pos_embedding_wrapper <- function(maxlen = 100, vocabulary_size = 4, load_r6 = FALSE, embed_dim = 64) {
  
  layer_pos_embedding <- keras::new_layer_class(
    "layer_pos_embedding",
    
    initialize = function(maxlen=100, vocabulary_size=4, embed_dim=64, ...) {
      super$initialize(...)
      if (embed_dim != 0) {
        self$token_emb <- tensorflow::tf$keras$layers$Embedding(input_dim = as.integer(vocabulary_size),
                                                                output_dim = as.integer(embed_dim))
        self$position_embeddings <- tensorflow::tf$keras$layers$Embedding(input_dim = as.integer(maxlen),
                                                                          output_dim = as.integer(embed_dim))
      } else {
        self$position_embeddings <- tensorflow::tf$keras$layers$Embedding(input_dim = as.integer(maxlen),
                                                                          output_dim = as.integer(vocabulary_size))
      }
      self$embed_dim <- as.integer(embed_dim)
      self$maxlen <- as.integer(maxlen)
      self$vocabulary_size <- as.integer(vocabulary_size)
    },
    
    call = function(inputs) {
      positions <- tensorflow::tf$range(self$maxlen, dtype = "int32") 
      embedded_positions <- self$position_embeddings(positions)
      if (self$embed_dim != 0) inputs <- self$token_emb(inputs)
      inputs + embedded_positions
    },
    
    get_config = function() {
      config <- super$get_config()
      config$maxlen <- self$maxlen
      config$vocabulary_size <- self$vocabulary_size
      config$embed_dim <- self$embed_dim
      config
    }
  )
  
  if (load_r6) {
    return(layer_pos_embedding)
  } else {
    return(layer_pos_embedding(maxlen=maxlen, vocabulary_size=vocabulary_size, embed_dim=embed_dim))
  }
  
}

#' Layer for positional encoding
#' 
#' Positional encoding layer with sine/cosine matrix of different frequencies.
#' 
#' @inheritParams create_model_transformer
#' @param load_r6 Whether to load the R6 layer class.
#' @export
layer_pos_sinusoid_wrapper <- function(maxlen = 100, vocabulary_size = 4, n = 10000, load_r6 = FALSE, embed_dim = 64) {
  
  layer_pos_sinusoid <- keras::new_layer_class(
    "layer_pos_sinusoid",
    initialize = function(maxlen, vocabulary_size, n, embed_dim, ...) {
      super$initialize(...)
      self$maxlen <- as.integer(maxlen)
      self$vocabulary_size <- vocabulary_size
      self$n <- as.integer(n)
      self$pe_matrix <- positional_encoding(seq_len = maxlen,
                                            d_model = ifelse(embed_dim == 0,
                                                             as.integer(vocabulary_size),
                                                             as.integer(embed_dim)),  
                                            n = n)
      
      if (embed_dim != 0) {
        self$token_emb <- tensorflow::tf$keras$layers$Embedding(input_dim = vocabulary_size, output_dim = as.integer(embed_dim))
      }
      self$embed_dim <- as.integer(embed_dim)
      
      # self$position_embedding <- tensorflow::tf$keras$layers$Embedding(
      #   input_dim = as.integer(maxlen),
      #   output_dim = ifelse(is.null(embed_dim),
      #                       as.integer(vocabulary_size),
      #                       as.integer(embed_dim)), 
      #   weights = list(pe_matrix),
      #   trainable = FALSE,
      #   name="position_embedding"
      # )(tensorflow::tf$range(start=0, limit=as.integer(maxlen), delta=1))
      
    },
    
    call = function(inputs) {
      if (self$embed_dim != 0) {
        inputs <- self$token_emb(inputs)
      } 
      inputs + self$pe_matrix
    },
    
    get_config = function() {
      config <- super$get_config()
      config$maxlen <- self$maxlen
      config$vocabulary_size <- self$vocabulary_size
      config$n <- self$n
      config$embed_dim <- self$embed_dim
      config$pe_matrix <- self$pe_matrix
      config
    }
  )
  
  if (load_r6) {
    return(layer_pos_sinusoid)
  } else {
    return(layer_pos_sinusoid(maxlen=maxlen, vocabulary_size=vocabulary_size, n=n,
                              embed_dim = embed_dim))
  }
  
}

#' @export
layer_transformer_block_wrapper <- function(num_heads = 2, head_size = 4, dropout_rate = 0, ff_dim = 64,  
                                            vocabulary_size = 4, load_r6 = FALSE, embed_dim = 64) {
  
  layer_transformer_block <- keras::new_layer_class(
    "layer_transformer_block",
    initialize = function(num_heads=2, head_size=4, dropout_rate=0, ff_dim=64L, vocabulary_size=4, embed_dim=64, ...) {
      super$initialize(...)
      self$num_heads <- num_heads
      self$head_size <- head_size
      self$dropout_rate <- dropout_rate
      self$ff_dim <- ff_dim
      self$embed_dim <- as.integer(embed_dim)
      self$vocabulary_size <- vocabulary_size
      self$att <- tensorflow::tf$keras$layers$MultiHeadAttention(num_heads=as.integer(num_heads),
                                                                 key_dim=as.integer(head_size))
      
      self$ffn <- keras::keras_model_sequential() %>% keras::layer_dense(units=as.integer(ff_dim), activation="relu") %>%
        keras::layer_dense(units=ifelse(embed_dim == 0, as.integer(vocabulary_size), as.integer(embed_dim)))
      
      self$layernorm1 <- keras::layer_layer_normalization(epsilon=1e-6)
      self$layernorm2 <- keras::layer_layer_normalization(epsilon=1e-6)
      self$dropout1 <- keras::layer_dropout(rate=dropout_rate)
      self$dropout2 <- keras::layer_dropout(rate=dropout_rate)
    },
    
    call = function(inputs) {
      attn_output <- self$att(inputs, inputs, inputs)
      attn_output <- self$dropout1(attn_output)
      out1 <- self$layernorm1(inputs + attn_output)
      ffn_output <- self$ffn(out1)
      ffn_output <- self$dropout2(ffn_output)
      seq_output <- self$layernorm2(out1 + ffn_output)
      return(seq_output)
    },
    
    get_config = function() {
      config <- super$get_config()
      config$num_heads <- self$num_heads
      config$head_size <- self$head_size
      config$dropout_rate <- self$dropout_rate
      config$ff_dim <- self$ff_dim
      config$vocabulary_size <- self$vocabulary_size
      config$embed_dim <- self$embed_dim
      config
    }
  )
  
  if (load_r6) {
    return(layer_transformer_block)
  } else {
    return(layer_transformer_block(num_heads=num_heads,
                                   head_size=head_size,
                                   dropout_rate=dropout_rate,
                                   vocabulary_size=vocabulary_size,
                                   embed_dim=embed_dim,
                                   ff_dim=ff_dim))
  }
  
}

#' Create transformer model
#'
#' Creates transformer network for classification. Model can consist of several stacked attention blocks.
#' 
#' @inheritParams keras::layer_multi_head_attention
#' @inheritParams create_model_lstm_cnn
#' @param pos_encoding Either `"sinusoid"` or `"embedding"`. How to add positional information.
#' If `"sinusoid"`, will add sine waves of different frequencies to input.
#' If `"embedding"`, model learns positional embedding.
#' @param embed_dim Dimension for token embedding. No embedding if set to 0. Should be used when input is not one-hot encoded
#' (integer sequence).
#' @param head_size Dimensions of attention key.
#' @param n Frequency of sine waves for positional encoding. Only applied if `pos_encoding = "sinusoid"`.
#' @param ff_dim Units of first dense layer after attention blocks.
#' @param dropout Vector of dropout rates after attention block(s). 
#' @param dropout_dense Dropout for dense layers.
#' @param flatten_method How to process output of last attention block. Can be `"gap_channels_last"`, `"gap_channels_first"`, `"none"`,
#' or `"flatten"`. If `"gap_channels_last"` or `"gap_channels_first"`, will apply global average pooling. If `"flatten"`, will flatten 
#' output after last attention block. If `"none"` no flattening applied.
#' @examples 
#' 
#' model <- create_model_transformer(maxlen = 50,
#'                                   head_size=c(10,12), num_heads=c(7,8), ff_dim=c(5,9), 
#'                                   dropout=c(0.3, 0.5))
#' 
#' @export
create_model_transformer <- function(maxlen,
                                     vocabulary_size = 4,
                                     embed_dim = 64,
                                     pos_encoding = "embedding",
                                     head_size = 4L,
                                     num_heads = 5L,
                                     ff_dim = 8,
                                     dropout=0,
                                     n = 10000, # pos emb frequency
                                     layer_dense = 2,
                                     dropout_dense = NULL,
                                     flatten_method = "flatten",
                                     last_layer_activation = "softmax",
                                     loss_fn = "categorical_crossentropy",
                                     solver = "adam",
                                     learning_rate = 0.01,
                                     label_noise_matrix = NULL,
                                     bal_acc = FALSE,
                                     f1_metric = FALSE,
                                     auc_metric = FALSE,
                                     label_smoothing = 0,
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
      model <- do.call(create_model_transformer, argg)
    })
    return(model)
  }
  
  stopifnot(length(head_size) == length(num_heads))
  stopifnot(length(head_size) == length(dropout))
  stopifnot(length(head_size) == length(ff_dim))
  stopifnot(flatten_method %in% c("gap_channels_last", "gap_channels_first", "flatten", "none"))
  stopifnot(pos_encoding %in% c("sinusoid", "embedding"))
  num_dense_layers <- length(layer_dense)
  head_size <- as.integer(head_size)
  num_heads <- as.integer(num_heads)
  maxlen <- as.integer(maxlen)
  num_attention_blocks <- length(num_heads)
  vocabulary_size <-  as.integer(vocabulary_size)
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  
  if (embed_dim == 0) {
    input_tensor <- keras::layer_input(shape = c(maxlen, vocabulary_size))
  } else {
    input_tensor <- keras::layer_input(shape = c(maxlen))
  }
  
  # positional encoding 
  if (pos_encoding == "sinusoid") { 
    pos_enc_layer <- layer_pos_sinusoid_wrapper(maxlen = maxlen, vocabulary_size = vocabulary_size,
                                                n = n, embed_dim = embed_dim)
  } 
  if (pos_encoding == "embedding") { 
    pos_enc_layer <- layer_pos_embedding_wrapper(maxlen = maxlen, vocabulary_size = vocabulary_size,
                                                 embed_dim = embed_dim)
  }
  output_tensor <- input_tensor %>% pos_enc_layer
  
  # attention blocks
  for (i in 1:num_attention_blocks) {
    attn_block <- layer_transformer_block_wrapper(
      num_heads = num_heads[i],
      head_size = head_size[i],
      dropout_rate = dropout[i],
      ff_dim = ff_dim[i],
      embed_dim = embed_dim,
      vocabulary_size = vocabulary_size,
      load_r6 = FALSE)
    output_tensor <- output_tensor %>% attn_block
  }
  
  # flatten
  if (flatten_method == "gap_channels_last") {
    output_tensor <- output_tensor %>% keras::layer_global_average_pooling_1d(data_format="channels_last") 
  }
  if (flatten_method == "gap_channels_first") {
    output_tensor <- output_tensor %>% keras::layer_global_average_pooling_1d(data_format="channels_first") 
  }
  if (flatten_method == "flatten") {
    output_tensor <- output_tensor %>% keras::layer_flatten()
  }
  
  # dense layers
  if (num_dense_layers > 1) {
    for (i in 1:(num_dense_layers - 1)) {
      output_tensor <- output_tensor %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
      if (!is.null(dropout_dense)) {
        output_tensor <- output_tensor %>% keras::layer_dropout(rate = dropout_dense[i])
      }
    }
  }  
  
  output_tensor <- output_tensor %>% keras::layer_dense(units = layer_dense[length(layer_dense)], activation = last_layer_activation, dtype = "float32")
  
  # create model
  model <- keras::keras_model(inputs = input_tensor, outputs = output_tensor)
  
  # compile
  model <- compile_model(model = model, label_smoothing = label_smoothing,
                         solver = solver, learning_rate = learning_rate, loss_fn = loss_fn, 
                         num_output_layers = 1, label_noise_matrix = label_noise_matrix,
                         bal_acc = bal_acc, f1_metric = f1_metric, auc_metric = auc_metric)
  
  if (verbose) print(model)
  
  argg <- c(as.list(environment()))
  model <- add_hparam_list(model, argg)
  reticulate::py_set_attr(x = model, name = "hparam", value = model$hparam)
  
  model
  
}

#' Compile model
#' 
#' @inheritParams create_model_lstm_cnn
#' @export
compile_model <- function(model, solver, learning_rate, loss_fn, label_smoothing = 0,
                          num_output_layers = 1, label_noise_matrix = NULL,
                          bal_acc = FALSE, f1_metric = FALSE, auc_metric = FALSE) {
  
  optimizer <- set_optimizer(solver, learning_rate) 
  if (num_output_layers == 1) num_targets <- model$output_shape[[2]]
  
  #add metrics
  if (loss_fn == "binary_crossentropy") {
    model_metrics <- c(tf$keras$metrics$BinaryAccuracy(name = "acc"))
  } else if (loss_fn == "sparse_categorical_crossentropy") {
    model_metrics <- tensorflow::tf$keras$metrics$SparseCategoricalAccuracy(name = "acc")
  } else {
    model_metrics <- c("acc")
  } 
  
  cm_dir <- NULL
  if (num_output_layers == 1) {
    cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
    while (dir.exists(cm_dir)) {
      cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
    }
    dir.create(cm_dir)
    model$cm_dir <- cm_dir
    
    if (loss_fn == "sparse_categorical_crossentropy") {
      
      if (length(model$outputs) == 1 & length(model$output_shape) == 3) {
        loss_fn <- tensorflow::tf$keras$losses$SparseCategoricalCrossentropy(
          reduction=tensorflow::tf$keras$losses$Reduction$NONE
        )
      }
    }
    
    if (loss_fn == "categorical_crossentropy") {
      
      if (length(model$outputs) == 1 & length(model$output_shape) == 3) {
        loss_fn <- tensorflow::tf$keras$losses$CategoricalCrossentropy(
          label_smoothing=label_smoothing,
          reduction=tensorflow::tf$keras$losses$Reduction$NONE,
          name='categorical_crossentropy'
        )
      }
    }  
    
    if (loss_fn == "categorical_crossentropy" | loss_fn == "sparse_categorical_crossentropy") {
      
      if (bal_acc) {
        macro_average_cb <- balanced_acc_wrapper(num_targets, cm_dir)
        model_metrics <- c(macro_average_cb, "acc")
      }
      
      if (f1_metric) {
        f1 <- f1_wrapper(num_targets)
        model_metrics <- c(model_metrics, f1)
      }
    }
    
    if (auc_metric) {
      auc <- auc_wrapper(model_output_size = layer_dense[length(layer_dense)],
                         loss = loss_fn)
      model_metrics <- c(model_metrics, auc)
    }
    
  }
  
  if (label_smoothing > 0 & !is.null(label_noise_matrix)) {
    stop("Can not use label smoothing and label noise at the same time. Either set label_smoothing = 0 or label_noise_matrix = NULL")
  }
  
  if (label_smoothing > 0) {
    if (loss_fn == "categorical_crossentropy") {
      smooth_loss <- tensorflow::tf$losses$CategoricalCrossentropy(label_smoothing = label_smoothing, name = "smooth_loss")
    }
    if (loss_fn == "binary_crossentropy") {
      smooth_loss <- tensorflow::tf$losses$BinaryCrossentropy(label_smoothing = label_smoothing, name = "smooth_loss")
    }
    model %>% keras::compile(loss = smooth_loss,
                             optimizer = optimizer, metrics = model_metrics)
  } else if (!is.null(label_noise_matrix)) {
    row_sums <- rowSums(label_noise_matrix)
    if (!all(row_sums == 1)) {
      warning("Sum of noise matrix rows don't add up to 1")
    }
    noisy_loss <- noisy_loss_wrapper(solve(label_noise_matrix))
    model %>% keras::compile(loss =  noisy_loss,
                             optimizer = optimizer, metrics = model_metrics)
  } else {
    model %>% keras::compile(loss = loss_fn,
                             optimizer = optimizer, metrics = model_metrics)
  }
  
  model
  
}

get_layer_names <- function(model) {
  
  n <- length(model$layers)
  layer_names <- vector("character", n)
  for (i in 1:n) {
    layer_names[i] <- model$layers[[i]]$name
  }
  layer_names
} 

layer_euc_dist_wrapper <- function(load_r6 = FALSE) {
  
  layer_euc_dist <- keras::new_layer_class(
    "layer_euc_dist",
    
    initialize = function(...) {
      super$initialize(...)
    },
    
    call = function(inputs) {
      euclidean_distance(vects=inputs)
    },
    
    get_config = function() {
      config <- super$get_config()
      config
    }
  )
  
  if (load_r6) {
    return(layer_euc_dist)
  } else {
    return(layer_euc_dist())
  }
  
}


layer_cosine_sim_wrapper <- function(load_r6 = FALSE) {
  
  layer_cosine_sim <- keras::new_layer_class(
    "layer_cosine_sim",
    
    initialize = function(...) {
      super$initialize(...)
    },
    
    call = function(inputs) {
      cosine_similarity(vects=inputs)
    },
    
    get_config = function() {
      config <- super$get_config()
      config
    }
  )
  
  if (load_r6) {
    return(layer_cosine_sim)
  } else {
    return(layer_cosine_sim())
  }
  
}

#' @title Create twin network with contrastive loss
#'
#' @description Twin network can be trained to maximize the distance
#' between embeddings of inputs.
#' Implements approach as described [here](https://keras.io/examples/vision/siamese_contrastive/).
#'
#' @inheritParams create_model_lstm_cnn
#' @param margin Integer, defines the baseline for distance for which pairs should be classified as dissimilar.
#' @param layer_dense Vector containing number of neurons per dense layer, before euclidean distance layer.
#' @param distance_method Either "euclidean" or "cosine".
#' @examples
#' model <- create_model_twin_network(
#'   maxlen = 50,
#'   layer_dense = 16,
#'   kernel_size = 12,
#'   filters = 4,
#'   pool_size = 3,
#'   learning_rate = 0.001,
#'   margin = 1) 
#' @export
create_model_twin_network <- function(
    maxlen = 50,
    dropout_lstm = 0,
    recurrent_dropout_lstm = 0,
    layer_lstm = NULL,
    layer_dense = c(4),
    dropout_dense = NULL,
    kernel_size = NULL,
    filters = NULL,
    strides = NULL,
    pool_size = NULL,
    solver = "adam",
    learning_rate = 0.001,
    vocabulary_size = 4,
    bidirectional = FALSE,
    compile = TRUE,
    padding = "same",
    dilation_rate = NULL,
    gap = FALSE,
    use_bias = TRUE,
    residual_block = FALSE,
    residual_block_length = 1,
    size_reduction_1Dconv = FALSE,
    zero_mask = FALSE,
    margin = 1,
    verbose = TRUE,
    batch_norm_momentum = 0.99,
    distance_method = "euclidean",
    last_layer_activation = "sigmoid",
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
      model <- do.call(create_model_twin_network, argg)
    })
    return(model)
  }
  
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  if (!is.null(dropout_dense)) stopifnot(length(dropout_dense) == length(layer_dense))
  stopifnot(distance_method %in% c("euclidean", "cosine"))
  
  model_base <- create_model_lstm_cnn_multi_input(
    maxlen = maxlen,
    dropout_lstm = dropout_lstm,
    recurrent_dropout_lstm = recurrent_dropout_lstm,
    layer_lstm = layer_lstm,
    solver = solver,
    learning_rate = learning_rate,
    vocabulary_size =  vocabulary_size,
    bidirectional = bidirectional,
    batch_size = NULL,
    compile = FALSE,
    kernel_size = kernel_size,
    filters = filters,
    strides = strides,
    pool_size = pool_size,
    padding = padding,
    dilation_rate = dilation_rate,
    gap = gap,
    use_bias = use_bias,
    zero_mask = zero_mask,
    samples_per_target = 2,
    batch_norm_momentum = batch_norm_momentum,
    verbose = FALSE,
    mixed_precision = mixed_precision,
    mirrored_strategy = FALSE,
    model_seed = model_seed)
  
  model_base <- model_base$layers[[3]]
  input_base <- model_base$input
  
  if (length(layer_dense) > 0) {
    for (i in 1:(length(layer_dense))) {
      if (!is.null(dropout_dense) & i == 1) {
        model_base <- model_base$output %>% keras::layer_dropout(dropout_dense[i])
        model_base <- model_base %>% keras::layer_dense(units = layer_dense[i], activation = "tanh")
      } 
      if (i == 1 & is.null(dropout_dense)) {
        model_base <- model_base$output %>% keras::layer_dense(units = layer_dense[i], activation = "tanh")
      }
      if (i > 1) {
        if (!is.null(dropout_dense)) model_base <- model_base %>% keras::layer_dropout(dropout_dense[i])
        model_base <- model_base %>% keras::layer_dense(units = layer_dense[i], activation = "tanh")
      }
    }
  }
  
  model_base <- keras::keras_model(inputs = input_base, outputs = model_base)
  
  input_1 <- keras::layer_input(shape = c(maxlen, vocabulary_size))
  input_2 <- keras::layer_input(shape = c(maxlen, vocabulary_size))
  tower_1 <- input_1 %>% model_base
  tower_2 <- input_2 %>% model_base
  
  if (distance_method == "euclidean") {
    euc_dist <- layer_euc_dist_wrapper(load_r6 = FALSE)
    outputs <- euc_dist(list(tower_1, tower_2))
  }
  
  if (distance_method == "cosine") {
    cosine_dist <- layer_cosine_sim_wrapper(load_r6 = FALSE)
    outputs <- cosine_dist(list(tower_1, tower_2))
  }
  
  outputs <- outputs %>% keras::layer_batch_normalization(momentum = batch_norm_momentum)
  outputs <- outputs %>% keras::layer_dense(units = 1, activation = last_layer_activation, dtype = "float32")
  model <- keras::keras_model(inputs = list(input_1, input_2), outputs = outputs)
  
  if (compile) {
    model %>% keras::compile(loss = loss_cl(margin=margin),
                             optimizer = set_optimizer(solver, learning_rate),
                             metrics="acc")
  }
  
  if (verbose) print(model)
  model
}

