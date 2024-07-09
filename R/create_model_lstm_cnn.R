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
#' @param dropout_dense Dropout rates between dense layers. No dropout if `NULL`.
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
#' @param last_layer_activation Activation function of output layer(s). For example `"sigmoid"` or `"softmax"`.
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
#' 
#' @returns A keras model, stacks CNN, LSTM and dense layers.   
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
    bal_acc = FALSE,
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
  #browser()
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



#' @title Create LSTM/CNN network to predict middle part of a sequence
#'
#' @description
#' Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers.
#' Function creates two sub networks consisting each of (optional) CNN layers followed by an arbitrary number of LSTM layers. Afterwards the last LSTM layers
#' get concatenated and followed by one or more dense layers. Last layer is a dense layer.
#' Network tries to predict target in the middle of a sequence. If input is AACCTAAGG, input tensors should correspond to x1 = AACC, x2 = GGAA and y = T.
#' 
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
#'  
#' @returns A keras model with two input layers. Consists of LSTN, CNN and dense layers.
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
    auc_metric = FALSE,
    bal_acc = FALSE,
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
  
  if (compile) {
    model <- compile_model(model = model, label_smoothing = label_smoothing, layer_dense = layer_dense,
                           solver = solver, learning_rate = learning_rate, loss_fn = loss_fn, 
                           num_output_layers = num_output_layers, label_noise_matrix = label_noise_matrix,
                           bal_acc = bal_acc, f1_metric = f1_metric, auc_metric = auc_metric)
  }
  
  argg <- c(as.list(environment()))
  model <- add_hparam_list(model, argg)
  reticulate::py_set_attr(x = model, name = "hparam", value = model$hparam)
  
  if (verbose) model$summary()
  return(model)
}
