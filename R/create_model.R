#' @title Creates LSTM/CNN network
#'
#' @description
#' Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers. 
#' Last layer is a dense layer with softmax activation.
#' @param maxlen Length of predictor sequence.
#' @param dropout Fraction of the units to drop for inputs.
#' @param recurrent_dropout Fraction of the units to drop for recurrent state.
#' @param layer_lstm Number of cells per network layer. Can be a scalar or vector.
#' @param layer_dense Dense layers of size layer_dense after last LSTM (or last CNN is \code{layers.lstm = 0}) layer.
#' @param solver Optimization method, options are "adam", "adagrad", "rmsprop" or "sgd".
#' @param learning.rate Learning rate for optimizer.
#' @param use.multiple.gpus If true, multi_gpu_model() will be used based on gpu_num.
#' @param gpu.num Number of GPUs to be used, only relevant if multiple_gpu is true.
#' @param merge.on.cpu True on default, false recommend if the server supports NVlink, only relevant if use.multiple.gpu is true.
#' @param bidirectional Use bidirectional wrapper for lstm layers.
#' @param vocabulary.size Number of unique character in vocabulary.
#' @param stateful Boolean. Whether to use stateful LSTM layer. 
#' @param batch.size Number of samples that are used for one network update. Only used if \code{stateful = TRUE}.
#' @param compile Whether to compile the model. 
#' @param kernel_size Size of 1d convolutional layers. For multiple layers, assign a vector. (e.g, rep(3,2) for two layers and kernel size 3)
#' @param filters Number of filters. For multiple layers, assign a vector.
#' @param strides Stride values. For multiple layers, assign a vector.
#' @param pool_size Integer, size of the max pooling windows. For multiple layers, assign a vector. 
#' @param padding Padding of CNN layers, e.g. "same", "valid" or "causal".
#' @param dilation_rate Integer, the dilation rate to use for dilated convolution.
#' @param gap Whether to apply global average pooling after last CNN layer.
#' @param use_bias Boolean. Usage of bias for CNN layers.
#' @param residual_block Boolean. If true, the residual connections are used in CNN. It is not used in the first convolutional layer.
#' @param residual_block_length Integer. Determines how many convolutional layers (or triplets when size_reduction_1D_conv is TRUE) exist
#  between the legs of a residual connection. e.g. if the length kernel_size/filters is 7 and residual_block_length is 2, there are 1+(7-1)*2 convolutional
#  layers in the model when size_reduction_1Dconv is FALSE and 1+(7-1)*2*3 convolutional layers when size_reduction_1Dconv is TRUE.
#' @param size_reduction_1Dconv Boolean. When TRUE, the number of filters in the convolutional layers is reduced to 1/4 of the number of filters of
#  the original layer by a convolution layer with kernel size 1, and number of filters are increased back to the original value by a convolution layer 
#  with kernel size 1 after the convolution with original kernel size with reduced number of filters.
#' @param label_input Integer or NULL. If not NULL, adds additional input layer of \code{label_input} size.
#' @param zero_mask Boolean, whether to apply zero masking before LSTM layer. Only used if model does not use any CNN layers. 
#' @param label_smoothing Float in [0, 1]. If 0, no smoothing is applied. If > 0, loss between the predicted 
#' labels and a smoothed version of the true labels, where the smoothing squeezes the labels towards 0.5.
#' The closer the argument is to 1 the more the labels get smoothed.
#' @param label_noise_matrix Matrix of label noises. Every row stands for one class and columns for percentage of labels in that class.
#' If first label contains 5% wrong labels and second label no noise, then
#' \code{label_noise_matrix <- matrix(c(0.95, 0.05, 0, 1), nrow = 2, by_row = TRUE )}  
#' @param last_layer_activation Either "sigmoid" or "softmax".  
#' @param loss_fn Either "categorical_crossentropy" or "binary_crossentropy". If label_noise_matrix given, will use custom "noisy_loss".
#' @param num_output_layers Number of output layers. 
#' @param auc_metric Whether to add AUC metric.
#' @param f1_metric Whether to add F1 metric.
#' @param bal_acc Whether to add balanced accuracy. 
#' @export
create_model_lstm_cnn <- function(
  maxlen = 50,
  dropout = 0,
  recurrent_dropout = 0,
  layer_lstm = NULL,
  layer_dense = c(4),
  solver = "adam",
  learning.rate = 0.001,
  use.multiple.gpus = FALSE,
  merge.on.cpu = TRUE,
  gpu.num = 2,
  num_targets = 4,  # for compatibility
  vocabulary.size = 4,
  bidirectional = FALSE,
  stateful = FALSE,
  batch.size = NULL,
  compile = TRUE,
  kernel_size = NULL,
  filters = NULL,
  strides = NULL,
  pool_size = NULL,
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
  bal_acc = TRUE
) {
  
  num_targets <- layer_dense[length(layer_dense)]
  layers.lstm <- length(layer_lstm)
  use.cnn <- ifelse(!is.null(kernel_size), TRUE, FALSE)
  
  if (!is.null(layer_lstm)) {
    stopifnot(length(layer_lstm) == 1 | (length(layer_lstm) ==  layers.lstm))
  }
  
  if (layers.lstm == 0 & !use.cnn) {
    stop("Model does not use LSTM or CNN layers. Set use.cnn to TRUE or layers.lstm > 0.")
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
  stopifnot(dropout <= 1 & dropout >= 0)
  stopifnot(recurrent_dropout <= 1 & recurrent_dropout >= 0)
  
  if (length(layer_lstm) == 1) {
    layer_lstm <- rep(layer_lstm, layers.lstm)
  }
  
  if (stateful) {
    input_tensor <- keras::layer_input(batch_shape = c(batch.size, maxlen, vocabulary.size))
  } else {
    input_tensor <- keras::layer_input(shape = c(maxlen, vocabulary.size))
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
            input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
            use_bias = use_bias
          )
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
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
                input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
                use_bias = use_bias
              )
            residual_layer <- residual_layer %>% keras::layer_batch_normalization(momentum = .8)
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
                    input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
                    use_bias = use_bias
                  )
                output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
                
                output_tensor <- output_tensor %>%
                  keras::layer_conv_1d(
                    kernel_size = kernel_size[i],
                    padding = padding,
                    activation = "relu",
                    filters = filters[i]/4,
                    strides = 1,
                    dilation_rate = dilation_rate[i],
                    input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
                    use_bias = use_bias
                  )
                output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
                
                output_tensor <- output_tensor %>%
                  keras::layer_conv_1d(
                    kernel_size = 1,
                    padding = padding,
                    activation = "relu",
                    filters = filters[i],
                    strides = 1,
                    dilation_rate = dilation_rate[i],
                    input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
                    use_bias = use_bias
                  )
                output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
                
              } else {
                output_tensor <- output_tensor %>%
                  keras::layer_conv_1d(
                    kernel_size = kernel_size[i],
                    padding = padding,
                    activation = "relu",
                    filters = filters[i],
                    strides = 1,
                    dilation_rate = dilation_rate[i],
                    input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
                    use_bias = use_bias
                  )
                output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)    
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              use_bias = use_bias
            )
          output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
          
          output_tensor <- output_tensor %>%
            keras::layer_conv_1d(
              kernel_size = kernel_size[i],
              padding = padding,
              activation = "relu",
              filters = filters[i]/4,
              strides = strides[i],
              dilation_rate = dilation_rate[i],
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              use_bias = use_bias
            )
          output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
          
          output_tensor <- output_tensor %>%
            keras::layer_conv_1d(
              kernel_size = 1,
              padding = padding,
              activation = "relu",
              filters = filters[i],
              strides = 1,
              dilation_rate = dilation_rate[i],
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              use_bias = use_bias
            )
          output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
          
        } else {
          output_tensor <- output_tensor %>%
            keras::layer_conv_1d(
              kernel_size = kernel_size[i],
              padding = padding,
              activation = "relu",
              filters = filters[i],
              strides = strides[i],
              dilation_rate = dilation_rate[i],
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              use_bias = use_bias
            )
          output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)    
        }
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        #output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              keras::layer_lstm(
                units = layer_lstm[i],
                return_sequences = TRUE,
                dropout = dropout,
                recurrent_dropout = recurrent_dropout,
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              return_sequences = TRUE,
              dropout = dropout,
              recurrent_dropout = recurrent_dropout,
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
          input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout, recurrent_dropout = recurrent_dropout, 
                            stateful = stateful, recurrent_activation = "sigmoid")
        )
    } else {
      output_tensor <- output_tensor %>%
        keras::layer_lstm(units = layer_lstm[length(layer_lstm)], 
                          input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
                          dropout = dropout, recurrent_dropout = recurrent_dropout, stateful = stateful,
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
        eval(parse(text = paste0("label_input_layer_", as.character(i), "<- keras::layer_input(batch_shape = c(batch.size, label_input[i]))")))
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
      keras::layer_dense(units = num_targets, activation = last_layer_activation) 
  } else {
    output_list <- list()
    for (i in 1:num_output_layers) {
      layer_name <- paste0("output_", i, "_", num_output_layers)
      output_list[[i]] <- output_tensor %>%
        keras::layer_dense(units = num_targets, activation = last_layer_activation, name = layer_name) 
    }
  }
  
  # print model layout to screen, should be done before multi_gpu_model
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
  
  if (use.multiple.gpus) {
    model <- keras::multi_gpu_model(model,
                                    gpus = gpu.num,
                                    cpu_merge = merge.on.cpu)
  }
  
  # choose optimization method
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(lr = learning.rate)
  }
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning.rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning.rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning.rate)
  }
  
  #add metrics
  metrics <- c("acc")
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
        metrics <- c(macro_average_cb, "acc")
      }
      
      if (f1_metric) {
        f1 <- f1_wrapper(num_targets)
        metrics <- c(metrics, f1)
      }
    }
    
    if (auc_metric) {
      # multi_label only for tf > 2.2
      multi_label <- ifelse(loss_fn == "binary_crossentropy", TRUE, FALSE)
      #auc <- tensorflow::tf$keras$metrics$AUC(multi_label = multi_label)
      auc <- tensorflow::tf$keras$metrics$AUC()
      metrics <- c(metrics, auc)
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
                             optimizer = optimizer, metrics = metrics)
  } else if (!is.null(label_noise_matrix)) {
    row_sums <- rowSums(label_noise_matrix)
    if (!all(row_sums == 1)) {
      warning("Sum of noise matrix rows don't add up to 1")
    }
    noisy_loss <- noisy_loss_wrapper(solve(label_noise_matrix))
    model %>% keras::compile(loss =  noisy_loss,
                             optimizer = optimizer, metrics = metrics)
  } else {
    model %>% keras::compile(loss = loss_fn,
                             optimizer = optimizer, metrics = metrics)
  }
  
  argg <- c(as.list(environment()))
  argg["metrics"] <- NULL
  argg["model"] <- NULL
  argg["i"] <- NULL
  argg["optimizer"] <- NULL
  argg["layer_lstm"] <- paste(as.character(layer_lstm), collapse = " ")
  argg["filters"] <- paste(as.character(filters), collapse = " ")
  argg["kernel_size"] <- paste(as.character(kernel_size), collapse = " ")
  argg["pool_size"] <- paste(as.character(pool_size), collapse = " ")
  argg["strides"] <- paste(as.character(strides), collapse = " ")
  argg["residual_block"] <- paste(as.character(residual_block), collapse = " ")
  argg["residual_block_length"] <- paste(as.character(residual_block_length), collapse = " ")
  argg["size_reduction_1Dconv"] <- paste(as.character(size_reduction_1Dconv), collapse = " ")
  argg["layer_dense"] <- paste(as.character(layer_dense), collapse = " ")
  argg["padding"] <- paste(as.character(padding), collapse = " ")
  argg["use_bias"] <- paste(as.character(use_bias), collapse = " ")
  argg["input_label_list"] <- paste(as.character(layer_dense), collapse = " ")
  argg["input_tensor"] <- NULL
  argg["label_inputs"] <- NULL
  argg["f1"] <- NULL
  argg["multi_acc"] <- NULL
  argg[["trainable_params"]] <- model$count_params()
  for (i in 1:length(label_input)) {
    argg[paste0("input_tensor_", i)] <- NULL
    argg[paste0("label_input_layer_", i)] <- NULL
  }
  argg["output_tensor"] <- NULL
  argg["output_list"] <- NULL
  argg["residual_layer"] <- NULL
  argg["label_noise_matrix"] <- NULL
  argg["smooth_loss"] <- NULL
  argg["noisy_loss"] <- NULL
  argg["col_sums"] <- NULL
  argg["auc"] <- NULL
  argg["multi_label"] <- NULL
  argg["macro_average_cb"] <- NULL
  model$hparam <- argg
  model$cm_dir <- cm_dir
  
  summary(model)
  return(model)
}

#' create wavenet model
#'  
#' @inheritParams wavenet::wavenet
#' @export
create_model_wavenet <- function(filters = 16, kernel_size = 2, residual_blocks, maxlen,
                                 input_tensor = NULL, initial_kernel_size = 32, initial_filters = 32,
                                 output_channels = 4, output_activation = "softmax", solver = "adam",
                                 learning.rate = 0.001, compile = TRUE) {
  
  model <- wavenet::wavenet(filters = filters, kernel_size = kernel_size, residual_blocks = residual_blocks,
                            input_shape = list(maxlen, output_channels), input_tensor = input_tensor, initial_kernel_size = initial_kernel_size, 
                            initial_filters = initial_filters, output_channels = output_channels, output_activation = "softmax")
  if (solver == "adam") {
    optimizer <- keras::optimizer_adam(lr = learning.rate)
  }  
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning.rate)
  }  
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning.rate)
  }  
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning.rate) 
  }
  
  if (compile) {
    model %>% keras::compile(loss = "categorical_crossentropy",
                             optimizer = optimizer, metrics = c("acc"))
  }
  
  argg <- c(as.list(environment()))
  argg["model"] <- NULL
  argg["optimizer"] <- NULL
  argg["residual_blocks"] <- paste(as.character(residual_blocks), collapse = " ")
  model$hparam <- argg
  reticulate::py_set_attr(x = model, name = "hparam", value = model$hparam)
  model
}

#######

#' @title Creates LSTM/CNN network
#'
#' @description
#' Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers. 
#' Last layer is a dense layer with softmax activation.
#' Network tries to predict target in the middle of a sequence. If input is AACCTAAGG, input tensors should correspond to x1 = AACC, x2 = GGAA and y = T 
#' (set \code{target_middle = TRUE} in \code{trainNetwork} function).  
#' Function creates two sub networks consisting each of an (optional) CNN layer followed by an arbitrary number of LSTM layers. Afterwards the last LSTM layers
#' get concatenated and followed by a dense layers. 
#' @param maxlen Length of predictor sequence.
#' @param dropout Fraction of the units to drop for inputs.
#' @param recurrent_dropout Fraction of the units to drop for recurrent state.
#' @param layer_lstm Number of cells per network layer. Can be a scalar or vector.
#' @param solver Optimization method, options are "adam", "adagrad", "rmsprop" or "sgd".
#' @param learning.rate Learning rate for optimizer.
#' @param use.multiple.gpus If true, multi_gpu_model() will be used based on gpu_num.
#' @param gpu.num Number of GPUs to be used, only relevant if multiple_gpu is true.
#' @param merge.on.cpu True on default, false recommend if the server supports NVlink, only relevant if use.multiple.gpu is true.
#' @param bidirectional Use bidirectional wrapper for lstm layers.
#' @param vocabulary.size Number of unique character in vocabulary.
#' @param stateful Boolean. Whether to use stateful LSTM layer. 
#' @param batch.size Number of samples that are used for one network update. Only used if \code{stateful = TRUE}.
#' @param padding Padding of CNN layers, e.g. "same", "valid" or "causal".
#' @param compile Whether to compile the model. 
#' @param layer_dense Dense layers of size layer_dense after last LSTM (or last CNN if \code{layers.lstm = 0}) layer.
#' @param kernel_size Size of 1d convolutional layer.
#' @param filters Number of filters.
#' @param pool_size Integer, size of the max pooling windows.
#' @param strides Stide length of convolution.
#' @param label_input Integer or NULL. If not NULL, adds additional input layer of \code{label_input} size.
#' @export
create_model_lstm_cnn_target_middle <- function(
  maxlen = 50,
  dropout = 0,
  recurrent_dropout = 0,
  layer_lstm = 128,
  solver = "adam",
  learning.rate = 0.001,
  use.multiple.gpus = FALSE,
  merge.on.cpu = TRUE,
  gpu.num = 2,
  vocabulary.size = 4,
  bidirectional = FALSE,
  stateful = FALSE,
  batch.size = NULL,
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
  num_output_layers = 1
) {
  
  use.cnn <- ifelse(!is.null(kernel_size), TRUE, FALSE)
  num_targets <- layer_dense[length(layer_dense)]
  layers.lstm <- length(layer_lstm)
  
  stopifnot(length(layer_lstm) == 1 | (length(layer_lstm) ==  layers.lstm))
  stopifnot(maxlen > 0)
  stopifnot(dropout <= 1 & dropout >= 0)
  stopifnot(recurrent_dropout <= 1 & recurrent_dropout >= 0)
  
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
    input_tensor_1 <- keras::layer_input(batch_shape = c(batch.size, maxlen_1, vocabulary.size))
  } else {
    input_tensor_1 <- keras::layer_input(shape = c(maxlen_1, vocabulary.size))
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
            input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL)
          ) 
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor_1 <- output_tensor_1 %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        output_tensor_1 <- output_tensor_1 %>% keras::layer_batch_normalization(momentum = .8)
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
        output_tensor_1 <- output_tensor_1 %>% keras::layer_batch_normalization(momentum = .8)
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              keras::layer_lstm(
                units = layer_lstm[i],
                return_sequences = TRUE,
                dropout = dropout,
                recurrent_dropout = recurrent_dropout,
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              return_sequences = TRUE,
              dropout = dropout,
              recurrent_dropout = recurrent_dropout,
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
          input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout, recurrent_dropout = recurrent_dropout, 
                            stateful = stateful, recurrent_activation = "sigmoid")
        )
    } else {
      output_tensor_1 <- output_tensor_1 %>%
        keras::layer_lstm(units = layer_lstm[length(layer_lstm)], 
                          input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
                          dropout = dropout, recurrent_dropout = recurrent_dropout, stateful = stateful,
                          recurrent_activation = "sigmoid")
    }
  }
  
  if (stateful) {
    input_tensor_2 <- keras::layer_input(batch_shape = c(batch.size, maxlen_2, vocabulary.size))
  } else {
    input_tensor_2 <- keras::layer_input(shape = c(maxlen_2, vocabulary.size))
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
            input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL)
          )  
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor_2 <- output_tensor_2 %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        output_tensor_2 <- output_tensor_2 %>% keras::layer_batch_normalization(momentum = .8)
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
        output_tensor_2 <- output_tensor_2 %>% keras::layer_batch_normalization(momentum = .8)
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              keras::layer_lstm(
                units = layer_lstm[i],
                return_sequences = TRUE,
                dropout = dropout,
                recurrent_dropout = recurrent_dropout,
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              return_sequences = TRUE,
              dropout = dropout,
              recurrent_dropout = recurrent_dropout,
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
          input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout, recurrent_dropout = recurrent_dropout, 
                            stateful = stateful, recurrent_activation = "sigmoid")
        )
    } else {
      output_tensor_2 <- output_tensor_2 %>%
        keras::layer_lstm(units = layer_lstm[length(layer_lstm)], 
                          input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
                          dropout = dropout, recurrent_dropout = recurrent_dropout, stateful = stateful,
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
        eval(parse(text = paste0("label_input_layer_", as.character(i), "<- keras::layer_input(batch_shape = c(batch.size, label_input[i]))")))
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
      keras::layer_dense(units = num_targets, activation = last_layer_activation) 
  }  else {
    output_list <- list()
    for (i in 1:num_output_layers) {
      layer_name <- paste0("output_", i, "_", num_output_layers)
      output_list[[i]] <- output_tensor %>%
        keras::layer_dense(units = num_targets, activation = last_layer_activation, name = layer_name) 
    }
  }
  
  # print model layout to screen, should be done before multi_gpu_model
  if (!is.null(label_input)) {
    label_inputs <- list()
    for (i in 1:length(label_input)) {
      eval(parse(text = paste0("label_inputs$label_input_layer_", as.character(i), "<- label_input_layer_", as.character(i))))
    }
    model <- keras::keras_model(inputs = c(label_inputs, input_tensor_1, input_tensor_2), outputs = output_tensor)
  } else {
    model <- keras::keras_model(inputs = list(input_tensor_1, input_tensor_2), outputs = output_tensor)
  }
  
  if (use.multiple.gpus) {
    model <- keras::multi_gpu_model(model,
                                    gpus = gpu.num,
                                    cpu_merge = merge.on.cpu)
  }
  
  # choose optimization method
  if (solver == "adam")
    optimizer <-
    keras::optimizer_adam(lr = learning.rate)
  if (solver == "adagrad")
    optimizer <-
    keras::optimizer_adagrad(lr = learning.rate)
  if (solver == "rmsprop")
    optimizer <-
    keras::optimizer_rmsprop(lr = learning.rate)
  if (solver == "sgd")
    optimizer <-
    keras::optimizer_sgd(lr = learning.rate)
  
  #add metrics
  cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  while (dir.exists(cm_dir)) {
    cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  }
  dir.create(cm_dir)
  model$cm_dir <- cm_dir
  metrics <- c("acc")
  
  if (loss_fn == "categorical_crossentropy") {
    
    macro_average_cb <- balanced_acc_wrapper(num_targets, cm_dir)
    metrics <- c(macro_average_cb, "acc")
    
    if (f1_metric) {
      f1 <- f1_wrapper(num_targets)
      metrics <- c(metrics, f1)
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
                             optimizer = optimizer, metrics = c("acc", f1))
  } else if (!is.null(label_noise_matrix)) {
    row_sums <- rowSums(label_noise_matrix)
    if (!all(row_sums == 1)) {
      warning("Sum of noise matrix rows don't add up to 1")
    }
    noisy_loss <- noisy_loss_wrapper(solve(label_noise_matrix))
    model %>% keras::compile(loss =  noisy_loss,
                             optimizer = optimizer, metrics = c("acc", f1))
  } else {
    model %>% keras::compile(loss = "categorical_crossentropy",
                             optimizer = optimizer, metrics = c("acc", f1))
  }
  
  argg <- c(as.list(environment()))
  argg["metrics"] <- NULL
  argg["i"] <- NULL
  argg["optimizer"] <- NULL
  argg["model"] <- NULL
  argg["input_tensor_1"] <- NULL
  argg["input_tensor_2"] <- NULL
  argg["input_label_list"] <- NULL
  for (i in 1:length(label_input)) {
    argg[paste0("input_tensor_", i)] <- NULL
    argg[paste0("label_input_layer_", i)] <- NULL
  }
  argg[["trainable_params"]] <- model$count_params()
  argg["label_input"] <- NULL
  argg["label_inputs"] <- NULL
  argg["maxlen_1"] <- NULL
  argg["maxlen_2"] <- NULL
  argg["f1"] <- NULL
  argg["output_tensor"] <- NULL
  argg["output_tensor_1"] <- NULL
  argg["output_tensor_2"] <- NULL
  argg["output_list"] <- NULL
  argg["optimizer"] <- NULL
  argg["label_noise_matrix"] <- NULL
  argg["smooth_loss"] <- NULL
  argg["noisy_loss"] <- NULL
  argg["col_sums"] <- NULL
  argg["layer_lstm"] <- paste(as.character(layer_lstm), collapse = " ")
  argg["filters"] <- paste(as.character(filters), collapse = " ")
  argg["kernel_size"] <- paste(as.character(kernel_size), collapse = " ")
  argg["strides"] <- paste(as.character(strides), collapse = " ")
  argg["layer_dense"] <- paste(as.character(layer_dense), collapse = " ")
  argg["macro_average_cb"] <- NULL
  model$hparam <- argg
  reticulate::py_set_attr(x = model, name = "hparam", value = model$hparam)
  
  summary(model)
  return(model)
}

#' Extract hyperparameters from model
#'
#' @param model A keras model.
#' @export
get_hyper_param <- function(model) {
  layers.lstm <- 0
  use.cudnn <- FALSE
  bidirectional <- FALSE
  use.codon.cnn <- FALSE
  learning.rate <- keras::k_eval(model$optimizer$lr)
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
      recurrent_dropout <- 0
      dropout <- 0
    }
    
    if (layer_class_name == "LSTM") {
      layers.lstm <- layers.lstm + 1
      lstm_layer_size <- layerList[[i]]$config$units
      recurrent_dropout <- layerList[[i]]$config$recurrent_dropout
      dropout <- layerList[[i]]$config$dropout
    }
    # TODO: wrong output since bidirectional is layer wrapper (?)
    if (layer_class_name == "Bidirectional") {
      bidirectional <- TRUE
      if (layerList[[i]]$config$layer$class_name == "LSTM") {
        use.cudnn <- FALSE
        layers.lstm <- layers.lstm + 1
        lstm_layer_size <- layerList[[i]]$config$layer$config$units
        recurrent_dropout <- layerList[[i]]$config$layer$config$recurrent_dropout
        dropout <- layerList[[i]]$config$layer$config$dropout
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
  
  list(dropout = dropout,
       recurrent_dropout = recurrent_dropout,
       lstm_layer_size =  lstm_layer_size,
       solver = solver,
       use.cudnn = use.cudnn,
       layers.lstm = layers.lstm,
       learning.rate = learning.rate,
       use.codon.cnn = use.codon.cnn,
       bidirectional = bidirectional
  )
}


#' Remove layers from model and add dense layers
#'
#' @inheritParams create_model_lstm_cnn
#' @param layer_name Name of last layer to use from old model.
#' @param model A keras model. If model and model_path are both not NULL, model_path will be used.
#' @param dense_layers List of vectors specifiying number of units for each dense layer. If this is a list of length > 1, model 
#' has multiple output layers.
#' @param last_activation List of activations for last entry for each list entry from \code{dense_layers}. Either "softmax" or "sigmoid".    
#' @param output_names List of names for each output layer.
#' @param losses List of loss function for each output.
#' @param verbose Boolean.
#' @param dropout List of vectors with dropout rates for each new dense layer.
#' @param freeze_base_model Whether to freeze all weights before new dense layers.  
#' @param compile Boolean, whether the new model is compiled or not
#' @param lr learning rate if compile == TRUE, default -> learning rate of the old model
#' @param solver "adam", "adagrad", "rmsprop" or "sgd" if compile == TRUE, default -> solver of the old model
#' @export
remove_add_layers <- function(model = NULL,
                              layer_name = NULL,
                              layer_depth = NULL,
                              dense_layers = NULL,
                              last_activation = list("softmax"),
                              output_names = NULL,
                              losses = NULL,
                              verbose = TRUE,
                              dropout = NULL,
                              freeze_base_model = FALSE,
                              compile = FALSE,
                              lr = 0.001,
                              solver = "adam") {
  
  if (!is.null(layer_name)) check_layer_name(model, layer_name)
  
  if (!is.list(dense_layers)) {
    dense_layers <- list(dense_layers)
  }
  
  if (is.null(losses)) {
    losses <- list()
    for (i in 1:length(last_activation)) {
      loss <- ifelse(last_activation[[i]] == "softmax", "categorical_crossentropy", "binary_crossentropy")
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
    stopifnot(length(dropout) == length(dense_layers))
  }
  
  if (verbose) {
    print("Original model: ")
    print(model)
  }
  
  is_sequential <- any(stringr::str_detect(class(model), "sequential"))
  if (!is_sequential & is.null(layer_name) & !is.null(layer_depth)) {
    stop("Model is not sequential, use layer_name argument instead of layer_depth.")
  }
  
  if (!is.null(layer_name)) {
    model_new <- tensorflow::tf$keras$Model(model$input, model$get_layer(layer_name)$output)
    
    if (freeze_base_model) {
      keras::freeze_weights(model_new)
    }
    
    output_list <- list()
    
    if (!is.null(dense_layers[[1]])) {
      for (output_num in 1:length(dense_layers)) {
        for (i in 1:length(dense_layers[[output_num]])) {
          if (i == length(dense_layers[[output_num]])) {
            activation <- last_activation[[output_num]]
          } else {
            activation <- "relu"
          }
          
          if (i == length(dense_layers[[output_num]])) {
            layer_name <- output_names[[output_num]]
          } else {
            layer_name <- NULL
          }
          
          
          if (is.null(dropout)) {
            if (i == 1) {
              output_list[[output_num]] <- model_new$output %>%
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name)
            } else {
              output_list[[output_num]] <- output_list[[output_num]] %>%
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name)
            } 
          } else {
            if (i == 1) {
              output_list[[output_num]] <- model_new$output %>%
                keras::layer_dropout(rate = dropout[[output_num]][i]) %>% 
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name)
            } else {
              output_list[[output_num]] <- output_list[[output_num]] %>%
                keras::layer_dropout(rate = dropout[[output_num]][i]) %>% 
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name)
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
    layerList <- keras::get_config(model_new)["layers"]
    for (i in 1:length(layerList)) {
      cat(layerList[[i]]$config$name , "trainable:" , layerList[[i]]$config$trainable, "\n")
    }
  }
  
  if (compile) {
    if (is.null(lr)) {
      learning.rate <- keras::k_eval(model$optimizer$lr)      
    } else {
      learning.rate <- lr  
    }
    
    if (!is.null(solver)) {
      solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
    } else {
      warning("You need to specify solver argument if compile is TRUE")
      return(NULL)
    }
    
    if (solver == "adam") {
      optimizer <-  keras::optimizer_adam(lr = learning.rate)
    }
    if (solver == "adagrad") {
      optimizer <- keras::optimizer_adagrad(lr = learning.rate)
    }
    if (solver == "rmsprop") {
      optimizer <- keras::optimizer_rmsprop(lr = learning.rate)
    }
    if (solver == "sgd") {
      optimizer <- keras::optimizer_sgd(lr = learning.rate)
    }
    
    model_new %>% keras::compile(loss = losses,
                                 optimizer = optimizer, metrics = c("acc"))
  }
  
  model_new
}


#' Merge two models
#'
#' @param models List of two models
#' @param layer_names Vector of length 2 with names of layers to merge.
#' @inheritParams create_model_cnn_lstm
#' @export
merge_models <- function(models, layer_names, layer_dense, solver = "adam", learning.rate = 0.0001,
                         freeze_base_model = c(FALSE, FALSE)) {
  
  model_1 <- remove_add_layers(model = models[[1]],
                               layer_name = layer_names[1],
                               layer_depth = NULL,
                               dense_layers = NULL,
                               verbose = FALSE,
                               dropout = NULL,
                               freeze_base_model = freeze_base_model[1],
                               compile = FALSE,
                               lr = NULL)
  
  model_2 <- remove_add_layers(model = models[[2]],
                               layer_name = layer_names[2],
                               layer_depth = NULL,
                               dense_layers = NULL,
                               verbose = FALSE,
                               dropout = NULL,
                               freeze_base_model = freeze_base_model[2],
                               compile = FALSE,
                               lr = NULL)
  
  # choose optimization method
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(lr = learning.rate)
  }
  
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning.rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning.rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning.rate)
  }
  
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
#' @export
check_layer_name <- function(model, layer_name) {
  num_layers <- length(model$get_config()$layers)
  layer_names <- vector("character")
  for (i in 1:num_layers) {
    layer_names[i] <- model$get_config()$layers[[i]]$name
  }
  if (!(layer_name %in% layer_names)) {
    message <- paste0("Model has no layer named ", "'", layer_name, "'")
    stop(message)
  }
}  


#' @title Creates LSTM/CNN network that can multiple samples for one target
#'
#' @inheritParams create_model_lstm_cnn
#' @param samples_per_target Number of samples to combine for one target.
#' @export
create_model_lstm_cnn_time_dist <- function(
  maxlen = 50,
  dropout = 0,
  recurrent_dropout = 0,
  layer_lstm = NULL,
  layer_dense = c(4),
  solver = "adam",
  learning.rate = 0.001,
  use.multiple.gpus = FALSE,
  merge.on.cpu = TRUE,
  gpu.num = 2,
  vocabulary.size = 4,
  bidirectional = FALSE,
  stateful = FALSE,
  batch.size = NULL,
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
  samples_per_target) {
  
  num_output_layers <- 1
  num_input_layers <- 1
  
  num_targets <- layer_dense[length(layer_dense)]
  layers.lstm <- length(layer_lstm)
  use.cnn <- ifelse(!is.null(kernel_size), TRUE, FALSE)
  
  if (!is.null(layer_lstm)) {
    stopifnot(length(layer_lstm) == 1 | (length(layer_lstm) ==  layers.lstm))
  }
  
  if (layers.lstm == 0 & !use.cnn) {
    stop("Model does not use LSTM or CNN layers. Set use.cnn to TRUE or layers.lstm > 0.")
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
  stopifnot(dropout <= 1 & dropout >= 0)
  stopifnot(recurrent_dropout <= 1 & recurrent_dropout >= 0)
  
  if (length(layer_lstm) == 1) {
    layer_lstm <- rep(layer_lstm, layers.lstm)
  }
  
  if (stateful) {
    input_tensor <- keras::layer_input(batch_shape = c(batch.size, maxlen, vocabulary.size))
  } else {
    input_tensor <- keras::layer_input(shape = c(samples_per_target, maxlen, vocabulary.size))
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
            input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
            use_bias = use_bias
          ))
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_max_pooling_1d(pool_size = pool_size[i]))
        }
        output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_batch_normalization(momentum = .8))
      } else {
        
        output_tensor <- output_tensor %>%
          keras::time_distributed(keras::layer_conv_1d(
            kernel_size = kernel_size[i],
            padding = padding,
            activation = "relu",
            filters = filters[i],
            strides = strides[i],
            dilation_rate = dilation_rate[i],
            input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
            use_bias = use_bias
          ))
        output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_batch_normalization(momentum = .8))
        
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::time_distributed(keras::layer_max_pooling_1d(pool_size = pool_size[i]))
        }
        #output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
        
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              keras::layer_lstm(
                units = layer_lstm[i],
                return_sequences = TRUE,
                dropout = dropout,
                recurrent_dropout = recurrent_dropout,
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
              input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
              return_sequences = TRUE,
              dropout = dropout,
              recurrent_dropout = recurrent_dropout,
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
          input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout, recurrent_dropout = recurrent_dropout,
                            stateful = stateful, recurrent_activation = "sigmoid")
        ))
    } else {
      output_tensor <- output_tensor %>%
        keras::time_distributed(keras::layer_lstm(units = layer_lstm[length(layer_lstm)],
                                                  input_shape = switch(stateful + 1, c(maxlen, vocabulary.size), NULL),
                                                  dropout = dropout, recurrent_dropout = recurrent_dropout, stateful = stateful,
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
  
  #output_tensor <- output_tensor %>% keras::layer_add()
  output_tensor <- keras::layer_lambda(output_tensor, f = function(x) {
    tensorflow::tf$math$reduce_sum(x, axis=1L)
  })
  
  if (length(layer_dense) > 1) {
    for (i in 1:(length(layer_dense) - 1)) {
      output_tensor <- output_tensor %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
    }
  }
  
  if (num_output_layers == 1) {
    output_tensor <- output_tensor %>%
      keras::layer_dense(units = num_targets, activation = last_layer_activation)
  } else {
    output_list <- list()
    for (i in 1:num_output_layers) {
      layer_name <- paste0("output_", i, "_", num_output_layers)
      output_list[[i]] <- output_tensor %>%
        keras::layer_dense(units = num_targets, activation = last_layer_activation, name = layer_name)
    }
  }
  
  if (num_output_layers == 1) {
    model <- keras::keras_model(inputs = input_tensor, outputs = output_tensor)
  } else {
    model <- keras::keras_model(inputs = input_tensor, outputs = output_list)
  }
  
  if (use.multiple.gpus) {
    model <- keras::multi_gpu_model(model,
                                    gpus = gpu.num,
                                    cpu_merge = merge.on.cpu)
  }
  
  # choose optimization method
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(lr = learning.rate)
  }
  
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning.rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning.rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning.rate)
  }
  
  #add metrics
  cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  while (dir.exists(cm_dir)) {
    cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  }
  dir.create(cm_dir)
  model$cm_dir <- cm_dir
  metrics <- c("acc")
  
  if (loss_fn == "categorical_crossentropy") {
    
    macro_average_cb <- balanced_acc_wrapper(num_targets, cm_dir)
    metrics <- c(macro_average_cb, "acc")
    
    if (f1_metric) {
      f1 <- f1_wrapper(num_targets)
      metrics <- c(metrics, f1)
    }
  }
  
  if (auc_metric) {
    # multi_label only for tf > 2.2
    multi_label <- ifelse(loss_fn == "binary_crossentropy", TRUE, FALSE)
    #auc <- tensorflow::tf$keras$metrics$AUC(multi_label = multi_label)
    auc <- tensorflow::tf$keras$metrics$AUC()
    metrics <- c(metrics, auc)
  }
  
  model %>% keras::compile(loss = loss_fn,
                           optimizer = optimizer, metrics = metrics)
  
  
  argg <- c(as.list(environment()))
  argg["metrics"] <- NULL
  argg["model"] <- NULL
  argg["i"] <- NULL
  argg["optimizer"] <- NULL
  argg["layer_lstm"] <- paste(as.character(layer_lstm), collapse = " ")
  argg["filters"] <- paste(as.character(filters), collapse = " ")
  argg["kernel_size"] <- paste(as.character(kernel_size), collapse = " ")
  argg["pool_size"] <- paste(as.character(pool_size), collapse = " ")
  argg["strides"] <- paste(as.character(strides), collapse = " ")
  argg["layer_dense"] <- paste(as.character(layer_dense), collapse = " ")
  argg["padding"] <- paste(as.character(padding), collapse = " ")
  argg["use_bias"] <- paste(as.character(use_bias), collapse = " ")
  argg["input_label_list"] <- paste(as.character(layer_dense), collapse = " ")
  argg["input_tensor"] <- NULL
  argg["label_inputs"] <- NULL
  argg["f1"] <- NULL
  argg["multi_acc"] <- NULL
  argg[["trainable_params"]] <- model$count_params()
  argg["output_tensor"] <- NULL
  argg["output_list"] <- NULL
  argg["label_noise_matrix"] <- NULL
  argg["smooth_loss"] <- NULL
  argg["noisy_loss"] <- NULL
  argg["col_sums"] <- NULL
  argg["auc"] <- NULL
  argg["multi_label"] <- NULL
  argg["macro_average_cb"] <- NULL
  argg["input_list"] <- NULL
  argg["metrics"] <- NULL
  argg["representation_list"] <- NULL
  argg["same_length"] <- NULL
  argg["y"] <- NULL
  model$hparam <- argg
  model$cm_dir <- cm_dir
  
  summary(model)
  return(model)
}


#' @title Creates LSTM/CNN network that can process multiple samples for one target
#'
#' @inheritParams create_model_lstm_cnn
#' @param samples_per_target Number of samples to combine for one target.
#' @export
create_model_lstm_cnn_multi_input <- function(
  maxlen = 50,
  dropout = 0,
  recurrent_dropout = 0,
  layer_lstm = NULL,
  layer_dense = c(4),
  solver = "adam",
  learning.rate = 0.001,
  use.multiple.gpus = FALSE,
  merge.on.cpu = TRUE,
  gpu.num = 2,
  vocabulary.size = 4,
  bidirectional = FALSE,
  batch.size = NULL,
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
  num_input_layers = 1,
  auc_metric = FALSE,
  f1_metric = FALSE,
  samples_per_target) {
  
  num_targets <- layer_dense[length(layer_dense)]
  layers.lstm <- length(layer_lstm)
  use.cnn <- ifelse(!is.null(kernel_size), TRUE, FALSE)
  num_output_layers = 1
  
  if (!is.null(layer_lstm)) {
    stopifnot(length(layer_lstm) == 1 | (length(layer_lstm) ==  layers.lstm))
  }
  
  if (layers.lstm == 0 & !use.cnn) {
    stop("Model does not use LSTM or CNN layers. Set use.cnn to TRUE or layers.lstm > 0.")
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
  stopifnot(dropout <= 1 & dropout >= 0)
  stopifnot(recurrent_dropout <= 1 & recurrent_dropout >= 0)
  
  if (length(layer_lstm) == 1) {
    layer_lstm <- rep(layer_lstm, layers.lstm)
  }
  
  input_tensor <- keras::layer_input(shape = c(maxlen, vocabulary.size))
  
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
            input_shape = c(maxlen, vocabulary.size),
            use_bias = use_bias
          )
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
      } else {
        
        output_tensor <- output_tensor %>%
          keras::layer_conv_1d(
            kernel_size = kernel_size[i],
            padding = padding,
            activation = "relu",
            filters = filters[i],
            strides = strides[i],
            dilation_rate = dilation_rate[i],
            input_shape = c(maxlen, vocabulary.size),
            use_bias = use_bias
          )
        output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8) 
        
        if (!is.null(pool_size) && pool_size[i] > 1) {
          output_tensor <- output_tensor %>% keras::layer_max_pooling_1d(pool_size = pool_size[i])
        }
        #output_tensor <- output_tensor %>% keras::layer_batch_normalization(momentum = .8)
        
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
              input_shape = c(maxlen, vocabulary.size),
              keras::layer_lstm(
                units = layer_lstm[i],
                return_sequences = TRUE,
                dropout = dropout,
                recurrent_dropout = recurrent_dropout,
                recurrent_activation = "sigmoid"
              )
            )
        }
      } else {
        for (i in 1:(layers.lstm - 1)) {
          output_tensor <- output_tensor %>%
            keras::layer_lstm(
              units = layer_lstm[i],
              input_shape = c(maxlen, vocabulary.size),
              return_sequences = TRUE,
              dropout = dropout,
              recurrent_dropout = recurrent_dropout,
              recurrent_activation = "sigmoid"
            )
        }
      }
    }
    # last LSTM layer
    if (bidirectional) {
      output_tensor <- output_tensor %>%
        keras::bidirectional(
          input_shape = c(maxlen, vocabulary.size),
          keras::layer_lstm(units = layer_lstm[length(layer_lstm)], dropout = dropout, recurrent_dropout = recurrent_dropout, 
                            recurrent_activation = "sigmoid")
        )
    } else {
      output_tensor <- output_tensor %>%
        keras::layer_lstm(units = layer_lstm[length(layer_lstm)], 
                          input_shape = c(maxlen, vocabulary.size),
                          dropout = dropout, recurrent_dropout = recurrent_dropout,
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
    input_list[[i]] <- keras::layer_input(shape = c(maxlen, vocabulary.size), name = paste0("input_", i))
    representation_list[[i]] <- feature_ext_model(input_list[[i]])
  }
  
  y <- keras::layer_add(representation_list) 
  #y <- keras::layer_average(representation_list) 
  
  if (length(layer_dense) > 1) {
    for (i in 1:(length(layer_dense) - 1)) {
      y <- y %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
    }
  }
  
  y <- y %>% keras::layer_dense(units = num_targets, activation = last_layer_activation) 
  
  # print model layout to screen, should be done before multi_gpu_model
  
  model <- keras::keras_model(inputs = input_list, outputs = y)
  
  if (use.multiple.gpus) {
    model <- keras::multi_gpu_model(model,
                                    gpus = gpu.num,
                                    cpu_merge = merge.on.cpu)
  }
  
  # choose optimization method
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(lr = learning.rate)
  }
  
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning.rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning.rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning.rate)
  }
  
  #add metrics
  cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  while (dir.exists(cm_dir)) {
    cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
  }
  dir.create(cm_dir)
  model$cm_dir <- cm_dir
  metrics <- c("acc")
  
  if (loss_fn == "categorical_crossentropy") {
    
    macro_average_cb <- balanced_acc_wrapper(num_targets, cm_dir)
    metrics <- c(macro_average_cb, "acc")
    
    if (f1_metric) {
      f1 <- f1_wrapper(num_targets)
      metrics <- c(metrics, f1)
    }
  }
  
  if (auc_metric) {
    # multi_label only for tf > 2.2
    multi_label <- ifelse(loss_fn == "binary_crossentropy", TRUE, FALSE)
    #auc <- tensorflow::tf$keras$metrics$AUC(multi_label = multi_label)
    auc <- tensorflow::tf$keras$metrics$AUC()
    metrics <- c(metrics, auc)
  }
  
  model %>% keras::compile(loss = loss_fn,
                           optimizer = optimizer, metrics = metrics)
  
  
  argg <- c(as.list(environment()))
  argg["metrics"] <- NULL
  argg["model"] <- NULL
  argg["i"] <- NULL
  argg["optimizer"] <- NULL
  argg["layer_lstm"] <- paste(as.character(layer_lstm), collapse = " ")
  argg["filters"] <- paste(as.character(filters), collapse = " ")
  argg["kernel_size"] <- paste(as.character(kernel_size), collapse = " ")
  argg["pool_size"] <- paste(as.character(pool_size), collapse = " ")
  argg["strides"] <- paste(as.character(strides), collapse = " ")
  argg["layer_dense"] <- paste(as.character(layer_dense), collapse = " ")
  argg["padding"] <- paste(as.character(padding), collapse = " ")
  argg["use_bias"] <- paste(as.character(use_bias), collapse = " ")
  argg["input_label_list"] <- paste(as.character(layer_dense), collapse = " ")
  argg["input_tensor"] <- NULL
  argg["label_inputs"] <- NULL
  argg["f1"] <- NULL
  argg["multi_acc"] <- NULL
  argg[["trainable_params"]] <- model$count_params()
  argg["output_tensor"] <- NULL
  argg["output_list"] <- NULL
  argg["label_noise_matrix"] <- NULL
  argg["smooth_loss"] <- NULL
  argg["noisy_loss"] <- NULL
  argg["col_sums"] <- NULL
  argg["auc"] <- NULL
  argg["multi_label"] <- NULL
  argg["macro_average_cb"] <- NULL
  argg["input_list"] <- NULL
  argg["metrics"] <- NULL
  argg["representation_list"] <- NULL
  argg["same_length"] <- NULL
  argg["y"] <- NULL
  argg["feature_ext_model"] <- NULL
  model$hparam <- argg
  model$cm_dir <- cm_dir
  
  summary(model)
  return(model)
}

#' Replace input layer
#' 
#' Replace first layer of model with new input layer of different shape. Only works for sequential models.
#' 
#' @param model A keras model.
#' @param input_shape The new input shape vector (without batch size). 
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
  return(new_model)
}
