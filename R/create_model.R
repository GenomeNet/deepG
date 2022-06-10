#' @title Creates LSTM/CNN network
#'
#' @description Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers.
#' Last layer is a dense layer with default softmax activation.
#' @param maxlen Length of predictor sequence.
#' @param dropout_lstm Fraction of the units to drop for inputs.
#' @param recurrent_dropout_lstm Fraction of the units to drop for recurrent state.
#' @param layer_lstm Number of cells per network layer. Can be a scalar or vector.
#' @param layer_dense Dense layers of size layer_dense after last LSTM (or last CNN is \code{layers.lstm = 0}) layer.
#' @param solver Optimization method, options are "adam", "adagrad", "rmsprop" or "sgd".
#' @param learning_rate Learning rate for optimizer.
#' @param use_multiple_gpus If true, multi_gpu_model() will be used based on gpu_num.
#' @param gpu_num Number of GPUs to be used, only relevant if multiple_gpu is true.
#' @param merge_on_cpu True on default, false recommend if the server supports NVlink, only relevant if use.multiple.gpu is true.
#' @param bidirectional Use bidirectional wrapper for lstm layers.
#' @param vocabulary_size Number of unique character in vocabulary.
#' @param stateful Boolean. Whether to use stateful LSTM layer.
#' @param batch_size Number of samples that are used for one network update. Only used if \code{stateful = TRUE}.
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
#' If first label contains 5 percent wrong labels and second label no noise, then
#' \code{label_noise_matrix <- matrix(c(0.95, 0.05, 0, 1), nrow = 2, byrow = TRUE )}
#' @param last_layer_activation Either "sigmoid" or "softmax".
#' @param loss_fn Either "categorical_crossentropy" or "binary_crossentropy". If label_noise_matrix given, will use custom "noisy_loss".
#' @param num_output_layers Number of output layers.
#' @param auc_metric Whether to add AUC metric.
#' @param f1_metric Whether to add F1 metric.
#' @param bal_acc Whether to add balanced accuracy.
#' @param batch_norm_momentum Momentum for the moving mean and the moving variance.
#' @export
create_model_lstm_cnn <- function(
  maxlen = 50,
  dropout_lstm = 0,
  recurrent_dropout_lstm = 0,
  layer_lstm = NULL,
  layer_dense = c(4),
  solver = "adam",
  learning_rate = 0.001,
  use_multiple_gpus = FALSE,
  merge_on_cpu = TRUE,
  gpu_num = 2,
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
  batch_norm_momentum = 0.99) {

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

  if (use_multiple_gpus) {
    model <- keras::multi_gpu_model(model,
                                    gpus = gpu_num,
                                    cpu_merge = merge_on_cpu)
  }

  # choose optimization method
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(lr = learning_rate)
  }
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning_rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning_rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning_rate)
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
      auc <- auc_wrapper(model_output_size = layer_dense[length(layer_dense)],
                         loss = loss_fn)
      metrics <- c(model$metrics, auc)
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

  if (verbose) summary(model)
  return(model)
}

#' create wavenet model
#'
#' @inheritParams wavenet::wavenet
#' @export
create_model_wavenet <- function(filters = 16, kernel_size = 2, residual_blocks, maxlen,
                                 input_tensor = NULL, initial_kernel_size = 32, initial_filters = 32,
                                 output_channels = 4, output_activation = "softmax", solver = "adam",
                                 learning_rate = 0.001, compile = TRUE, verbose = TRUE) {

  model <- wavenet::wavenet(filters = filters, kernel_size = kernel_size, residual_blocks = residual_blocks,
                            input_shape = list(maxlen, output_channels), input_tensor = input_tensor, initial_kernel_size = initial_kernel_size,
                            initial_filters = initial_filters, output_channels = output_channels, output_activation = "softmax")
  if (solver == "adam") {
    optimizer <- keras::optimizer_adam(lr = learning_rate)
  }
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning_rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning_rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning_rate)
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
  if (verbose) summary(model)
  model
}

#######

#' @title Creates LSTM/CNN network
#'
#' @description
#' Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers.
#' Last layer is a dense layer with softmax activation.
#' Network tries to predict target in the middle of a sequence. If input is AACCTAAGG, input tensors should correspond to x1 = AACC, x2 = GGAA and y = T
#' (set \code{target_middle = TRUE} in \code{train_model} function).
#' Function creates two sub networks consisting each of an (optional) CNN layer followed by an arbitrary number of LSTM layers. Afterwards the last LSTM layers
#' get concatenated and followed by a dense layers.
#' @param maxlen Length of predictor sequence.
#' @param dropout_lstm Fraction of the units to drop for inputs.
#' @param recurrent_dropout_lstm Fraction of the units to drop for recurrent state.
#' @param layer_lstm Number of cells per network layer. Can be a scalar or vector.
#' @param solver Optimization method, options are "adam", "adagrad", "rmsprop" or "sgd".
#' @param learning_rate Learning rate for optimizer.
#' @param use_multiple_gpus If true, multi_gpu_model() will be used based on gpu_num.
#' @param gpu_num Number of GPUs to be used, only relevant if multiple_gpu is true.
#' @param merge_on_cpu True on default, false recommend if the server supports NVlink, only relevant if use.multiple.gpu is true.
#' @param bidirectional Use bidirectional wrapper for lstm layers.
#' @param vocabulary_size Number of unique character in vocabulary.
#' @param stateful Boolean. Whether to use stateful LSTM layer.
#' @param batch_size Number of samples that are used for one network update. Only used if \code{stateful = TRUE}.
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
  dropout_lstm = 0,
  recurrent_dropout_lstm = 0,
  layer_lstm = 128,
  solver = "adam",
  learning_rate = 0.001,
  use_multiple_gpus = FALSE,
  merge_on_cpu = TRUE,
  gpu_num = 2,
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
  batch_norm_momentum = 0.99) {

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

  if (use_multiple_gpus) {
    model <- keras::multi_gpu_model(model,
                                    gpus = gpu_num,
                                    cpu_merge = merge_on_cpu)
  }

  # choose optimization method
  if (solver == "adam")
    optimizer <-
    keras::optimizer_adam(lr = learning_rate)
  if (solver == "adagrad")
    optimizer <-
    keras::optimizer_adagrad(lr = learning_rate)
  if (solver == "rmsprop")
    optimizer <-
    keras::optimizer_rmsprop(lr = learning_rate)
  if (solver == "sgd")
    optimizer <-
    keras::optimizer_sgd(lr = learning_rate)

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
    model %>% keras::compile(loss = "categorical_crossentropy",
                             optimizer = optimizer, metrics = metrics)
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

  if (verbose) summary(model)
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
#' @param layer_name Name of last layer to use from old model.
#' @param model A keras model. If model and path_model are both not NULL, path_model will be used.
#' @param dense_layers List of vectors specifiying number of units for each dense layer. If this is a list of length > 1, model
#' has multiple output layers.
#' @param last_activation List of activations for last entry for each list entry from \code{dense_layers}. Either "softmax", "sigmoid" or "linear".
#' @param output_names List of names for each output layer.
#' @param losses List of loss function for each output.
#' @param verbose Boolean.
#' @param dropout List of vectors with dropout rates for each new dense layer.
#' @param freeze_base_model Whether to freeze all weights before new dense layers.
#' @param compile Boolean, whether the new model is compiled or not
#' @param lr learning_rate if compile == TRUE, default -> learning_rate of the old model
#' @param solver "adam", "adagrad", "rmsprop" or "sgd" if compile == TRUE, default -> solver of the old model
#' @export
remove_add_layers <- function(model = NULL,
                              layer_name = NULL,
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
      learning_rate <- keras::k_eval(model$optimizer$lr)
    } else {
      learning_rate <- lr
    }

    if (is.null(solver)) {
      solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
    } #else {
    #   warning("You need to specify solver argument if compile is TRUE")
    #   return(NULL)
    # }

    if (solver == "adam") {
      optimizer <-  keras::optimizer_adam(lr = learning_rate)
    }
    if (solver == "adagrad") {
      optimizer <- keras::optimizer_adagrad(lr = learning_rate)
    }
    if (solver == "rmsprop") {
      optimizer <- keras::optimizer_rmsprop(lr = learning_rate)
    }
    if (solver == "sgd") {
      optimizer <- keras::optimizer_sgd(lr = learning_rate)
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
#' @inheritParams create_model_lstm_cnn
#' @export
merge_models <- function(models, layer_names, layer_dense, solver = "adam", learning_rate = 0.0001,
                         freeze_base_model = c(FALSE, FALSE)) {

  model_1 <- remove_add_layers(model = models[[1]],
                               layer_name = layer_names[1],
                               dense_layers = NULL,
                               verbose = FALSE,
                               dropout = NULL,
                               freeze_base_model = freeze_base_model[1],
                               compile = FALSE,
                               lr = NULL)

  model_2 <- remove_add_layers(model = models[[2]],
                               layer_name = layer_names[2],
                               dense_layers = NULL,
                               verbose = FALSE,
                               dropout = NULL,
                               freeze_base_model = freeze_base_model[2],
                               compile = FALSE,
                               lr = NULL)

  # choose optimization method
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(lr = learning_rate)
  }

  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning_rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning_rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning_rate)
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
  dropout_lstm = 0,
  recurrent_dropout_lstm = 0,
  layer_lstm = NULL,
  layer_dense = c(4),
  solver = "adam",
  learning_rate = 0.001,
  use_multiple_gpus = FALSE,
  merge_on_cpu = TRUE,
  gpu_num = 2,
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
  verbose = TRUE) {

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

  if (use_multiple_gpus) {
    model <- keras::multi_gpu_model(model,
                                    gpus = gpu_num,
                                    cpu_merge = merge_on_cpu)
  }

  # choose optimization method
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(lr = learning_rate)
  }

  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning_rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning_rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning_rate)
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
    auc <- auc_wrapper(model_output_size = layer_dense[length(layer_dense)],
                       loss = loss_fn)
    metrics <- c(model$metrics, auc)
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

  if (verbose) summary(model)
  return(model)
}


#' @title Creates LSTM/CNN network that can process multiple samples for one target
#'
#' @inheritParams create_model_lstm_cnn
#' @param samples_per_target Number of samples to combine for one target.
#' @export
create_model_lstm_cnn_multi_input <- function(
  maxlen = 50,
  dropout_lstm = 0,
  recurrent_dropout_lstm = 0,
  layer_lstm = NULL,
  layer_dense = c(4),
  solver = "adam",
  learning_rate = 0.001,
  use_multiple_gpus = FALSE,
  merge_on_cpu = TRUE,
  gpu_num = 2,
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
  num_input_layers = 1,
  auc_metric = FALSE,
  f1_metric = FALSE,
  samples_per_target,
  batch_norm_momentum = 0.99,
  verbose = TRUE) {

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
      y <- y %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
    }
  }

  y <- y %>% keras::layer_dense(units = num_targets, activation = last_layer_activation)

  # print model layout to screen, should be done before multi_gpu_model

  model <- keras::keras_model(inputs = input_list, outputs = y)

  if (use_multiple_gpus) {
    model <- keras::multi_gpu_model(model,
                                    gpus = gpu_num,
                                    cpu_merge = merge_on_cpu)
  }

  # choose optimization method
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(lr = learning_rate)
  }

  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning_rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning_rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning_rate)
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
    auc <- auc_wrapper(model_output_size = layer_dense[length(layer_dense)],
                       loss = loss_fn)
    metrics <- c(model$metrics, auc)
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

  if (verbose) summary(model)
  return(model)
}

#' Replace input layer
#'
#' Replace first layer of model with new input layer of different shape. Only works for sequential models that
#' uses CNN and LSTM layers.
#'
#' @param model A keras model.
#' @param input_shape The new input shape vector (without batch_size).
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


#' Get solver and learning_rate from model.
#'
get_optimizer <- function(model) {
  solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
  learning_rate <- keras::k_eval(model$optimizer$lr)
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(lr = learning_rate)
  }
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(lr = learning_rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(lr = learning_rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(lr = learning_rate)
  }
  return(optimizer)
}

#' Self-genomenet model 
#'
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
  num_targets = 2) {

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
    keras::layer_dense(units = num_targets, activation = "softmax")

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

  keras_optimizer <- switch(optimizer,
                            adam = keras::optimizer_adam(lr = learning_rate),
                            adagrad = keras::optimizer_adagrad(lr = learning_rate),
                            rmsprop = keras::optimizer_rmsprop(lr = learning_rate),
                            sgd = keras::optimizer_sgd(lr = learning_rate)
  )

  model %>% keras::compile(loss = "categorical_crossentropy", optimizer = keras_optimizer, metrics = "acc")

  argg <- c(as.list(environment()))
  argg["metrics"] <- NULL
  argg["model"] <- NULL
  argg["i"] <- NULL
  argg["optimizer"] <- NULL
  argg["number_of_cnn_layers"] <- paste(as.character(number_of_cnn_layers), collapse = " ")
  argg["filters"] <- paste(as.character(filters), collapse = " ")
  argg["kernel_size"] <- paste(as.character(kernel_size), collapse = " ")
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

  model
}

#' Load pretrained self-genomenet model
#'
#' Classification model with labels "bacteria", "virus-no-phage","virus-phage".
#' TODO: add link to paper
#'
#' @inheritParams create_model_lstm_cnn
#' @param maxlen Model input size. Either 150 or 10000.
#' @param learning_rate Learning rate for optimizer. If compile is TRUE and learning_rate is NULL,
#' will use learning rate from previous training.
#' @export
load_model_self_genomenet <- function(maxlen, compile = FALSE, optimizer = "adam",
                                      learning_rate = NULL) {

  stopifnot(any(maxlen == c(150,10000)))

  if (maxlen == 150) {
    data(model_self_genomenet_maxlen_150)
    model <- keras::unserialize_model(model_self_genomenet_maxlen_150, compile = FALSE)
  }

  if (maxlen == 10000) {
    data(model_self_genomenet_maxlen_10k)
    model <- keras::unserialize_model(model_self_genomenet_maxlen_10k, compile = FALSE)
  }

  if (is.null(learning_rate)) {
    if (maxlen == 150) learning_rate <- 0.00039517784549691
    if (maxlen == 10000) learning_rate <- 8.77530464905713e-05
  }

  if (compile) {
    keras_optimizer <- switch(optimizer,
                              adam = keras::optimizer_adam(lr = learning_rate),
                              adagrad = keras::optimizer_adagrad(lr = learning_rate),
                              rmsprop = keras::optimizer_rmsprop(lr = learning_rate),
                              sgd = keras::optimizer_sgd(lr = learning_rate)
    )

    model %>% keras::compile(loss = "categorical_crossentropy", optimizer = keras_optimizer, metrics = "acc")
  }

  return(model)
}
