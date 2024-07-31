#' @title Create GenomeNet Model with Given Architecture Parameters
#'
#' @param maxlen (integer `numeric(1)`)\cr
#'   Input sequence length.
#' @param learning_rate (`numeric(1)`)\cr
#'   Used by the `keras` optimizer that is specified by `optimizer`.
#' @param number_of_cnn_layers (integer `numeric(1)`)\cr
#'   Target number of CNN-layers to use in total. If `number_of_cnn_layers` is
#'   greater than `conv_block_count`, then the effective number of CNN layers
#'   is set to the closest integer that is divisible by `conv_block_count`.
#' @param conv_block_count (integer `numeric(1)`)\cr
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
#' @param kernel_size_0 (`numeric(1)`)\cr
#'   Target CNN kernel size of the first CNN-layer. Although CNN kernel size is
#'   always an integer, this value can be non-integer, potentially affecting
#'   the kernel-sizes of intermediate layers (which are geometrically
#'   interpolated between `kernel_size_0` and `kernel_size_end`).
#' @param kernel_size_end (`numeric(1)`)\cr
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
#'   the filter-numbers of intermediate dilation_rates layers (which are geometrically
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
#' @param vocabulary_size (integer `numeric(1)`)\cr
#'   Vocabulary size of (one-hot encoded) input strings. This determines the
#'   input tensor shape, together with `maxlen`.
#' @param last_layer_activation Either `"sigmoid"` or `"softmax"`.
#' @param loss_fn Either `"categorical_crossentropy"` or `"binary_crossentropy"`. If `label_noise_matrix` given, will use custom `"noisy_loss"`.
#' @param auc_metric Whether to add AUC metric.
#' @param num_targets (integer `numeric(1)`)\cr
#'   Number of output units to create.
#' @return A keras model.
#' @inheritParams create_model_lstm_cnn
#' @examplesIf reticulate::py_module_available("tensorflow")
#' model <- create_model_genomenet()
#' model
#' 
#' @returns A keras model implementing genomenet architecture.
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
    dropout = 0,
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
    auc_metric = FALSE,
    num_targets = 2,
    model_seed = NULL,
    bal_acc = FALSE,
    f1_metric = FALSE,
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
    
    output_collection <- utils::tail(output_collection, use_blocks)
    
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
    layer <- keras::layer_dropout(rate = dropout)
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
  
  model %>% keras::compile(loss = loss_fn, optimizer = keras_optimizer, metrics = model_metrics)
  
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
#   return(model)
# }
