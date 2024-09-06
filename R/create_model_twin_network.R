#' @title Create twin network 
#'
#' @description Twin network can be trained to maximize the distance
#' between embeddings of inputs.
#' Implements approach as described [here](https://keras.io/examples/vision/siamese_contrastive/).
#'
#' @param layer_dense Vector containing number of neurons per dense layer, before euclidean distance layer.
#' @param distance_method Either "euclidean" or "cosine".
#' @param metrics Vector or list of metrics.
#' @inheritParams create_model_lstm_cnn
#' @inheritParams create_model_lstm_cnn_multi_input
#' @examplesIf reticulate::py_module_available("tensorflow")
#' model <- create_model_twin_network(
#'   maxlen = 50,
#'   layer_dense = 16,
#'   kernel_size = 12,
#'   filters = 4,
#'   pool_size = 3,
#'   learning_rate = 0.001)
#'   
#' @returns A keras model implementing twin network architecture.   
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
    pool_method = 'max',
    solver = "adam",
    learning_rate = 0.001,
    vocabulary_size = 4,
    bidirectional = FALSE,
    compile = TRUE,
    padding = "same",
    dilation_rate = NULL,
    gap_inputs = NULL,
    use_bias = TRUE,
    residual_block = FALSE,
    residual_block_length = 1,
    size_reduction_1Dconv = FALSE,
    zero_mask = FALSE,
    verbose = TRUE,
    batch_norm_momentum = 0.99,
    distance_method = "euclidean",
    last_layer_activation = "sigmoid",
    loss_fn = loss_cl(margin=1),
    metrics = "acc",
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
    pool_method == pool_method,
    padding = padding,
    dilation_rate = dilation_rate,
    gap_inputs = gap_inputs,
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
    model %>% keras::compile(loss = loss_fn,
                             optimizer = set_optimizer(solver, learning_rate),
                             metrics = metrics)
  }
  
  if (verbose) model$summary()
  model
}

