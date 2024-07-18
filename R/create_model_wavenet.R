#' Create wavenet model
#'
#' Creates network architecture as described [here](https://arxiv.org/abs/1609.03499). Implementation 
#' uses code from [here](https://github.com/r-tensorflow/wavenet).
#' 
#' @inheritParams wavenet::wavenet
#' @inheritParams create_model_lstm_cnn
#' @examples 
#'
#' model <- create_model_wavenet(residual_blocks = 2^rep(1:4, 2), maxlen = 1000)
#'  
#' @returns A keras model implementing wavenet architecture.
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
  if (verbose) model$summary()
  model
}
