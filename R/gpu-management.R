#' Limit Tensorflow proccess to one or multiple GPUs.
#'
#' @param gpus GPU ids
#' @param growth memory consumption of GPU grows and will not be blocked
#' @export
startGPUSession <- function(gpus = "0", growth = T){
  
  require(tensorflow)
  tf$reset_default_graph()
  sess_config <- list()
  Sys.setenv(CUDA_VISIBLE_DEVICES = gpus)
  sess_config$device_count <- list(GPU = 1L, CPU = 1L)
  if (growth)
    sess_config$gpu_options <- tf$GPUOptions(allow_growth = TRUE)
  session_conf <- do.call(tf$ConfigProto, sess_config)
  sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)
  sess$run(tf$global_variables_initializer())
  return(sess)
}

#' Ends the GPU session
#'
#' @export
end_gpu_session <- function(){
  sess$close()
}

list.local.devices <- function(){
  message(tensorflow::tf$python$client$device_lib$list_local_devices())
}

is.gpu.available <- function() {
  res <- tryCatch({
    tensorflow::tf$test$is_gpu_available()
  }, error = function(e) {
    warning("Can not determine if GPU is configured.", call. = FALSE);
    NA
  })
  res
}

is.cuda.build <- function() {
  res <- tryCatch({
    tensorflow::tf$test$is_built_with_cuda()
  }, error = function(e) {
    warning("Can not determine if TF is build with CUDA", call. = FALSE);
    NA
  })
  res
}
