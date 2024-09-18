#' deepG for GenomeNet
#'
#' deepG is an open source software library for building deep neural
#' networks for genomic modeling.
#'
#' This package generates \href{https://deepg.de}{deepG}
#'
#' For additional documentation on the deepG package see
#' \href{https://genomenet.de}{https://genomenet.de}
#'
#' @import reticulate
#' @import keras
#' @import tensorflow
#' @import dplyr
#'
"_PACKAGE"

# globals
.globals <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  
  # Minimal HDF5 operation to prevent errors from hdf5r package
  temp_file <- tempfile(fileext = ".h5")
  h5_file <- hdf5r::H5File$new(temp_file, mode = "w")
  h5_file$close_all()
  file.remove(temp_file)
  
  # # Set TensorFlow log level to suppress excessive logs
  # if (!exists("tf_initialized", envir = .globals)) {
  #   Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3)
  #   if (reticulate::py_module_available("tensorflow")) {
  #     .globals$tf_initialized <- TRUE
  #     tensorflow::tf$constant(1)
  #   } else {
  #     .globals$tf_initialized <- FALSE
  #   }
  # }
  
  
  # if (!exists("tf_initialized", envir = .globals)) {
  #   Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3)
  # 
  #   tryCatch({
  #     if (requireNamespace('tensorflow')) {
  #       # Perform a minimal TensorFlow operation to ensure it's working
  #       tensorflow::tf$constant(1)
  #       .globals$tf_initialized <- TRUE  # Mark TensorFlow as initialized
  #     }
  # 
  #   }, error = function(e) {
  #     .globals$tf_initialized <- FALSE
  #   })
  # 
  # }
  
  return(NULL)
  
}

.onAttach <- function(libname, pkgname) {
  # Check if TensorFlow is available globally and store the result
  .globals$tf_available <- reticulate::py_module_available("tensorflow")
  
  if (!.globals$tf_available) {
    packageStartupMessage("TensorFlow is not available. Some examples will be skipped.")
  }
}

