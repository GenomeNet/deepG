#' deepG for GenomeNet
#'
#' deepG is a is an open source software library for building deep neuronal 
#' networks for genomic modeling
#' 
#' This package generates \href{http://www.genomenet.de}{GenomeNet}
#'
#' For additional documentation on the deepG package see
#' \href{https://genomenet.de}{https://genomenet.de}
#'
#' @import reticulate
#' @import keras
#' @import tensorflow
#' @import dplyr
#'
#' @docType package
#' @name deepG
NULL

# globals
.globals <- new.env(parent = emptyenv())
.globals$tensorboard <- NULL

.onLoad <- function(libname, pkgname) {  
  # call hdf5r function to avoid error message 
  temp_file <- tempfile(fileext = ".h5")
  h5_file <- hdf5r::H5File$new(temp_file, mode = "w") 
  h5_file$close_all()
  file.remove(temp_file)

  #usethis::use_pipe(export = TRUE)
  
  packageStartupMessage("The deepG package has been successfully loaded.")
}
