#' generates hdf5 file containing the character id for each time step
#'
#' @param dat character
#' @param filename filename where hdf5 file is written to
#' @param verbose TRUE/FALSE
#' @export
writeHdf5 <- function(dat, 
                      filename = "train.hdf5", 
                      verbose = F) {

	tokenized <- tokenizers::tokenize_characters(dat, strip_non_alphanum = FALSE,
																		simplify = TRUE)
	# get sorted vocabulary
	charset <- sort(unique(unlist(dat_tokenized)))

	# get corresponding character numbers
	charset_index <- 1:length(charset)
	names(charset_index) <- charset

	# replace character by character index
	tokenized_index <- plyr::mapvalues(tokenized, from = charset,
															 to = charset_index)

	if (verbose) message("Saving states...")
	file <- hdf5r::H5File$new(filename, mode = "a")
	file.grp <- hdf5r::file.h5$create_group("words")
	file.grp <- array(tokenized_index)
	hdf5r::h5close(file)
}


#' Generates dictionary for LSTMVis
#'
#' @param dat filename where dict file is written to
#' @param filename file name
#' @export
writeDict <- function(dat, filename = "train.dict") {

	tokenized <- tokenizers::tokenize_characters(dat, strip_non_alphanum = FALSE,
																		simplify = TRUE)
	# get sorted vocabulary
	charset <- sort(unique(unlist(dat_tokenized)))

	dict <- data.frame(char = charset,
										 index = 1:length(charset))
	write.table(dict, file = filename, quote = F, row.names = F,
							col.names = F, sep = " ")
}

#' Wrapper for message(sprintf)
#' @param ... text input
#' @param newline print in new line 
#' @export
messagef <- function (..., newline = TRUE) 
{
  message(sprintf(...), appendLF = newline)
}


print.tf.version <- function(x) 
{
  message(paste("Tensorflow", tensorflow::tf$`__version__`, "found."))
}

#' Calculate the steps per epoch
#' 
#' Do one full preprocessing iteration to the FASTA file to figure out what the observed
#' steps_per_epoch value is.
#' @param dir Input directory where .fasta files are located
#' @param batch.size Number of samples  
#' @param maxlen Length of the semi-redundant sequences
#' @param format File format
#' @param verbose TRUE/FALSE
#' @export
calculateStepsPerEpoch <-
  function(dir,
           batch.size = 256,
           maxlen = 250,
           format = "fasta",
           verbose = F) {

    steps.per.epoch <- 0
    fasta.files <- list.files(
      path = xfun::normalize_path(dir),
      pattern = paste0("*.", format),
      full.names = TRUE
    )
    for (file in fasta.files) {
      fasta.file <- Biostrings::readDNAStringSet(file)
      seq <- paste0(paste(fasta.file, collapse = "\n"), "\n")
      steps.per.epoch <-
        steps.per.epoch + ceiling((nchar(seq) - maxlen - 1) / batch.size)
    }
    return(steps.per.epoch)
  }


#' Tests if Keras is available for the unit tests
#' @param version required keras version 
skip_if_no_keras <- function(version = NULL){
  if(!is_keras_available(version))
    skip("Required keras version not avaible for testing!")
}
