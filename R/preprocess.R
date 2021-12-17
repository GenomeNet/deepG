#' Returns the vocabulary from character string
#'
#' Use this function with a character string.
#'
#' @param char character string of text with the length of one
#' @param verbose TRUE/FALSE
#' @examples
#' getVocabulary(data(crispr_sample))
#' getVocabulary("abcd")
#' @export
getVocabulary <- function(char, verbose = F) {

  stopifnot(!is.null(char))
  stopifnot(nchar(char) > 0)

  vocabulary <- sort(unique(tokenizers::tokenize_characters(
    stringr::str_c(stringr::str_to_lower(char),collapse = "\n"), strip_non_alphanum = FALSE, simplify = TRUE)))

  if (verbose)
    message("The vocabulary:", vocabulary)
  return(vocabulary)
}

#' Preprocess string to semi-redundant one-hot vector
#'
#' @description
#' Outputs semi-redundant set of input character string.
#' Collapse, tokenize, and vectorize the character.
#' Use this function with a character string as input. For example,
#' if the input text is ABCDEFGHI and the length(maxlen) is 5, the generating chunks would be:
#' X(1): ABCDE and Y(1): F;
#' X(2): BCDEF and Y(2): G;
#' X(3): CDEFG and Y(3): H;
#' X(4): DEFGH and Y(4): I
#'
#' @param char character input string of text with the length of one
#' @param maxlen length of the semi-redundant sequences
#' @param vocabulary char contains the vocabulary from the input char
#' If no vocabulary exists, it is generated from the input char
#' @param verbose TRUE/FALSE
#' @export
preprocessSemiRedundant <- function(char,
                                    maxlen = 250,
                                    vocabulary = c("l", "p", "a", "c", "g", "t"),
                                    verbose = F) {

  stopifnot(!is.null(char))
  stopifnot(nchar(char) > 0)
  stopifnot(maxlen > 0)

  # Load, collapse, and tokenize text ("ACGT" -> "a" "c" "g" "t")
  text <- tokenizers::tokenize_characters(stringr::str_c(stringr::str_to_lower(char), collapse = "\n"), strip_non_alphanum = FALSE, simplify = TRUE)

  # Generating vocabulary from input char with the function getVocabulary()
  if (missing(vocabulary)) {
    if (verbose)
      message("Finding the vocabulary ...")
    vocabulary <- getVocabulary(char)
  }

  if(verbose)
    message("Vocabulary size:", length(vocabulary))
  # Cut the text in semi-redundant sequences of maxlen characters

  if (verbose)
    message("Generation of semi-redundant sequences ...")

  dataset <- purrr::map(seq(1, length(text) - maxlen, by = 1),
                        ~ list(sentece = text[.x:(.x + maxlen - 1)],
                               next_char = text[.x + maxlen]))
  dataset <- purrr::transpose(dataset)
  x <-
    array(0, dim = c(length(dataset$sentece), maxlen, length(vocabulary)))
  y <- array(0, dim = c(length(dataset$sentece), length(vocabulary)))
  # Vectorization

  if (verbose)
    message("Vectorization ...")
  if (verbose)
    pb <-  txtProgressBar(min = 0,
                          max = length(dataset$sentece),
                          style = 3)
  for (i in 1:length(dataset$sentece)) {
    if (verbose)
      setTxtProgressBar(pb, i)
    # generate one-hot encoding for one subset
    x[i, ,] <- sapply(vocabulary, function(x) {
      as.integer(x == dataset$sentece[[i]])
    })
    # target (next nucleotide in sequence)
    y[i,] <- as.integer(vocabulary == dataset$next_char[[i]])
  }

  results <- list("X" = x, "Y" = y)
  return(results)
}

#' Wrapper of the preprocessSemiRedundant()-function
#'
#' @description
#' Is called on the genomic contents of one
#' FASTA file. Multiple entries are combined with newline characters.
#' @param path path to the FASTA file
#' @param maxlen length of the semi-redundant sequences
#' @param vocabulary char contains the vocabulary from the input char
#' If no vocabulary exists, it is generated from the input char
#' @param verbose TRUE/FALSE
#' @export
preprocessFasta <- function(path,
                            maxlen = 250,
                            vocabulary = c("l", "p", "a", "c", "g", "t"),
                            verbose = F) {


  # process corpus
  if (endsWith(path, "fasta")) {
    fasta.file <- microseq::readFasta(path)
  }
  if (endsWith(path, "fastq")) {
    fasta.file <- microseq::readFastq(path)
  }
  seq <- paste(fasta.file$Sequence, collapse = "")

  if(verbose)
    message("Preprocessing the data ...")

  seq.processed <-
    preprocessSemiRedundant(char = seq, maxlen = maxlen, vocabulary = vocabulary,
                            verbose = F)
  return(seq.processed)
}

#' One-hot-encodes integer sequence
#'
#' \code{sequenceToArray} Helper function for \code{\link{{fastaFileGenerator}}, returns one hot encoding for sequence
#'
#' @param sequence Sequence of integers.
#' @param maxlen Length of one sample
#' @param vocabulary Set of characters to encode.
#' @param startInd Start positions of samples in \code{sequence}.
#' @param wavenet_format Boolean.
#' @param cnn_format Boolean. If true, nucleotides on the left and right side of the predicted nucleotide (language model) are
#' concatenated in the first layer, which make the structure more suitable to CNN architectures (without RNN).
#' @param target_middle Boolean, target is in middle of sequence.
#' @param ambiguous_nuc How to handle nucleotides outside vocabulary, either "zero", "discard" or "equal". If "zero", input gets encoded as zero vector;
#' if "equal" input is 1/length(vocabulary) x length(vocabulary). If "discard" samples containing nucleotides outside vocabulary get discarded.
#' @param nuc_dist Nucleotide distribution.
#' @param use_quality Use quality scores.
#' @param quality_vector Vector of quality probabilities.
#' @export
sequenceToArray <- function(sequence, maxlen, vocabulary, startInd, wavenet_format = FALSE, target_middle = FALSE,
                            ambiguous_nuc = "zero", nuc_dist = NULL, use_quality = FALSE, quality_vector = NULL,
                            cnn_format, target_len = 1, use_coverage = FALSE, max_cov = NULL, cov_vector = NULL,
                            n_gram = NULL, n_gram_stride = 1) {

  voc_len <- length(vocabulary)
  if (!is.null(n_gram)) {
    if (target_len < n_gram) stop("target_len needs to be at least as big as n_gram")
    #   sequence <- int_to_n_gram(int_seq = sequence, n = n_gram, voc_size = length(vocabulary))
    #   target_len <- target_len - n_gram + 1
    #   maxlen <- maxlen - n_gram + 1
    #   voc_len <- length(vocabulary)^n_gram
  }

  startInd <- startInd - startInd[1] + 1
  numberOfSamples <- length(startInd)

  # every row in z one-hot encodes one character in sequence, oov is zero-vector
  num_classes <- voc_len + 2
  z  <- keras::to_categorical(sequence, num_classes = num_classes)[ , -c(1, num_classes)]

  if (use_quality) {
    ones_pos <- apply(z, 1, which.max)
    is_zero_row <- apply(z == 0, 1, all)
    z <- purrr::map(1:length(quality_vector), ~create_quality_vector(pos = ones_pos[.x], prob = quality_vector[.x],
                                                                     voc_length = length(vocabulary))) %>% unlist() %>%
      matrix(ncol = length(vocabulary), byrow = TRUE)
    z[is_zero_row, ] <- 0
  }

  if (ambiguous_nuc == "equal") {
    amb_nuc_pos <- which(sequence == (voc_len + 1))
    z[amb_nuc_pos, ] <- matrix(rep(1/voc_len, ncol(z) * length(amb_nuc_pos)), ncol = ncol(z))
  }
  if (ambiguous_nuc == "empirical") {
    if (!is.null(n_gram)) stop("Can only use equal, zero or discard option for ambiguous_nuc when using n_gram encoding")
    amb_nuc_pos <- which(sequence == (voc_len + 1))
    z[amb_nuc_pos, ] <- matrix(rep(nuc_dist, length(amb_nuc_pos)), nrow = length(amb_nuc_pos), byrow = TRUE)
  }

  if (use_coverage) {
    z <- z * (cov_vector/max_cov)
  }

  if (target_len == 1) {
    if (!target_middle) {
      if (!wavenet_format) {
        # target right
        x <- array(0, dim = c(numberOfSamples, maxlen, voc_len))
        for (i in 1:numberOfSamples) {
          start <- startInd[i]
          x[i, , ] <- z[start : (start + maxlen - 1), ]
        }
        #if (!is.null(n_gram)) x <- x[ , 1:(dim(x)[2] - n_gram + 1), ]

        y <- z[startInd + maxlen, ]
      } else if (!cnn_format) {
        # wavenet
        if (!is.null(n_gram)) stop("Wavenet format not implemented for n_gram.")
        x <- array(0, dim = c(numberOfSamples, maxlen, voc_len))
        y <- array(0, dim = c(numberOfSamples, maxlen, voc_len))
        for (i in 1:numberOfSamples) {
          start <- startInd[i]
          x[i, , ] <- z[start : (start + maxlen - 1), ]
          y[i, , ] <- z[(start + 1) : (start + maxlen), ]
        }
      } else {
        # target middle cnn
        x <- array(0, dim = c(numberOfSamples, maxlen + 1, voc_len))
        for (i in 1:numberOfSamples) {
          start <- startInd[i]
          x[i, , ] <- z[start : (start + maxlen), ]
        }
        missing_val <- ceiling(maxlen/2)
        y <- z[startInd + missing_val, ]
        #if (is.null(n_gram)) {
        x <- x[ , -(missing_val + 1), ]
        #} else {
        #  x <- x[ , -((missing_val - n_gram):(missing_val + target_len + n_gram - 1)), ]
        #}
      }
    } else {
      if (!wavenet_format) {
        len_input_1 <- ceiling(maxlen/2)
        len_input_2 <- floor(maxlen/2)
        input_tensor_1 <- array(0, dim = c(numberOfSamples, len_input_1, voc_len))
        input_tensor_2 <- array(0, dim = c(numberOfSamples, len_input_2, voc_len))
        for (i in 1:numberOfSamples) {
          start <- startInd[i]
          input_tensor_1[i, , ] <- z[start : (start + len_input_1 - 1), ]
          input_tensor_2[i, , ] <- z[(start + maxlen) : (start + len_input_1 + 1), ]
        }
        if (!is.null(n_gram)) {
          input_tensor_1 <- input_tensor_1[ , 1:(dim(input_tensor_1) - n_gram + 1), ]
          input_tensor_2 <- input_tensor_2[ , 1:(dim(input_tensor_2) - n_gram + 1), ]
        }
        x <- list(input_tensor_1, input_tensor_2)
        y <- z[startInd + len_input_1, ]
      } else {
        stop("Target middle not implemented for wavenet format.")
      }
    }
    # target_len > 1
  } else {
    if (!target_middle) {
      if (!cnn_format) {
        # target right
        x <- array(0, dim = c(numberOfSamples, maxlen - target_len + 1, voc_len))
        for (i in 1:numberOfSamples) {
          start <- startInd[i]
          x[i, , ] <- z[start : (start + maxlen - target_len), ]
        }
        y <- list()
        for (i in 1:target_len) {
          y[[i]] <- z[startInd + maxlen - target_len + i, ]
        }
        #if (!is.null(n_gram)) x <- x[ , 1:(dim(x)[2] - n_gram + 1), ]

      } else {
        # cnn format
        x <- array(0, dim = c(numberOfSamples, maxlen + 1, voc_len))
        for (i in 1:numberOfSamples) {
          start <- startInd[i]
          x[i, , ] <- z[start : (start + maxlen), ]
        }
        missing_val <- ceiling((maxlen - target_len)/2)
        y <- list()
        for (i in 1:target_len) {
          y[[i]] <- z[startInd + missing_val + i - 1, ]
        }
        #if (is.null(n_gram)) {
        x <- x[ , -((missing_val + 1):(missing_val + target_len)), ]
        #} else {
        #  x <- x[ , -((missing_val - n_gram):(missing_val + target_len + n_gram - 1)), ]
        #}
      }
    } else {
      if (!wavenet_format) {
        # target middle lstm
        len_input_1 <- ceiling((maxlen - target_len + 1)/2)
        len_input_2 <- maxlen + 1 - len_input_1 - target_len
        input_tensor_1 <- array(0, dim = c(numberOfSamples, len_input_1, voc_len))
        input_tensor_2 <- array(0, dim = c(numberOfSamples, len_input_2, voc_len))
        for (i in 1:numberOfSamples) {
          start <- startInd[i]
          input_tensor_1[i, , ] <- z[start : (start + len_input_1 - 1), ]
          input_tensor_2[i, , ] <- z[(start + maxlen) : (start + maxlen - len_input_2 + 1), ]
        }
        if (!is.null(n_gram)) {
          input_tensor_1 <- input_tensor_1[ , 1:(dim(input_tensor_1) - n_gram + 1), ]
          input_tensor_2 <- input_tensor_2[ , 1:(dim(input_tensor_2) - n_gram + 1), ]
        }
        x <- list(input_tensor_1, input_tensor_2)
        y <- list()
        for (i in 1:target_len) {
          y[[i]] <- z[startInd + len_input_1 - 1 + i, ]
        }
      } else {
        stop("Multi target not implemented for wavenet format.")
      }
    }
  }
  if (is.matrix(x)) {
    x <- array(x, dim = c(1, dim(x)))
  }

  if (!is.null(n_gram)) {
    y <- do.call(rbind, y)
    y_list <- list()
    for (i in 1:numberOfSamples) {
      index <- (i-1)  + (1 + (0:(target_len-1))*numberOfSamples)
      input_matrix <- y[index, ]
      if (length(index) == 1) input_matrix <- matrix(input_matrix, nrow = 1)
      n_gram_matrix <- n_gram_of_matrix(input_matrix = input_matrix, n = n_gram)
      y_list[[i]] <- n_gram_matrix # tensorflow::tf$expand_dims(n_gram_matrix, axis = 0L)
    }
    y_tensor <- keras::k_stack(y_list, axis = 1L) %>% keras::k_eval()
    y <- vector("list", dim(y_tensor)[2])

    for (i in 1:dim(y_tensor)[2]) {
      y_subset <- y_tensor[ , i, ]
      if (numberOfSamples == 1) y_subset <- matrix(y_subset, nrow = 1)
      y[[i]] <- y_subset
    }

    if (is.list(y) & length(y) == 1) {
      y <- y[[1]]
    }

    if (n_gram_stride > 1 & is.list(y)) {
      stride_index <- 0:(length(y)-1) %% n_gram_stride == 0
      y <- y[stride_index]
    }
  }

  return(list(x, y))
}

#' One-hot-encodes integer
#'
#' \code{sequenceToArrayLabel} Helper function for \code{\link{{fastaLabelGenerator}}, returns one hot encoding for sequence and returns samples from
#' specified positions
#'
#' @param sequence Sequence of integers.
#' @param maxlen Length of predictor sequence.
#' @param vocabulary Set of characters to encode.
#' @param startInd Start positions of samples in \code{sequence}.
#' @param ambiguous_nuc How to handle nucleotides outside vocabulary, either "zero", "discard" or "equal". If "zero", input gets encoded as zero vector;
#' if "equal" input is 1/length(vocabulary) x length(vocabulary). If "discard" samples containing nucleotides outside vocabulary get discarded.
#' @param nuc_dist Nucleotide distribution.
#' @param use_quality Use quality scores.
#' @param quality_vector Vector of quality probabilities.
#' @export
sequenceToArrayLabel <- function(sequence, maxlen, vocabulary, startInd, ambiguous_nuc = "zero", nuc_dist = NULL,
                                 use_quality = FALSE, quality_vector = NULL, use_coverage = FALSE, max_cov = NULL,
                                 cov_vector = NULL, n_gram = NULL) {

  voc_len <- length(vocabulary)
  if (!is.null(n_gram)) {
    sequence <- int_to_n_gram(int_seq = sequence, n = n_gram, voc_size = length(vocabulary))
    maxlen <- maxlen - n_gram + 1
    voc_len <- length(vocabulary)^n_gram
  }
  startInd <- startInd - startInd[1] + 1
  numberOfSamples <- length(startInd)

  # every row in z one-hot encodes one character in sequence, oov is zero-vector
  z  <- keras::to_categorical(sequence, num_classes = voc_len + 2)[ , -c(1, voc_len + 2)]
  z <- matrix(z, ncol = voc_len)

  if (use_quality) {
    ones_pos <- apply(z, 1, which.max)
    is_zero_row <- apply(z == 0, 1, all)
    z <- purrr::map(1:length(quality_vector), ~create_quality_vector(pos = ones_pos[.x], prob = quality_vector[.x],
                                                                     voc_length = voc_len)) %>% unlist() %>% matrix(ncol = voc_len, byrow = TRUE)
    z[is_zero_row, ] <- 0
  }

  if (ambiguous_nuc == "equal") {
    amb_nuc_pos <- which(sequence == (voc_len + 1))
    z[amb_nuc_pos, ] <- matrix(rep(1/voc_len, ncol(z) * length(amb_nuc_pos)), ncol = ncol(z))
  }
  if (ambiguous_nuc == "empirical") {
    amb_nuc_pos <- which(sequence == (voc_len + 1))
    z[amb_nuc_pos, ] <- matrix(rep(nuc_dist, length(amb_nuc_pos)), nrow = length(amb_nuc_pos), byrow = TRUE)
  }

  if (use_coverage) {
    z <- z * (cov_vector/max_cov)
  }

  x <- array(0, dim = c(numberOfSamples, maxlen, voc_len))
  for (i in 1:numberOfSamples) {
    start <- startInd[i]
    x[i, , ] <- z[start : (start + maxlen - 1), ]
  }
  return(x)
}

#' Computes start position of samples
#'
#' Helper function for \code{\link{{fastaLabelGenerator}} and \code{\link{{fastaFileGenerator}}. Computes positions in sequence where samples can be extracted
#'
#' @param seq_vector Vector of character sequences.
#' @param length_vector Length of sequences in \code{seq_vector}.
#' @param maxlen Length of one predictor sequence.
#' @param step Distance between samples from one entry in \code{seq_vector}.
#' @param train_mode Either "lm" for language model or "label" for label classification. Language models need one character more
#' (the target) for one sample.
#' @param ignore_amb_nuc Discard all samples that contain characters outside vocabulary.
#' @export
getStartInd <- function(seq_vector, length_vector, maxlen,
                        step, train_mode = "label", discard_amb_nuc = FALSE, vocabulary = c("A", "C", "G", "T")) {
  stopifnot(train_mode == "lm" | train_mode == "label")
  if (!discard_amb_nuc) {
    if (length(length_vector) > 1) {
      startNewEntry <- cumsum(c(1, length_vector[-length(length_vector)]))
      if (train_mode == "label") {
        indexVector <- purrr::map(1:(length(length_vector) - 1), ~seq(startNewEntry[.x], startNewEntry[.x + 1] - maxlen, by = step))
      } else {
        indexVector <- purrr::map(1:(length(length_vector) - 1), ~seq(startNewEntry[.x], startNewEntry[.x + 1] - maxlen - 1, by = step))
      }
      indexVector <- unlist(indexVector)
      last_seq <- length(seq_vector)
      if (!(startNewEntry[last_seq] > (sum(length_vector) - maxlen + 1))) {
        if (train_mode == "label") {
          indexVector <- c(indexVector, seq(startNewEntry[last_seq], sum(length_vector) - maxlen + 1, by = step))
        } else {
          indexVector <- c(indexVector, seq(startNewEntry[last_seq], sum(length_vector) - maxlen, by = step))
        }
      }
      return(indexVector)
    } else {
      if (train_mode == "label") {
        indexVector <- seq(1, length_vector - maxlen + 1, by = step)
      } else {
        indexVector <- seq(1, length_vector - maxlen, by = step)
      }
    }
  } else {
    indexVector <- startIndicesIgnoreAmbNuc(seq_vector = seq_vector, length_vector = length_vector,
                                            maxlen = maxlen, step = step, vocabulary = c(vocabulary, "0"), train_mode = train_mode)
  }
  return(indexVector)
}


#' Helper function for getStartInd, extracts the start positions of all potential samples (considering step size and vocabulary)
#'
#' @param seq Sequences.
#' @param length_vector Length of sequences in \code{seq_vector}.
#' @param maxlen Length of one sample.
#' @param step How often to take a sample.
#' @param vocabulary Vector of allowed characters in samples.
#' @param train_mode "lm" or "label".
#' @export
startIndicesIgnoreAmbNucSingleSeq <- function(seq, maxlen, step, vocabulary, train_mode = "lm") {

  vocabulary <- stringr::str_to_lower(vocabulary)
  vocabulary <- c(vocabulary, "0")
  seq <- stringr::str_to_lower(seq)
  len_seq <- nchar(seq)
  if (train_mode != "label") maxlen <- maxlen + 1
  stopifnot(len_seq >= maxlen)
  # regular expressions for allowed characters
  voc_pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")

  pos_of_amb_nucleotides <- stringr::str_locate_all(seq, pattern = voc_pattern)[[1]][ , 1]
  non_start_index <-  pos_of_amb_nucleotides - maxlen + 1

  # define range of unallowed start indices
  non_start_index <- purrr::map(non_start_index, ~(.x:(.x + maxlen - 1))) %>%
    unlist() %>% union((len_seq - maxlen + 2):len_seq) %>% unique()
  # drop non-positive values
  if (length(non_start_index[non_start_index < 1])) {
    non_start_index <- unique(c(1, non_start_index[non_start_index >= 1]))
  }

  non_start_index <- non_start_index %>% sort()
  allowed_start <- setdiff(1:len_seq, non_start_index)
  len_start_vector <- length(allowed_start)


  if (len_start_vector < 1) {
    # message("Can not extract a single sampling point with current settings.")
    return(NULL)
  }

  # only keep indices with sufficient distance, as defined by step
  start_indices <- vector("integer")
  index <- allowed_start[1]
  start_indices[1] <- index
  count <- 1
  if (length(allowed_start) > 1) {
    for (j in 1:(length(allowed_start) - 1)) {
      if (allowed_start[j + 1] - index >= step) {
        count <- count + 1
        start_indices[count] <- allowed_start[j + 1]
        index <- allowed_start[j + 1]
      }
    }
  }

  start_indices
}


#' Helper function for getStartInd, extracts the start positions of all potential samples (considering step size and vocabulary)
#'
#' @param seq_vector Vector of character sequences.
#' @param length_vector Length of sequences in \code{seq_vector}.
#' @param maxlen Length of one sample.
#' @param step How often to take a sample.
#' @param vocabulary Vector of allowed characters in samples.
#' @param train_mode "lm" or "label".
#' @export
startIndicesIgnoreAmbNuc <- function(seq_vector, length_vector, maxlen, step, vocabulary, train_mode = "lm") {
  startInd <- purrr::map(1:length(seq_vector), ~startIndicesIgnoreAmbNucSingleSeq(seq = seq_vector[.x],
                                                                                  maxlen = maxlen,
                                                                                  step = step,
                                                                                  vocabulary = vocabulary,
                                                                                  train_mode = train_mode))

  cum_sum_length <- cumsum(length_vector)
  if (length(startInd) > 1) {
    for (i in 2:length(startInd)) {
      startInd[[i]] <- startInd[[i]] + cum_sum_length[i - 1]
    }
  }
  startInd <- unlist(startInd)
  startInd
}

#' convert fastq quality score to probability
#'
#' @export
quality_to_probability <- function(quality_vector) {
  Q <- utf8ToInt(quality_vector) - 33
  1 - 10^(-Q/10)
}

#' @export
create_quality_vector <- function(pos, prob, voc_length = 4) {
  vec <- rep(0, voc_length)
  vec[pos] <- prob
  vec[-pos] <- (1 - prob)/(voc_length - 1)
  vec
}

#' @export
remove_amb_nuc_entries <- function(fasta.file, skip_amb_nuc, pattern) {
  chars_per_row <- nchar(fasta.file$Sequence)
  amb_per_row <- stringr::str_count(stringr::str_to_lower(fasta.file$Sequence), pattern)
  threshold_index <- (amb_per_row/chars_per_row) > skip_amb_nuc
  fasta.file <- fasta.file[!threshold_index, ]
  fasta.file
}

#' Estimate frequency of different classes
#'
#' @inheritParams fastaFileGenerator
#' @inheritParams fastaLabelGenerator
#' @param file_proportion Proportion of files to randomly sample for estimating class distributions.
#' @param train_type Either "label_folder", "label_header" or "label_csv".
#' @export
get_class_weight <- function(path,
                             labelVocabulary = NULL,
                             format = "fasta",
                             # estimate class distribution from subset
                             file_proportion = 1,
                             train_type = "label_folder",
                             csv_path = NULL) {

  classes <- count_nuc(path = path,
                       labelVocabulary = labelVocabulary,
                       format = format,
                       file_proportion = file_proportion,
                       train_type = train_type,
                       csv_path = csv_path)

  zero_entry <- classes == 0
  if (sum(zero_entry) > 0) {
    warning_message <- paste("The following classes have no samples:", paste(labelVocabulary[zero_entry]),
                             "\n Try bigger file_proportion size or check labelVocabulary.")
    warning(warning_message)
  }

  if (!is.list(classes)) {
    weight_list <- list()
    weight_list[["0"]] <- 1
    for (i in 2:(length(classes))) {
      weight_list[[as.character(i-1)]] <- classes[1]/classes[i]
    }
    classes <- weight_list
  } else {
    weight_collection <- list()
    for (j in 1:length(classes)) {
      weight_list <- list()
      weight_list[["0"]] <- 1
      for (i in 2:(length(classes[[j]]))) {
        weight_list[[as.character(i-1)]] <- classes[[j]][1]/classes[[j]][i]
      }
      weigth_collection[[j]] <- weight_list
    }
    classes <- weight_collection
  }

  classes
}

#' Count nucleotides per class
#'
#' @inheritParams get_class_weight
#' @export
count_nuc <- function(path,
                      labelVocabulary = NULL,
                      format = "fasta",
                      # estimate class distribution from subset
                      file_proportion = 1,
                      train_type = "label_folder",
                      csv_path = NULL) {

  classes <- rep(0, length(labelVocabulary))
  names(classes) <- labelVocabulary

  # label by folder
  if (train_type == "label_folder") {
    for (j in 1:length(path)) {
      files <- list.files(path[[j]], full.names = TRUE)
      if (file_proportion < 1) {
        files <- sample(files, floor(file_proportion * length(files)))
      }
      for (i in files) {
        if (format == "fasta") {
          fasta.file <- microseq::readFasta(i)
        }
        if (format == "fastq") {
          fasta.file <- microseq::readFastq(i)
        }
        freq <- sum(nchar(fasta.file$Sequence))
        classes[i] <- classes[i] + freq
      }
    }
  }

  # label header
  if (train_type == "label_header") {
    files <- list.files(unlist(path), full.names = TRUE)
    if (file_proportion < 1) {
      files <- sample(files, floor(file_proportion * length(files)))
    }
    for (i in files) {
      if (format == "fasta") {
        fasta.file <- microseq::readFasta(i)
      }
      if (format == "fastq") {
        fasta.file <- microseq::readFastq(i)
      }
      df <- data.frame(Header = fasta.file$Header, freq = nchar(fasta.file$Sequence))
      df <- aggregate(df$freq, by = list(Category = df$Header), FUN = sum)
      freq <- df$x
      names(freq) <- df$Category
      for (k in names(freq)) {
        classes[k] <- classes[k] + freq[k]
      }
    }
  }

  # label csv
  if (train_type == "label_csv") {

    label_csv <- read.csv2(csv_path, header = TRUE, stringsAsFactors = FALSE)
    if (dim(label_csv)[2] == 1) {
      label_csv <- read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)
    }
    if (!("file" %in% names(label_csv))) {
      stop('csv file needs one column named "file"')
    }

    row_sums <- label_csv %>% dplyr::select(-file) %>% rowSums()
    if (!(all(row_sums == 1))) {
      stop("Can only estimate class weights if labels are mutually exclusive.")
    }

    if (is.null(labelVocabulary) || missing(labelVocabulary)) {
      labelVocabulary <-  names(label_csv)[!names(label_csv) == "file"]
    } else {
      label_csv <- label_csv %>% dplyr::select(c(dplyr::all_of(labelVocabulary), "file"))
    }

    classes <- rep(0, length(labelVocabulary))
    names(classes) <- labelVocabulary

    path <- unlist(path)
    single_file_index <- stringr::str_detect(path, "fasta$|fastq$")
    files <- c(list.files(path[!single_file_index], full.names = TRUE), path[single_file_index])
    if (file_proportion < 1) {
      files <- sample(files, floor(file_proportion * length(files)))
    }
    for (i in files) {
      if (format == "fasta") {
        fasta.file <- microseq::readFasta(i)
      }
      if (format == "fastq") {
        fasta.file <- microseq::readFastq(i)
      }
      count_nuc <- sum(nchar(fasta.file$Sequence))
      df <- label_csv %>% dplyr::filter(file == basename(i))
      if (nrow(df) == 0) next
      index <- df[1, ] == 1
      current_label <- names(df)[index]
      classes[current_label] <- classes[current_label] + count_nuc
    }
  }
  return(classes)
}

#' Take random subset
#'
#' @export
random_subset <- function(proportion_per_file, fasta.file, use_beta_dist = FALSE, use_quality_score) {
  # take random subset
  if (!is.null(proportion_per_file)) {
    fasta_width <- nchar(fasta.file$Sequence)
    sample_range <- floor(fasta_width - (proportion_per_file * fasta_width))
    if (use_beta_dist) {
      interval_list <- purrr::map(1:length(sample_range), ~seq(from = 0, to = 1, len = sample_range[.x])[-1])
      beta <- rbeta(length(sample_range), shape1 = 0.95, shape2 = 0.95, ncp = 0)
      start <- purrr::map(1:length(sample_range), ~min(which(beta[.x] < interval_list[[i]]))) %>% unlist()
    } else {
      start <- mapply(sample_range, FUN = sample, size = 1)
    }
    perc_length <- floor(fasta_width * proportion_per_file)
    stop <- start + perc_length
    seq_vector <- mapply(fasta.file$Sequence, FUN = substr, start = start, stop = stop)
    if (use_quality_score) {
      quality_scores <- mapply(fasta.file$Quality, FUN = substr, start = start, stop = stop)
    }
  } else {
    seq_vector <- fasta.file$Sequence
    if (use_quality_score) {
      quality_scores <- fasta.file$Quality
    }
    if (use_quality_score) {
      return(list(seq_vector = seq_vector, quality_scores = quality_scores))
    } else {
      return(seq_vector = seq_vector)
    }
  }
}

#' @export
read_fasta_fastq <- function(format, skip_amb_nuc, file_index, pattern, shuffleFastaEntries,
                             reverseComplements, fasta.files, use_coverage = FALSE, proportion_entries = NULL,
                             labelVocabulary = NULL, filter_header = FALSE, target_from_csv = NULL) {
  if (format == "fasta") {
    if (is.null(skip_amb_nuc)) {
      fasta.file <- microseq::readFasta(fasta.files[file_index])
    } else {
      fasta.file <- remove_amb_nuc_entries(microseq::readFasta(fasta.files[file_index]), skip_amb_nuc = skip_amb_nuc,
                                           pattern = pattern)
    }

    if (filter_header & is.null(target_from_csv)) {
      label_vector <- trimws(stringr::str_to_lower(fasta.file$Header))
      label_filter <- label_vector %in% labelVocabulary
      fasta.file <- fasta.file[label_filter, ]
    }

    if (!is.null(proportion_entries) && proportion_entries < 1) {
      index <- sample(nrow(fasta.file), max(1, floor(nrow(fasta.file) * proportion_entries)))
      fasta.file <- fasta.file[index, ]
    }

    if (shuffleFastaEntries) {
      fasta.file <- fasta.file[sample(nrow(fasta.file)), ]
    }

    if (reverseComplements & sample(c(TRUE, FALSE), 1)) {
      fasta.file$Sequence <- microseq::reverseComplement(fasta.file$Sequence)
    }

  }
  if (format == "fastq") {
    if (is.null(skip_amb_nuc)) {
      fasta.file <- microseq::readFastq(fasta.files[file_index])
    } else {
      fasta.file <- remove_amb_nuc_entries(microseq::readFastq(fasta.files[file_index]), skip_amb_nuc = skip_amb_nuc,
                                           pattern = pattern)
    }

    if (filter_header & is.null(target_from_csv)) {
      label_vector <- trimws(stringr::str_to_lower(fasta.file$Header))
      label_filter <- label_vector %in% labelVocabulary
      fasta.file <- fasta.file[label_filter, ]
    }

    if (!is.null(proportion_entries) && proportion_entries < 1) {
      index <- sample(nrow(fasta.file), max(1, floor(nrow(fasta.file) * proportion_entries)))
      fasta.file <- fasta.file[index, ]
    }

    if (shuffleFastaEntries) {
      fasta.file <- fasta.file[sample(nrow(fasta.file)), ]
    }

    if (reverseComplements & sample(c(TRUE, FALSE), 1)) {
      fasta.file$Sequence <- microseq::reverseComplement(fasta.file$Sequence)
    }
  }
  return(fasta.file)
}

#' @export
input_from_csv <- function(added_label_path) {
  .datatable.aware = TRUE
  label_csv <- read.csv2(added_label_path, header = TRUE, stringsAsFactors = FALSE)
  if (dim(label_csv)[2] == 1) {
    label_csv <- read.csv(added_label_path, header = TRUE, stringsAsFactors = FALSE)
  }
  label_csv <- data.table::as.data.table(label_csv)
  label_csv$file <- stringr::str_to_lower(as.character(label_csv$file))
  data.table::setkey(label_csv, file)
  #   data.table::setkey(label_csv, file)
  # if ("file" %in% names(label_csv) & "header" %in% names(label_csv)) {
  #   stop('names in added_label_path should contain "header" or "file" not both')
  # } else if ("header" %in% names(label_csv)) {
  #   added_label_by_header <- TRUE
  #   label_csv$file <- stringr::str_to_lower(as.character(label_csv$header))
  #   data.table::setkey(label_csv, header)
  # } else if ("file" %in% names(label_csv)) {
  #   added_label_by_header <- FALSE
  #   label_csv$file <- stringr::str_to_lower(as.character(label_csv$file))
  #   data.table::setkey(label_csv, file)
  # } else {
  #   stop('file in added_label_path must contain one column named "header" or "file"')
  # }
  #
  # if (!added_label_by_header) {
  #   fasta.file$Header <- rep(basename(fasta.files[file_index]), nrow(fasta.file))
  # }
  added_label_by_header <- FALSE

  if (!("file" %in% names(label_csv))) {
    stop('names in added_label_path should contain one column named "file" ')
  }
  col_name <- ifelse(added_label_by_header, "header", "file")
  # header_vector <- fasta.file$Header
  # return(list(label_csv = label_csv, col_name = col_name, header_vector = header_vector))
  return(list(label_csv = label_csv, col_name = col_name))
}

#' @import data.table
#' @export
csv_to_tensor <- function(label_csv, added_label_vector, added_label_by_header, batch.size,
                          start_index_list) {
  .datatable.aware = TRUE
  label_tensor <- matrix(0, ncol = ncol(label_csv) - 1, nrow = batch.size, byrow = TRUE)

  if (added_label_by_header) {
    header_unique <- unique(added_label_vector)
    for (i in header_unique) {
      label_from_csv <- label_csv[ .(i), -"header"]
      index_label_vector <- added_label_vector == i
      if (nrow(label_from_csv) > 0) {
        label_tensor[index_label_vector, ] <- matrix(as.matrix(label_from_csv[1, ]),
                                                     nrow = sum(index_label_vector), ncol = ncol(label_tensor), byrow = TRUE)
      }
    }
  } else {
    row_index <- 1
    for (i in 1:length(added_label_vector)) {
      row_filter <- added_label_vector[i]
      label_from_csv <- label_csv[data.table(row_filter), -"file"]
      samples_per_file <- length(start_index_list[[i]])
      assign_rows <-  row_index:(row_index + samples_per_file - 1)

      if (nrow(na.omit(label_from_csv)) > 0) {
        label_tensor[assign_rows, ] <- matrix(as.matrix(label_from_csv[1, ]),
                                              nrow = samples_per_file, ncol = ncol(label_tensor), byrow = TRUE)
      }
      row_index <- row_index + samples_per_file
    }
  }
  return(label_tensor)
}

#' Devide tensor to list of subsets
#'
#' @export
slice_tensor <- function(tensor, target_split) {

  num_row <- nrow(tensor)
  l <- vector("list", length = length(target_split))
  for (i in 1:length(target_split)) {
    if (length(target_split[[i]]) == 1 | num_row == 1) {
      l[[i]] <- matrix(tensor[ , target_split[[i]]], ncol = length(target_split[[i]]))
    } else {
      l[[i]] <- tensor[ , target_split[[i]]]
    }
  }
  return(l)
}

check_header_names <- function(target_split, labelVocabulary) {
  target_split <- unlist(target_split)
  if (!all(target_split %in% labelVocabulary)) {
    stop_text <- paste("Your csv file has no columns named",
                       paste(target_split[!(target_split %in% labelVocabulary)], collapse = " "))
    stop(stop_text)
  }
  if (!all(labelVocabulary %in% target_split)) {
    warning_text <- paste("target_split does not cover the following columns:",
                          paste(labelVocabulary[!(labelVocabulary %in% target_split)], collapse = " "))
    warning(warning_text)
  }
}

count_files <- function(path, format = "fasta", train_type) {
  num_files <- rep(0, length(path))
  for (i in 1:length(path)) {
    for (k in 1:length(path[[i]])) {
      if (endsWith(path[[i]][[k]], paste0(".", format))) {
        current_files <- 1
      } else {
        current_files <- length(list.files(path[[i]][[k]], pattern = paste0(".", format)))
      }
      num_files[i] <- num_files[i] + current_files

      if (current_files == 0) {
        stop(paste0(path[[i]][[k]], " is empty or no files with .", format, " ending in this directory"))
      }
    }
  }
  # return number of files per class for "label_folder"
  if (train_type == "label_folder") {
    return(num_files)
  } else {
    return(sum(num_files))
  }
}

list_fasta_files <- function(corpus.dir, format, file_filter) {
  #if (is.list(corpus.dir)) {
  fasta.files <- list()
  for (i in 1:length(corpus.dir)) {

    if (endsWith(corpus.dir[[i]], paste0(".", format))) {
      fasta.files[[i]] <- corpus.dir[[i]]

    } else {

      fasta.files[[i]] <- list.files(
        path = xfun::normalize_path(corpus.dir[[i]]),
        pattern = paste0("\\.", format, "$"),
        full.names = TRUE)
    }
  }
  fasta.files <- unlist(fasta.files)
  num_files <- length(fasta.files)
  # } else {
  #   fasta.files <- list()
  #   # single file
  #   if (endsWith(corpus.dir, paste0(".", format))) {
  #     num_files <- 1
  #     fasta.files <- corpus.dir
  #   } else {
  #
  #     fasta.files <- list.files(
  #       path = xfun::normalize_path(corpus.dir),
  #       pattern = paste0("\\.", format, "$"),
  #       full.names = TRUE)
  #     num_files <- length(fasta.files)
  #   }
  # }

  if (!is.null(file_filter)) {
    fasta.files <- fasta.files[basename(fasta.files) %in% file_filter]
    if (length(fasta.files) < 1) {
      stop_text <- paste0("None of the files from ", unlist(corpus.dir),
                          " are present in train_val_split_csv table for either train or validation. \n")
      stop(stop_text)
    }
  }

  return(fasta.files)
}

get_coverage <- function(fasta.file) {
  header <- fasta.file$Header
  cov <- stringr::str_extract(header, "cov_\\d+") %>%
    stringr::str_extract("\\d+") %>% as.integer()
  cov[is.na(cov)] <- 1
  return(cov)
}

get_coverage_concat <- function(fasta.file, concat_seq) {
  header <- fasta.file$Header
  cov <- stringr::str_extract(header, "cov_\\d+") %>%
    stringr::str_extract("\\d+") %>% as.integer()
  cov[is.na(cov)] <- 1
  len_vec <- nchar(fasta.file$Sequence)
  cov <- purrr::map(1:nrow(fasta.file), ~rep(cov[.x], times = len_vec[.x]))
  cov <- lapply(cov, append, rep(1, nchar(concat_seq)))
  cov <- unlist(cov)
  cov <- cov[-((length(cov) - nchar(concat_seq)) : length(cov))]
  return(cov)
}


reshape_tensor <- function(x, y, new_batch_size,
                           samples_per_target,
                           batch.size, path, voc_len,
                           buffer_len = NULL,
                           maxlen, reshape_mode = "time_dist",
                           concat_maxlen = NULL) {

  if (reshape_mode == "time_dist") {

    x_new <- array(0, dim = c(new_batch_size, samples_per_target, maxlen, voc_len))
    y_new <- array(0, dim = c(new_batch_size, length(path)))
    for (i in 1:new_batch_size) {
      index <- (1:samples_per_target) + (i-1)*samples_per_target
      x_new[i, , , ] <- x[index, , ]
      y_new[i, ]  <- y[index[1], ]
    }

    return(list(x = x_new, y = y_new))
  }

  if (reshape_mode == "multi_input") {

    x_list <- vector("list", samples_per_target)
    for (i in 1:samples_per_target) {
      x_index <- base::seq(i, batch.size, samples_per_target)
      x_list[[i]] <- x[x_index, , ]
    }
    y <- y[base::seq(1, batch.size, samples_per_target), ]
    return(list(x = x_list, y = y))
  }

  if (reshape_mode == "concat") {

    x_new <- array(0, dim = c(new_batch_size, concat_maxlen, voc_len))
    y_new <- array(0, dim = c(new_batch_size, length(path)))
    use_buffer <- !is.null(buffer_len) && buffer_len > 0
    if (use_buffer) {
      buffer_tensor <- array(0, dim = c(buffer_len, voc_len))
      buffer_tensor[ , voc_len] <- 1
    }

    for (i in 1:new_batch_size) {
      index <- (1:samples_per_target) + (i-1)*samples_per_target
      if (!use_buffer) {
        x_temp <- x[index, , ]
        x_temp <- reticulate::array_reshape(x_temp, dim = c(1, dim(x_temp)[1] * dim(x_temp)[2], voc_len))
      } else {
        # create list of subsequences interspaced with buffer tensor
        x_list <- vector("list", (2*samples_per_target) - 1)
        x_list[seq(2, length(x_list), by = 2)] <- list(buffer_tensor)
        for (k in 1:length(index)) {
          x_list[[(2*k) - 1]] <- x[index[k], , ]
        }
        x_temp <- do.call(rbind, x_list)
        # # choose random region of size concat_maxlen
        # start_index <- sample(nrow(x_temp) - concat_maxlen + 1, 1)
        # end_index <- start_index + concat_maxlen - 1
        # x_temp <- x_temp[start_index : end_index,  ]
      }

      x_new[i, , ] <- x_temp
      y_new[i, ]  <- y[index[1], ]
    }
    return(list(x = x_new, y = y_new))
  }
}

#' Transform confusion matrix with total numbers to matrix with percentages.
#'
#' @export
cm_perc <- function(cm, round_dig = 2) {
  col_sums <- colSums(cm)
  for (i in 1:ncol(cm)) {
    if (col_sums[i] == 0) {
      cm[ , i] <- 0
    } else {
      cm[ , i] <- cm[ , i]/col_sums[i]
    }
  }
  cm <- round(cm, round_dig)
}

create_conf_mat_obj <- function(m, confMatLabels) {
  dimnames(m) <- list(Truth = confMatLabels, Prediction = confMatLabels)
  l <- list()
  m <- as.table(m)
  l[["table"]] <- m
  l[["dots"]] <- list()
  class(l) <- "conf_mat"
  return(l)
}

int_to_n_gram <- function(int_seq, n, voc_size = 4) {
  encoding_len <- length(int_seq) - n + 1
  n_gram_encoding <- vector("numeric", encoding_len)
  oov_token <- voc_size^n + 1
  for (i in 1:encoding_len) {
    int_seq_subset <- int_seq[i:(i + n - 1)]

    if (any(int_seq_subset > voc_size) | prod(int_seq_subset) == 0) {
      # zero encoding for amb nuc
      n_gram_encoding[i] <- oov_token
    } else {
      int_seq_subset <- int_seq_subset - 1
      n_gram_encoding[i] <- 1 + sum(voc_size^((n-1):0) * (int_seq_subset))
    }
  }
  n_gram_encoding
}

n_gram_of_matrix <- function(input_matrix, n = 3) {
  voc_len <- ncol(input_matrix)^n
  oov_index <- apply(input_matrix, 1, max) != 1
  max_index <- apply(input_matrix, 1, which.max)
  max_index[oov_index] <- voc_len + 1
  int_enc <- int_to_n_gram(int_seq = max_index, n = n, voc_size = ncol(input_matrix))
  if (length(int_enc) == 1) {
    n_gram_matrix <- matrix(keras::to_categorical(int_enc, num_classes = voc_len + 2), nrow = 1)[ , -c(1, voc_len + 2)]
  } else {
    n_gram_matrix <- keras::to_categorical(int_enc, num_classes = voc_len + 2)[ , -c(1, voc_len + 2)]
  }
  n_gram_matrix <- matrix(n_gram_matrix, ncol = voc_len)
  return(n_gram_matrix)
}

n_gram_of_3d_tensor <- function(tensor_3d, n) {
  new_dim <- dim(tensor_3d)
  new_dim[2] <- new_dim[2] - n + 1
  new_dim[3] <- new_dim[3]^n
  new_tensor <- array(0, dim = new_dim)
  for (i in 1:dim(tensor_3d)[1]) {
    new_tensor[i, , ] <- n_gram_of_matrix(tensor_3d[i, , ], n = n)
  }
  new_tensor
}
