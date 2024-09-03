#' Encodes integer sequence for language model
#'
#' Helper function for \code{\link{generator_fasta_lm}}. 
#' Encodes integer sequence to input/target list according to \code{output_format} argument. 
#'
#' @inheritParams generator_fasta_lm
#' @param sequence Sequence of integers.
#' @param start_ind Start positions of samples in \code{sequence}.
#' @param ambiguous_nuc How to handle nucleotides outside vocabulary, either `"zero"`, `"empirical"` or `"equal"`.
#' See \code{\link{train_model}}. Note that `"discard"` option is not available for this function.
#' @param nuc_dist Nucleotide distribution.
#' @param max_cov Biggest coverage value. Only applies if `use_coverage = TRUE`.
#' @param cov_vector Vector of coverage values associated to the input. 
#' @param adjust_start_ind Whether to shift values in \code{start_ind} to start at 1: for example (5,11,25) becomes (1,7,21).
#' @param quality_vector Vector of quality probabilities.
#' @param tokenizer A keras tokenizer.
#' @param char_sequence A character string.
#' @examplesIf reticulate::py_module_available("tensorflow")
#' # use integer sequence as input 
#' 
#' z <- seq_encoding_lm(sequence = c(1,0,5,1,3,4,3,1,4,1,2),
#' maxlen = 5,
#' vocabulary = c("a", "c", "g", "t"),
#' start_ind = c(1,3),
#' ambiguous_nuc = "equal",
#' target_len = 1,
#' output_format = "target_right")
#' 
#' x <- z[[1]]
#' y <- z[[2]]
#' 
#' x[1,,] # 1,0,5,1,3
#' y[1,] # 4
#' 
#' x[2,,] # 5,1,3,4,
#' y[2,] # 1
#' 
#' # use character string as input
#' z <- seq_encoding_lm(sequence = NULL,
#' maxlen = 5,
#' vocabulary = c("a", "c", "g", "t"),
#' start_ind = c(1,3),
#' ambiguous_nuc = "zero",
#' target_len = 1,
#' output_format = "target_right",
#' char_sequence = "ACTaaTNTNaZ")
#' 
#' 
#' x <- z[[1]]
#' y <- z[[2]]
#' 
#' x[1,,] # actaa
#' y[1,] # t
#' 
#' x[2,,] # taatn
#' y[2,] # t
#' 
#' @returns A list of 2 tensors.
#' @export
seq_encoding_lm <- function(sequence = NULL, maxlen, vocabulary, start_ind, ambiguous_nuc = "zero",
                            nuc_dist = NULL, quality_vector = NULL, return_int = FALSE,
                            target_len = 1, use_coverage = FALSE, max_cov = NULL, cov_vector = NULL,
                            n_gram = NULL, n_gram_stride = 1, output_format = "target_right",
                            char_sequence = NULL, adjust_start_ind = FALSE,
                            tokenizer = NULL) {
  
  use_quality <- ifelse(is.null(quality_vector), FALSE, TRUE)
  discard_amb_nt <- FALSE
  ## TODO: add discard_amb_nt
  if (!is.null(char_sequence)) {
    
    vocabulary <- stringr::str_to_lower(vocabulary)
    pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
    
    
    # token for ambiguous nucleotides
    for (i in letters) {
      if (!(i %in% stringr::str_to_lower(vocabulary))) {
        amb_nuc_token <- i
        break
      }
    }
    
    if (is.null(tokenizer)) {
      tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), c(vocabulary, amb_nuc_token))
    }
    
    sequence <- stringr::str_to_lower(char_sequence)
    sequence <- stringr::str_replace_all(string = sequence, pattern = pattern, amb_nuc_token)
    sequence <- keras::texts_to_sequences(tokenizer, sequence)[[1]] - 1
  }
  
  voc_len <- length(vocabulary)
  if (target_len == 1) {
    n_gram <- NULL
  }
  if (!is.null(n_gram)) {
    if (target_len < n_gram) stop("target_len needs to be at least as big as n_gram")
  }
  
  if (adjust_start_ind) start_ind <- start_ind - start_ind[1] + 1
  numberOfSamples <- length(start_ind)
  
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
    
    if (output_format == "target_right") {
      x <- array(0, dim = c(numberOfSamples, maxlen, voc_len))
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        x[i, , ] <- z[start : (start + maxlen - 1), ]
      }
      y <- z[start_ind + maxlen, ]
    }
    
    if (output_format == "wavenet") {
      if (!is.null(n_gram)) stop("Wavenet format not implemented for n_gram.")
      x <- array(0, dim = c(numberOfSamples, maxlen, voc_len))
      y <- array(0, dim = c(numberOfSamples, maxlen, voc_len))
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        x[i, , ] <- z[start : (start + maxlen - 1), ]
        y[i, , ] <- z[(start + 1) : (start + maxlen), ]
      }
    }
    
    if (output_format == "target_middle_cnn") {
      x <- array(0, dim = c(numberOfSamples, maxlen + 1, voc_len))
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        x[i, , ] <- z[start : (start + maxlen), ]
      }
      missing_val <- ceiling(maxlen/2)
      y <- z[start_ind + missing_val, ]
      x <- x[ , -(missing_val + 1), ]
    }
    
    if (output_format == "target_middle_lstm") {
      len_input_1 <- ceiling(maxlen/2)
      len_input_2 <- floor(maxlen/2)
      input_tensor_1 <- array(0, dim = c(numberOfSamples, len_input_1, voc_len))
      input_tensor_2 <- array(0, dim = c(numberOfSamples, len_input_2, voc_len))
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        input_tensor_1[i, , ] <- z[start : (start + len_input_1 - 1), ]
        input_tensor_2[i, , ] <- z[(start + maxlen) : (start + len_input_1 + 1), ]
      }
      if (!is.null(n_gram)) {
        input_tensor_1 <- input_tensor_1[ , 1:(dim(input_tensor_1) - n_gram + 1), ]
        input_tensor_2 <- input_tensor_2[ , 1:(dim(input_tensor_2) - n_gram + 1), ]
      }
      x <- list(input_tensor_1, input_tensor_2)
      y <- z[start_ind + len_input_1, ]
    }
    
  }
  
  if (target_len > 1) {
    
    if (output_format == "target_right") {
      x <- array(0, dim = c(numberOfSamples, maxlen - target_len + 1, voc_len))
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        x[i, , ] <- z[start : (start + maxlen - target_len), ]
      }
      y <- list()
      for (i in 1:target_len) {
        y[[i]] <- z[start_ind + maxlen - target_len + i, ]
      }
    }
    
    if (output_format == "target_middle_cnn") {
      x <- array(0, dim = c(numberOfSamples, maxlen + 1, voc_len))
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        x[i, , ] <- z[start : (start + maxlen), ]
      }
      missing_val <- ceiling((maxlen - target_len)/2)
      y <- list()
      for (i in 1:target_len) {
        y[[i]] <- z[start_ind + missing_val + i - 1, ]
      }
      x <- x[ , -((missing_val + 1):(missing_val + target_len)), ]
    }
    
    if (output_format == "target_middle_lstm") {
      len_input_1 <- ceiling((maxlen - target_len + 1)/2)
      len_input_2 <- maxlen + 1 - len_input_1 - target_len
      input_tensor_1 <- array(0, dim = c(numberOfSamples, len_input_1, voc_len))
      input_tensor_2 <- array(0, dim = c(numberOfSamples, len_input_2, voc_len))
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        input_tensor_1[i, , ] <- z[start : (start + len_input_1 - 1), ]
        input_tensor_2[i, , ] <- z[(start + maxlen) : (start + maxlen - len_input_2 + 1), ]
      }
      
      x <- list(input_tensor_1, input_tensor_2)
      y <- list()
      for (i in 1:target_len) {
        y[[i]] <- z[start_ind + len_input_1 - 1 + i, ]
      }
    }
    
    if (output_format == "wavenet") {
      stop("Multi target not implemented for wavenet format.")
    }
  }
  
  if (is.matrix(x)) {
    x <- array(x, dim = c(1, dim(x)))
  }
  
  if (!is.null(n_gram)) {
    if (is.list(y)) y <- do.call(rbind, y)
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

#' Encodes integer sequence for label classification.
#'
#' Returns encoding for integer or character sequence.
#'
#' @inheritParams seq_encoding_lm
#' @inheritParams generator_fasta_lm
#' @inheritParams train_model
#' @param return_int Whether to return integer encoding or one-hot encoding.
#' @examplesIf reticulate::py_module_available("tensorflow")
#' # use integer sequence as input
#' x <- seq_encoding_label(sequence = c(1,0,5,1,3,4,3,1,4,1,2),
#'                         maxlen = 5,
#'                         vocabulary = c("a", "c", "g", "t"),
#'                         start_ind = c(1,3),
#'                         ambiguous_nuc = "equal")
#' 
#' x[1,,] # 1,0,5,1,3
#' 
#' x[2,,] # 5,1,3,4,
#' 
#' # use character string as input
#' x <- seq_encoding_label(maxlen = 5,
#'                         vocabulary = c("a", "c", "g", "t"),
#'                         start_ind = c(1,3),
#'                         ambiguous_nuc = "equal",
#'                         char_sequence = "ACTaaTNTNaZ")
#' 
#' x[1,,] # actaa
#' 
#' x[2,,] # taatn
#' 
#' @returns A list of 2 tensors.
#' @export
seq_encoding_label <- function(sequence = NULL, maxlen, vocabulary, start_ind, ambiguous_nuc = "zero", nuc_dist = NULL,
                               quality_vector = NULL, use_coverage = FALSE, max_cov = NULL,
                               cov_vector = NULL, n_gram = NULL, n_gram_stride = 1, masked_lm = NULL,
                               char_sequence = NULL, tokenizer = NULL, adjust_start_ind = FALSE,
                               return_int = FALSE) {
  
  ## TODO: add discard_amb_nt, add conditions for return_int
  use_quality <- ifelse(is.null(quality_vector), FALSE, TRUE)
  discard_amb_nt <- FALSE
  maxlen_original <- maxlen
  if (return_int) ambiguous_nuc <- "zero"
  
  if (!is.null(char_sequence)) {
    
    vocabulary <- stringr::str_to_lower(vocabulary)
    pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
    
    # token for ambiguous nucleotides
    for (i in letters) {
      if (!(i %in% stringr::str_to_lower(vocabulary))) {
        amb_nuc_token <- i
        break
      }
    }
    
    if (is.null(tokenizer)) {
      tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), c(vocabulary, amb_nuc_token))
    }
    
    sequence <- stringr::str_to_lower(char_sequence)
    sequence <- stringr::str_replace_all(string = sequence, pattern = pattern, amb_nuc_token)
    sequence <- keras::texts_to_sequences(tokenizer, sequence)[[1]] - 1
  }
  
  if (adjust_start_ind) start_ind <- start_ind - start_ind[1] + 1
  numberOfSamples <- length(start_ind)
  
  if (is.null(n_gram_stride)) n_gram_stride <- 1
  voc_len <- length(vocabulary)
  if (!is.null(n_gram)) {
    sequence <- int_to_n_gram(int_seq = sequence, n = n_gram, voc_size = length(vocabulary))
    maxlen <- ceiling((maxlen - n_gram + 1)/n_gram_stride)
    voc_len <- length(vocabulary)^n_gram
  }
  
  if (!is.null(masked_lm)) {
    l <- mask_seq(int_seq = sequence,
                  mask_rate = masked_lm$mask_rate,
                  random_rate = masked_lm$random_rate,
                  identity_rate = masked_lm$identity_rate,
                  start_ind = start_ind,
                  block_len = masked_lm$block_len,
                  voc_len = voc_len)
    masked_seq <- l$masked_seq
    sample_weight_seq <- l$sample_weight_seq
  }
  
  if (!return_int) {
    if (!is.null(masked_lm)) {
      # every row in z one-hot encodes one character in sequence, oov is zero-vector
      z_masked <- keras::to_categorical(masked_seq, num_classes = voc_len + 2)[ , -c(1)]
      z_masked <- matrix(z_masked, ncol = voc_len + 1)
      z <- keras::to_categorical(sequence, num_classes = voc_len + 2)[ , -c(1)]
      z <- matrix(z, ncol = voc_len + 1)
    } else {
      # every row in z one-hot encodes one character in sequence, oov is zero-vector
      z  <- keras::to_categorical(sequence, num_classes = voc_len + 2)[ , -c(1, voc_len + 2)]
      z <- matrix(z, ncol = voc_len)
    }
  }
  
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
  
  remove_end_of_seq <- ifelse(is.null(n_gram), 1, n_gram) 
  
  if (!return_int) {
    if (is.null(masked_lm)) {
      
      x <- array(0, dim = c(numberOfSamples, maxlen, voc_len))
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        subset_index <- seq(start, (start + maxlen_original - remove_end_of_seq), by = n_gram_stride)
        x[i, , ] <- z[subset_index, ]
      }
      return(x)
      
    } else {
      
      x <- array(0, dim = c(numberOfSamples, maxlen, voc_len + 1))
      y <- array(0, dim = c(numberOfSamples, maxlen, voc_len + 1))
      sw <- array(0, dim = c(numberOfSamples, maxlen))
      
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        subset_index <- seq(start, (start + maxlen - remove_end_of_seq), by = n_gram_stride)
        x[i, , ] <- z_masked[subset_index, ]
        y[i, , ] <- z[subset_index, ]
        sw[i, ] <- sample_weight_seq[subset_index]
      }
      return(list(x=x, y=y, sample_weight=sw))
      
    }
  }
  
  if (return_int) {
    if (is.null(masked_lm)) {
      
      x <- array(0, dim = c(numberOfSamples, maxlen))
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        subset_index <- seq(start, (start + maxlen_original - remove_end_of_seq), by = n_gram_stride)
        x[i, ] <- sequence[subset_index]
      }
      return(x)
      
    } else {
      x <- array(0, dim = c(numberOfSamples, maxlen))
      y <- array(0, dim = c(numberOfSamples, maxlen))
      sw <- array(0, dim = c(numberOfSamples, maxlen))
      for (i in 1:numberOfSamples) {
        start <- start_ind[i]
        subset_index <- seq(start, (start + maxlen_original - remove_end_of_seq), by = n_gram_stride)
        x[i, ] <- masked_seq[subset_index]
        y[i, ] <- sequence[subset_index]
        sw[i, ] <- sample_weight_seq[subset_index]
      }
      return(list(x=x, y=y, sample_weight=sw))
      
    }
  }
  
}

#' Computes start position of samples
#'
#' Helper function for data generators. 
#' Computes start positions in sequence where samples can be extracted, given maxlen, step size and ambiguous nucleotide constraints.
#'
#' @inheritParams train_model
#' @param seq_vector Vector of character sequences.
#' @param length_vector Length of sequences in \code{seq_vector}.
#' @param maxlen Length of one predictor sequence.
#' @param step Distance between samples from one entry in \code{seq_vector}.
#' @param train_mode Either `"lm"` for language model or `"label"` for label classification. 
#' @param discard_amb_nuc Whether to discard all samples that contain characters outside vocabulary.
#' @examples
#' seq_vector <- c("AAACCCNNNGGGTTT")
#' get_start_ind(
#'   seq_vector = seq_vector,
#'   length_vector = nchar(seq_vector),
#'   maxlen = 4,
#'   step = 2,
#'   train_mode = "label",
#'   discard_amb_nuc = TRUE,
#'   vocabulary = c("A", "C", "G", "T"))
#'   
#' @returns A numeric vector.   
#' @export
get_start_ind <- function(seq_vector, length_vector, maxlen,
                          step, train_mode = "label", 
                          discard_amb_nuc = FALSE,
                          vocabulary = c("A", "C", "G", "T")) {
  
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
    indexVector <- start_ind_ignore_amb(seq_vector = seq_vector, length_vector = length_vector,
                                        maxlen = maxlen, step = step, vocabulary = c(vocabulary, "0"), train_mode = train_mode)
  }
  return(indexVector)
}


#' Helper function for get_start_ind, extracts the start positions of all potential samples (considering step size and vocabulary)
#'
#' @param seq Sequences.
#' @param maxlen Length of one sample.
#' @param step How often to take a sample.
#' @param vocabulary Vector of allowed characters in samples.
#' @param train_mode "lm" or "label".
#' @noRd
start_ind_ignore_amb_single_seq <- function(seq, maxlen, step, vocabulary, train_mode = "lm") {
  
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


#' Helper function for get_start_ind, extracts the start positions of all potential samples (considering step size and vocabulary)
#'
#' @param seq_vector Vector of character sequences.
#' @param maxlen Length of one sample.
#' @param step How often to take a sample.
#' @param vocabulary Vector of allowed characters in samples.
#' @param train_mode "lm" or "label".
#' @noRd
start_ind_ignore_amb <- function(seq_vector, length_vector, maxlen, step, vocabulary, train_mode = "lm") {
  start_ind <- purrr::map(1:length(seq_vector), ~start_ind_ignore_amb_single_seq(seq = seq_vector[.x],
                                                                                 maxlen = maxlen,
                                                                                 step = step,
                                                                                 vocabulary = vocabulary,
                                                                                 train_mode = train_mode))
  
  cum_sum_length <- cumsum(length_vector)
  if (length(start_ind) > 1) {
    for (i in 2:length(start_ind)) {
      start_ind[[i]] <- start_ind[[i]] + cum_sum_length[i - 1]
    }
  }
  start_ind <- unlist(start_ind)
  start_ind
}

quality_to_probability <- function(quality_vector) {
  Q <- utf8ToInt(quality_vector) - 33
  1 - 10^(-Q/10)
}

create_quality_vector <- function(pos, prob, voc_length = 4) {
  vec <- rep(0, voc_length)
  vec[pos] <- prob
  vec[-pos] <- (1 - prob)/(voc_length - 1)
  vec
}

remove_amb_nuc_entries <- function(fasta.file, skip_amb_nuc, pattern) {
  chars_per_row <- nchar(fasta.file$Sequence)
  amb_per_row <- stringr::str_count(stringr::str_to_lower(fasta.file$Sequence), pattern)
  threshold_index <- (amb_per_row/chars_per_row) > skip_amb_nuc
  fasta.file <- fasta.file[!threshold_index, ]
  fasta.file
}

#' Estimate frequency of different classes
#' 
#' Count number of nucleotides for each class and use as estimation for relation of class distribution.
#' Outputs list of class relations. Can be used as input for \code{class_weigth} in \code{\link{train_model}} function.   
#'
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_header_csv
#' @inheritParams train_model
#' @param file_proportion Proportion of files to randomly sample for estimating class distributions.
#' @param csv_path If `train_type = "label_csv"`, path to csv file containing labels.
#' @param named_list Whether to give class weight list names `"0", "1", ...` or not.
#' @examples 
#' 
#' # create dummy data
#' path_1 <- tempfile()
#' path_2 <- tempfile()
#' 
#' for (current_path in c(path_1, path_2)) {
#'   
#'   dir.create(current_path)
#'   # create twice as much data for first class
#'   num_files <- ifelse(current_path == path_1, 6, 3)
#'   create_dummy_data(file_path = current_path,
#'                     num_files = num_files,
#'                     seq_length = 10,
#'                     num_seq = 5,
#'                     vocabulary = c("a", "c", "g", "t"))
#' }
#' 
#' 
#' class_weight <- get_class_weight(
#'   path = c(path_1, path_2),
#'   vocabulary_label = c("A", "B"),
#'   format = "fasta",
#'   file_proportion = 1,
#'   train_type = "label_folder",
#'   csv_path = NULL)
#' 
#' class_weight
#' 
#' @returns A list of numeric values (class weights).
#' @export
get_class_weight <- function(path,
                             vocabulary_label = NULL,
                             format = "fasta",
                             file_proportion = 1, 
                             train_type = "label_folder",
                             named_list = FALSE,
                             csv_path = NULL) {
  
  classes <- count_nuc(path = path,
                       vocabulary_label = vocabulary_label,
                       format = format,
                       file_proportion = file_proportion,
                       train_type = train_type,
                       csv_path = csv_path)
  
  zero_entry <- classes == 0
  if (sum(zero_entry) > 0) {
    warning_message <- paste("The following classes have no samples:", paste(vocabulary_label[zero_entry]),
                             "\n Try bigger file_proportion size or check vocabulary_label.")
    warning(warning_message)
  }
  
  if (!is.list(classes)) {
    num_classes <- length(classes)
    total <- sum(classes)
    weight_list <- list()
    for (i in 1:(length(classes))) {
      weight_list[[as.character(i-1)]] <- total/(classes[i] * num_classes)
    }
    if (!named_list) names(classes) <- NULL # no list names in tf version > 2.8
    classes <- weight_list
  } else {
    weight_collection <- list()
    for (j in 1:length(classes)) {
      num_classes <- length(classes[[j]])
      total <- sum(classes[[j]])
      weight_list <- list()
      for (i in 1:(length(classes[[j]]))) {
        weight_list[[as.character(i-1)]] <- total/(classes[[j]][i] * num_classes)
      }
      if (!named_list) names(classes) <- NULL
      weigth_collection[[j]] <- weight_list
    }
    classes <- weight_collection
  }
  
  classes
}

#' Count nucleotides per class
#'
#' @inheritParams get_class_weight
#' @noRd
count_nuc <- function(path,
                      vocabulary_label = NULL,
                      format = "fasta",
                      # estimate class distribution from subset
                      file_proportion = 1,
                      train_type = "label_folder",
                      csv_path = NULL) {
  
  classes <- rep(0, length(vocabulary_label))
  names(classes) <- vocabulary_label
  
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
        classes[j] <- classes[j] + freq
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
      df <- stats::aggregate(df$freq, by = list(Category = df$Header), FUN = sum)
      freq <- df$x
      names(freq) <- df$Category
      for (k in names(freq)) {
        classes[k] <- classes[k] + freq[k]
      }
    }
  }
  
  # label csv
  if (train_type == "label_csv") {
    
    label_csv <- utils::read.csv2(csv_path, header = TRUE, stringsAsFactors = FALSE)
    if (dim(label_csv)[2] == 1) {
      label_csv <- utils::read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)
    }
    if (!("file" %in% names(label_csv))) {
      stop('csv file needs one column named "file"')
    }
    
    row_sums <- label_csv %>% dplyr::select(-file) %>% rowSums()
    if (!(all(row_sums == 1))) {
      stop("Can only estimate class weights if labels are mutually exclusive.")
    }
    
    if (is.null(vocabulary_label) || missing(vocabulary_label)) {
      vocabulary_label <-  names(label_csv)[!names(label_csv) == "file"]
    } else {
      label_csv <- label_csv %>% dplyr::select(c(dplyr::all_of(vocabulary_label), "file"))
    }
    
    classes <- rep(0, length(vocabulary_label))
    names(classes) <- vocabulary_label
    
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

read_fasta_fastq <- function(format, skip_amb_nuc, file_index, pattern, shuffle_input,
                             reverse_complement, fasta.files, use_coverage = FALSE, proportion_entries = NULL,
                             vocabulary_label = NULL, filter_header = FALSE, target_from_csv = NULL) {
  
  if (stringr::str_detect(format, "fasta")) {
    if (is.null(skip_amb_nuc)) {
      fasta.file <- microseq::readFasta(fasta.files[file_index])
    } else {
      fasta.file <- remove_amb_nuc_entries(microseq::readFasta(fasta.files[file_index]), skip_amb_nuc = skip_amb_nuc,
                                           pattern = pattern)
    }
    
    if (filter_header & is.null(target_from_csv)) {
      label_vector <- trimws(stringr::str_to_lower(fasta.file$Header))
      label_filter <- label_vector %in% vocabulary_label
      fasta.file <- fasta.file[label_filter, ]
    }
    
    if (!is.null(proportion_entries) && proportion_entries < 1) {
      index <- sample(nrow(fasta.file), max(1, floor(nrow(fasta.file) * proportion_entries)))
      fasta.file <- fasta.file[index, ]
    }
    
    if (shuffle_input) {
      fasta.file <- fasta.file[sample(nrow(fasta.file)), ]
    }
    
    if (reverse_complement) {
      index <- sample(c(TRUE, FALSE), nrow(fasta.file), replace = TRUE)
      fasta.file$Sequence[index] <- microseq::reverseComplement(fasta.file$Sequence[index])
    }
    
  }
  
  if (stringr::str_detect(format, "fastq")) {
    if (is.null(skip_amb_nuc)) {
      fasta.file <- microseq::readFastq(fasta.files[file_index])
    } else {
      fasta.file <- remove_amb_nuc_entries(microseq::readFastq(fasta.files[file_index]), skip_amb_nuc = skip_amb_nuc,
                                           pattern = pattern)
    }
    
    if (filter_header & is.null(target_from_csv)) {
      label_vector <- trimws(stringr::str_to_lower(fasta.file$Header))
      label_filter <- label_vector %in% vocabulary_label
      fasta.file <- fasta.file[label_filter, ]
    }
    
    if (!is.null(proportion_entries) && proportion_entries < 1) {
      index <- sample(nrow(fasta.file), max(1, floor(nrow(fasta.file) * proportion_entries)))
      fasta.file <- fasta.file[index, ]
    }
    
    if (shuffle_input) {
      fasta.file <- fasta.file[sample(nrow(fasta.file)), ]
    }
    
    if (reverse_complement & sample(c(TRUE, FALSE), 1)) {
      fasta.file$Sequence <- microseq::reverseComplement(fasta.file$Sequence)
    }
  }
  return(fasta.file)
}

input_from_csv <- function(added_label_path) {
  .datatable.aware = TRUE
  label_csv <- utils::read.csv2(added_label_path, header = TRUE, stringsAsFactors = FALSE)
  if (dim(label_csv)[2] == 1) {
    label_csv <- utils::read.csv(added_label_path, header = TRUE, stringsAsFactors = FALSE)
  }
  label_csv <- data.table::as.data.table(label_csv)
  label_csv$file <- stringr::str_to_lower(as.character(label_csv$file))
  data.table::setkey(label_csv, file)
  added_label_by_header <- FALSE
  
  if (!("file" %in% names(label_csv))) {
    stop('names in added_label_path should contain one column named "file" ')
  }
  col_name <- ifelse(added_label_by_header, "header", "file")
  return(list(label_csv = label_csv, col_name = col_name))
}

#' @rawNamespace import(data.table, except = c(first, last, between))
#' @noRd
csv_to_tensor <- function(label_csv, added_label_vector, added_label_by_header, batch_size,
                          start_index_list) {
  .datatable.aware = TRUE
  label_tensor <- matrix(0, ncol = ncol(label_csv) - 1, nrow = batch_size, byrow = TRUE)
  
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
      
      if (nrow(stats::na.omit(label_from_csv)) > 0) {
        label_tensor[assign_rows, ] <- matrix(as.matrix(label_from_csv[1, ]),
                                              nrow = samples_per_file, ncol = ncol(label_tensor), byrow = TRUE)
      }
      row_index <- row_index + samples_per_file
    }
  }
  return(label_tensor)
}

#' Divide tensor to list of subsets
#'
#' @noRd
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

check_header_names <- function(target_split, vocabulary_label) {
  target_split <- unlist(target_split)
  if (!all(target_split %in% vocabulary_label)) {
    stop_text <- paste("Your csv file has no columns named",
                       paste(target_split[!(target_split %in% vocabulary_label)], collapse = " "))
    stop(stop_text)
  }
  if (!all(vocabulary_label %in% target_split)) {
    warning_text <- paste("target_split does not cover the following columns:",
                          paste(vocabulary_label[!(vocabulary_label %in% target_split)], collapse = " "))
    warning(warning_text)
  }
}

count_files <- function(path, format = "fasta", train_type,
                        target_from_csv = NULL, train_val_split_csv = NULL) {
  
  num_files <- rep(0, length(path))
  if (!is.null(target_from_csv) & train_type == "label_csv") {
    target_files <- utils::read.csv(target_from_csv)
    if (ncol(target_files) == 1) target_files <- utils::read.csv2(target_from_csv)
    target_files <- target_files$file
    # are files given with absolute path
    full.names <- ifelse(dirname(target_files[1]) == ".", FALSE, TRUE) 
  }  
  if (!is.null(train_val_split_csv)) {
    tvt_files <- utils::read.csv(train_val_split_csv)
    if (ncol(tvt_files) == 1) tvt_files <- utils::read.csv2(train_val_split_csv)
    train_index <- tvt_files$type == "train"
    tvt_files <- tvt_files$file
    target_files <- intersect(tvt_files[train_index], target_files)
  }  
  
  for (i in 1:length(path)) {
    for (k in 1:length(path[[i]])) {
      current_path <- path[[i]][[k]]
      
      if (!is.null(train_val_split_csv)) {
        if (!(current_path %in% target_files)) next
      }
      
      if (endsWith(current_path, paste0(".", format))) {
        # remove files not in csv file 
        if (!is.null(target_from_csv)) {
          current_files <- length(intersect(basename(target_files), basename(current_path)))
        } else {
          current_files <- 1
        }
      } else {
        # remove files not in csv file 
        if (!is.null(target_from_csv)) {
          current_files <- list.files(current_path, pattern = paste0(".", format, "$"), full.names = full.names) %>%
            intersect(target_files) %>% length()
        } else {
          current_files <- list.files(current_path, pattern = paste0(".", format, "$")) %>% length()
        }
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

list_fasta_files <- function(path_corpus, format, file_filter) {
  
  fasta.files <- list()
  path_corpus <- unlist(path_corpus)
  
  for (i in 1:length(path_corpus)) {
    
    if (endsWith(path_corpus[[i]], paste0(".", format))) {
      fasta.files[[i]] <- path_corpus[[i]]
      
    } else {
      
      fasta.files[[i]] <- list.files(
        path = path_corpus[[i]],
        pattern = paste0("\\.", format, "$"),
        full.names = TRUE)
    }
  }
  fasta.files <- unlist(fasta.files)
  num_files <- length(fasta.files)
  
  if (!is.null(file_filter)) {
    
    # file filter files given with/without absolute path
    if (all(basename(file_filter) == file_filter)) {
      fasta.files <- fasta.files[basename(fasta.files) %in% file_filter]
    } else {
      fasta.files <- fasta.files[fasta.files %in% file_filter]
    }
    
    if (length(fasta.files) < 1) {
      stop_text <- paste0("None of the files from ", unlist(path_corpus),
                          " are present in train_val_split_csv table for either train or validation. \n")
      stop(stop_text)
    }
  }
  
  fasta.files <- gsub(pattern="/+", replacement="/", x = fasta.files)
  fasta.files <- gsub(pattern="/$", replacement="", x = fasta.files)
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

#' Reshape tensors for set learning
#' 
#' Reshape input x and target y. Aggregates multiple samples from x and y into single input/target batches.  
#' 
#' @param x 3D input tensor.
#' @param y 2D target tensor.
#' @param samples_per_target How many samples to use for one target
#' @param reshape_mode `"time_dist", "multi_input"` or `"concat"` 
#' \itemize{
#' \item If `"multi_input"`, will produce `samples_per_target` separate inputs, each of length `maxlen`.
#' \item If `"time_dist"`, will produce a 4D input array. The dimensions correspond to
#' `(new_batch_size, samples_per_target, maxlen, length(vocabulary))`.
#' \item If `"concat"`, will concatenate `samples_per_target` sequences of length `maxlen` to one long sequence
#' }
#' @param buffer_len Only applies if `reshape_mode = "concat"`. If `buffer_len` is an integer, the subsequences are interspaced with `buffer_len` rows. The reshaped x has
#' new maxlen: (`maxlen` \eqn{*} `samples_per_target`) + `buffer_len` \eqn{*} (`samples_per_target` - 1).
#' @param new_batch_size Size of first axis of input/targets after reshaping.
#' @param check_y Check if entries in `y` are consistent with reshape strategy (same label when aggregating).   
#' @examplesIf reticulate::py_module_available("tensorflow")
#' # create dummy data
#' batch_size <- 8
#' maxlen <- 11
#' voc_len <- 4 
#' x <- sample(0:(voc_len-1), maxlen*batch_size, replace = TRUE)
#' x <- keras::to_categorical(x, num_classes = voc_len)
#' x <- array(x, dim = c(batch_size, maxlen, voc_len))
#' y <- rep(0:1, each = batch_size/2)
#' y <- keras::to_categorical(y, num_classes = 2)
#' y
#' 
#' # reshape data for multi input model
#' reshaped_data <- reshape_tensor(
#'   x = x,
#'   y = y,
#'   new_batch_size = 2,
#'   samples_per_target = 4,
#'   reshape_mode = "multi_input")
#' 
#' length(reshaped_data[[1]])
#' dim(reshaped_data[[1]][[1]])
#' reshaped_data[[2]]
#' 
#' @returns A list of 2 tensors.
#' @export
reshape_tensor <- function(x, y, new_batch_size,
                           samples_per_target,
                           buffer_len = NULL,
                           reshape_mode = "time_dist",
                           check_y = FALSE) {
  
  batch_size <- dim(x)[1]
  maxlen <- dim(x)[2]
  voc_len <- dim(x)[3]
  num_classes <- dim(y)[2]
  
  if (check_y) {
    targets <- apply(y, 1, which.max)
    test_y_dist <- all(targets == rep(1:num_classes, each = batch_size/num_classes))
    if (!test_y_dist) {
      stop("y must have same number of samples for each class")
    }
  }
  
  if (reshape_mode == "time_dist") {
    
    x_new <- array(0, dim = c(new_batch_size, samples_per_target, maxlen, voc_len))
    y_new <- array(0, dim = c(new_batch_size, num_classes))
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
      x_index <- base::seq(i, batch_size, samples_per_target)
      x_list[[i]] <- x[x_index, , ]
    }
    y <- y[base::seq(1, batch_size, samples_per_target), ]
    return(list(x = x_list, y = y))
  }
  
  if (reshape_mode == "concat") {
    
    use_buffer <- !is.null(buffer_len) && buffer_len > 0
    if (use_buffer) {
      buffer_tensor <- array(0, dim = c(buffer_len, voc_len))
      buffer_tensor[ , voc_len] <- 1
      concat_maxlen <- (maxlen * samples_per_target) + (buffer_len * (samples_per_target - 1))
    } else {
      concat_maxlen <- maxlen * samples_per_target 
    }
    
    x_new <- array(0, dim = c(new_batch_size, concat_maxlen, voc_len))
    y_new <- array(0, dim = c(new_batch_size, num_classes))
    
    
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
      }
      
      x_new[i, , ] <- x_temp
      y_new[i, ]  <- y[index[1], ]
    }
    return(list(x = x_new, y = y_new))
  }
}

#' Transform confusion matrix with total numbers to matrix with percentages.
#'
#' @noRd
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
  cm
}

create_conf_mat_obj <- function(m, confMatLabels) {
  dimnames(m) <- list(Prediction = confMatLabels, Truth = confMatLabels)
  l <- list()
  m <- as.table(m)
  l[["table"]] <- m
  l[["dots"]] <- list()
  class(l) <- "conf_mat"
  return(l)
}

#' Encode sequence of integers to sequence of n-gram 
#' 
#' Input is sequence of integers from vocabulary of size \code{voc_size}. 
#' Returns vector of integers corresponding to n-gram encoding.
#' Integers greater than `voc_size` get encoded as `voc_size^n + 1`.
#' 
#' @param int_seq Integer sequence
#' @param n Length of n-gram aggregation
#' @param voc_size Size of vocabulary.
#' @examples
#' int_to_n_gram(int_seq = c(1,1,2,4,4), n = 2, voc_size = 4)
#' 
#' @returns A numeric vector.
#' @export
int_to_n_gram <- function(int_seq, n, voc_size = 4) {
  
  encoding_len <- length(int_seq) - n + 1
  n_gram_encoding <- vector("numeric", encoding_len)
  oov_token <- voc_size^n + 1
  padding_token <- 0
  
  for (i in 1:encoding_len) {
    int_seq_subset <- int_seq[i:(i + n - 1)]
    
    if (prod(int_seq_subset) == 0) {
      n_gram_encoding[i] <- padding_token
    } else {
      # encoding for amb nuc
      if (any(int_seq_subset > voc_size)) {
        n_gram_encoding[i] <- oov_token
      } else {
        int_seq_subset <- int_seq_subset - 1
        n_gram_encoding[i] <- 1 + sum(voc_size^((n-1):0) * (int_seq_subset))
      }
    }
  }
  n_gram_encoding
}

#' One-hot encoding matrix to n-gram encoding matrix
#' 
#' @param input_matrix Matrix with one 1 per row and zeros otherwise.
#' @param n Length of one n-gram.   
#' @examplesIf reticulate::py_module_available("tensorflow")
#' x <- c(0,0,1,3,3) 
#' input_matrix <- keras::to_categorical(x, 4)
#' n_gram_of_matrix(input_matrix, n = 2) 
#' 
#' @returns Matrix of one-hot encodings. 
#' @export
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

n_gram_vocabulary <- function(n_gram = 3, vocabulary = c("A", "C", "G", "T")) {
  l <- list()
  for (i in 1:n_gram) {
    l[[i]] <- vocabulary
  }
  df <- expand.grid(l)
  df <- df[ , ncol(df) : 1]
  n_gram_nuc <- apply(df, 1, paste, collapse = "") 
  n_gram_nuc
}


#' Split fasta file into smaller files.
#'
#' Returns smaller files with same file name and "_x" (where x is an integer). For example,
#' assume we have input file called "abc.fasta" with 100 entries and `split_n = 50`. Function will
#' create two files called "abc_1.fasta" and "abc_2.fasta" in `target_path`.
#'
#' @param path_input Fasta file to split into smaller files
#' @param split_n Maximum number of entries to use in smaller file.
#' @param target_folder Directory for output.
#' @param shuffle_entries Whether to shuffle fasta entries before split.
#' @param delete_input Whether to delete the original file.
#' @examples
#' path_input <- tempfile(fileext = '.fasta')
#' create_dummy_data(file_path = path_input,
#'                   num_files = 1,
#'                   write_to_file_path = TRUE,
#'                   seq_length = 7,
#'                   num_seq = 25,
#'                   vocabulary = c("a", "c", "g", "t"))
#' target_folder <- tempfile()
#' dir.create(target_folder)
#' 
#' # split 25 entries into 5 files
#' split_fasta(path_input = path_input,
#'             target_folder = target_folder,
#'             split_n = 5)
#' length(list.files(target_folder)) 
#' 
#' @returns None. Writes files to output.
#' @export
split_fasta <- function(path_input,
                        target_folder,
                        split_n = 500,
                        shuffle_entries = TRUE,
                        delete_input = FALSE) {
  
  fasta_file <- microseq::readFasta(path_input)
  
  base_name <- basename(stringr::str_remove(path_input, ".fasta"))
  new_path <- paste0(target_folder, "/", base_name)
  count <- 1
  start_index <- 1
  end_index <- 1
  
  if (nrow(fasta_file) == 1) {
    fasta_name <- paste0(new_path, "_", count, ".fasta")
    microseq::writeFasta(fasta_file, fasta_name)
    if (delete_input) {
      file.remove(path_input)
    }
    return(NULL)
  }
  
  if (shuffle_entries) {
    fasta_file <- fasta_file[sample(nrow(fasta_file)), ]
  }
  
  while (end_index < nrow(fasta_file)) {
    end_index <- min(start_index + split_n - 1, nrow(fasta_file))
    index <- start_index : end_index
    sub_df <- fasta_file[index, ]
    fasta_name <- paste0(new_path, "_", count, ".fasta")
    microseq::writeFasta(sub_df, fasta_name)
    start_index <- start_index + split_n
    count <- count + 1
  }
  
  if (delete_input) {
    file.remove(path_input)
  }
}

#' Add noise to tensor
#'
#' @param noise_type "normal" or "uniform".
#' @param ... additional arguments for rnorm or runif call.
#' @noRd
add_noise_tensor <- function(x, noise_type, ...) {
  
  stopifnot(noise_type %in% c("normal", "uniform"))
  random_fn <- ifelse(noise_type == "normal", "rnorm", "runif")
  
  if (is.list(x)) {
    for (i in 1:length(x)) {
      x_dim <- dim(x[[i]])
      noise_tensor <- do.call(random_fn, list(n = prod(x_dim[-1]), ...))
      noise_tensor <- array(noise_tensor, dim = x_dim)
      x[[i]] <- x[[i]] + noise_tensor
    }
  } else {
    x_dim <- dim(x)
    stopifnot(noise_type %in% c("normal", "uniform"))
    random_fn <- ifelse(noise_type == "normal", "rnorm", "runif")
    noise_tensor <- do.call(random_fn, list(n = prod(x_dim[-1]), ...))
    noise_tensor <- array(noise_tensor, dim = x_dim)
    x <- x + noise_tensor
  }
  
  return(x)
}

reverse_complement_tensor <- function(x) {
  stopifnot(dim(x)[3] == 4)
  x_rev_comp <- x[ ,  dim(x)[2]:1, 4:1]
  x_rev_comp <- array(x_rev_comp, dim = dim(x))
  x_rev_comp
}


get_pos_enc <- function(pos, i, d_model, n = 10000) {
  
  pw <- (2 * floor(i/2)) / d_model
  angle_rates <- 1 / (n ^ pw)
  angle <- pos * angle_rates
  pos_enc <- ifelse(i %% 2 == 0, sin(angle), cos(angle))
  return(pos_enc)
}  

positional_encoding <- function(seq_len, d_model, n=10000) {
  
  P = matrix(0, nrow = seq_len, ncol = d_model)
  
  for (pos in 0:(seq_len - 1)) {
    for (i in 0:(d_model - 1)) {
      P[pos + 1, i + 1] <- get_pos_enc(pos, i, d_model, n)
    }
  }
  
  return(P)
}


subset_tensor_list <- function(tensor_list, dim_list, subset_index, dim_n_list) {
  
  for (i in 1:length(tensor_list)) {
    tensor_list[[i]] <- subset_tensor(tensor = tensor_list[[i]],
                                      subset_index = subset_index,
                                      dim_n = dim_n_list[[i]])
  }
  
}

subset_tensor <- function(tensor, subset_index, dim_n) {
  
  if (dim_n == 1) {
    subset_tensor <- tensor[subset_index]
  }
  
  if (dim_n == 2) {
    subset_tensor <- tensor[subset_index, ]
  }
  
  if (dim_n == 3) {
    subset_tensor <- tensor[subset_index, , ]
  }
  
  if (dim_n == 4) {
    subset_tensor <- tensor[subset_index, , , ]
  }
  
  if (length(subset_index) == 1 & dim_n > 1) {
    subset_tensor <- tensorflow::tf$expand_dims(subset_tensor, axis = 0L)
  }
}


mask_seq <- function(int_seq,
                     mask_rate = NULL,
                     random_rate = NULL,
                     identity_rate = NULL,
                     block_len = NULL,
                     start_ind = NULL,
                     voc_len) {
  
  mask_token <- voc_len + 1
  if (is.null(mask_rate)) mask_rate <- 0
  if (is.null(random_rate)) random_rate <- 0
  if (is.null(identity_rate)) identity_rate <- 0
  mask_perc <- mask_rate + random_rate + identity_rate
  if (mask_perc > 1) {
    stop("Sum of mask_rate, random_rate, identity_rate bigger than 1")
  } 
  # don't mask padding or oov positions 
  valid_pos <- which(int_seq != 0 & int_seq != mask_token) 
  
  # randomly decide whether to round up or down
  ceiling_floor <- sample(c(TRUE, FALSE), 3, replace = TRUE)
  # adjust for block len
  block_len_adjust <- ifelse(is.null(block_len), 1, block_len) 
  
  num_mask_pos <- (mask_rate * length(valid_pos))/block_len_adjust
  num_mask_pos <- ifelse(ceiling_floor[1], floor(num_mask_pos), ceiling(num_mask_pos))
  num_random_pos <- (random_rate * length(valid_pos))/block_len_adjust
  num_random_pos <- ifelse(ceiling_floor[2], floor(num_random_pos), ceiling(num_random_pos))
  num_identity_pos <- (identity_rate * length(valid_pos))/block_len_adjust
  num_identity_pos <- ifelse(ceiling_floor[3], floor(num_identity_pos), ceiling(num_identity_pos))
  num_all_pos <- num_mask_pos + num_random_pos + num_identity_pos
  if (is.null(block_len)) {
    all_pos <- sample(valid_pos, num_all_pos)
  } else {
    valid_pos_block_len <- seq(from = sample(1:(block_len - 1), 1), to = length(valid_pos), by = block_len)
    valid_pos <- intersect(valid_pos_block_len, valid_pos)
    all_pos <- sample(valid_pos, min(num_all_pos, length(valid_pos)))
  }
  
  sample_weight_seq <- rep(0, length(int_seq))
  if (is.null(block_len)) {
    sample_weight_seq[all_pos] <- 1
  } else {
    all_pos_blocks <- purrr::map(all_pos, ~seq(.x, .x + block_len - 1, by = 1))
    sample_weight_seq[unlist(all_pos_blocks)] <- 1
  }
  
  if (num_mask_pos > 0) {
    mask_index <- sample(all_pos, num_mask_pos)
    all_pos <- setdiff(all_pos, mask_index)
    if (!is.null(block_len)) {
      mask_index <- purrr::map(mask_index, ~seq(.x, .x + block_len - 1, by = 1)) %>% 
        unlist()
    }
    int_seq[mask_index] <- mask_token
  }
  
  if (num_random_pos > 0) {
    random_index <- sample(all_pos, num_random_pos)
    all_pos <- setdiff(all_pos, random_index)
    if (!is.null(block_len)) {
      random_index <- purrr::map(random_index, ~seq(.x, .x + block_len - 1, by = 1)) %>% 
        unlist()
    }
    int_seq[random_index] <- sample(1:voc_len, length(random_index), replace = TRUE)
  }
  
  # mask oov tokens
  sample_weight_seq[int_seq == mask_token] <- 1
  
  return(list(masked_seq = int_seq, sample_weight_seq = sample_weight_seq))
  
}

#' Char sequence corresponding to one-hot matrix.
#'
#' Return character sequence corresponding to one-hot elements in matrix or tensor.
#'
#' @inheritParams generator_fasta_lm
#' @param m One-hot encoding matrix or 3d array where each element of first axis is one-hot matrix.
#' @param amb_enc Either `"zero"` or `"equal"`. How oov tokens where treated for one-hot encoding. 
#' @param amb_char Char to use for oov positions.
#' @param paste_chars Whether to return vector or single sequence.
#' @examples 
#' m <- matrix(c(1,0,0,0,0,1,0,0), 2)
#' one_hot_to_seq(m)
#' 
#' @returns A string.
#' @export
one_hot_to_seq <- function(m, vocabulary = c("A", "C", "G", "T"), amb_enc = "zero",
                           amb_char = "N", paste_chars = TRUE) {
  
  if (length(dim(m)) == 3) {
    seq_list <- list()
    for (i in 1:dim(m)[1]) {
      seq_list[[i]] <- one_hot_to_seq(m = m[i, , ], vocabulary = vocabulary, amb_enc = amb_enc,
                                      amb_char = amb_char, paste_chars = paste_chars)
    }
    return(seq_list)
  }
  
  if (amb_enc == "zero") {
    amb_row <- which(rowSums(m) == 0)
  }
  
  if (amb_enc == "equal") {
    amb_row <- which(rowSums[ , 1] == 1/length(vocabulary))
  }
  
  nt_seq <- vocabulary[apply(m, 1, which.max)]
  nt_seq[amb_row] <- amb_char
  
  if (paste_chars) {
    nt_seq <- paste(nt_seq, collapse = "")
  } 
  
  return(nt_seq)
  
}
