#' Language model generator for fasta/fastq files
#'
#' @description Iterates over folder containing fasta/fastq files and produces encoding of predictor sequences
#' and target variables. Will take a sequence of fixed size and use some part of sequence as input and other part as target. 
#'
#' @inheritParams train_model
#' @param path_corpus Input directory where fasta files are located or path to single file ending with fasta or fastq
#' (as specified in format argument). Can also be a list of directories and/or files.
#' @param format File format, either `"fasta"` or `"fastq"`.
#' @param batch_size Number of samples in one batch.
#' @param maxlen Length of predictor sequence.
#' @param max_iter Stop after `max_iter` number of iterations failed to produce a new batch.
#' @param shuffle_file_order Logical, whether to go through files randomly or sequentially.
#' @param step How often to take a sample.
#' @param seed Sets seed for `set.seed` function for reproducible results.
#' @param shuffle_input Whether to shuffle entries in every fasta/fastq file before extracting samples.
#' @param verbose Whether to show messages.
#' @param path_file_log Write name of files to csv file if path is specified.
#' @param reverse_complement Boolean, for every new file decide randomly to use original data or its reverse complement.
#' @param ambiguous_nuc How to handle nucleotides outside vocabulary, either `"zero"`, `"discard"`, `"empirical"` or `"equal"`.
#' \itemize{
#' \item If `"zero"`, input gets encoded as zero vector.
#' \item If `"equal"`, input is repetition of `1/length(vocabulary)`.
#' \item If `"discard"`, samples containing nucleotides outside vocabulary get discarded.
#' \item If `"empirical"`, use nucleotide distribution of current file.
#' }
#' @param proportion_per_seq Numerical value between 0 and 1. Proportion of sequence to take samples from (use random subsequence).
#' @param use_quality_score Whether to use fastq quality scores. If TRUE input is not one-hot-encoding but corresponds to probabilities.
#' For example (0.97, 0.01, 0.01, 0.01) instead of (1, 0, 0, 0).
#' @param padding Whether to pad sequences too short for one sample with zeros.
#' @param added_label_path Path to file with additional input labels. Should be a csv file with one column named "file". Other columns should correspond to labels.
#' @param add_input_as_seq Boolean vector specifying for each entry in \code{added_label_path} if rows from csv should be encoded as a sequence or used directly.
#' If a row in your csv file is a sequence this should be `TRUE`. For example you may want to add another sequence, say ACCGT. Then this would correspond to 1,2,2,3,4 in
#' csv file (if vocabulary = c("A", "C", "G", "T")).  If \code{add_input_as_seq} is `TRUE`, 12234 gets one-hot encoded, so added input is a 3D tensor.  If \code{add_input_as_seq} is
#' `FALSE` this will feed network just raw data (a 2D tensor).
#' @param skip_amb_nuc Threshold of ambiguous nucleotides to accept in fasta entry. Complete entry will get discarded otherwise.
#' @param max_samples Maximum number of samples to use from one file. If not `NULL` and file has more than \code{max_samples} samples, will randomly choose a
#' subset of \code{max_samples} samples.
#' @param concat_seq Character string or `NULL`. If not `NULL` all entries from file get concatenated to one sequence with `concat_seq` string between them.
#' Example: If 1.entry AACC, 2. entry TTTG and `concat_seq = "ZZZ"` this becomes AACCZZZTTTG.
#' @param target_len Number of nucleotides to predict at once for language model.
#' @param file_filter Vector of file names to use from path_corpus.
#' @param use_coverage Integer or `NULL`. If not `NULL`, use coverage as encoding rather than one-hot encoding and normalize.
#' Coverage information must be contained in fasta header: there must be a string `"cov_n"` in the header, where `n` is some integer.
#' @param proportion_entries Proportion of fasta entries to keep. For example, if fasta file has 50 entries and `proportion_entries = 0.1`,
#' will randomly select 5 entries.
#' @param sample_by_file_size Sample new file weighted by file size (bigger files more likely).
#' @param n_gram Integer, encode target not nucleotide wise but combine n nucleotides at once. For example for `n=2, "AA" ->  (1, 0,..., 0),`
#' `"AC" ->  (0, 1, 0,..., 0), "TT" -> (0,..., 0, 1)`, where the one-hot vectors have length `length(vocabulary)^n`.
#' @param add_noise `NULL` or list of arguments. If not `NULL`, list must contain the following arguments: \code{noise_type} can be `"normal"` or `"uniform"`;
#' optional arguments `sd` or `mean` if noise_type is `"normal"` (default is `sd=1` and `mean=0`) or `min, max` if `noise_type` is `"uniform"`
#' (default is `min=0, max=1`).
#' @import data.table
#' @importFrom magrittr %>%
#' @export
generator_fasta_lm <- function(path_corpus,
                               format = "fasta",
                               batch_size = 256,
                               maxlen = 250,
                               max_iter = 10000,
                               vocabulary = c("a", "c", "g", "t"),
                               verbose = FALSE,
                               shuffle_file_order = FALSE,
                               step = 1,
                               seed = 1234,
                               shuffle_input = FALSE,
                               file_limit = NULL,
                               path_file_log = NULL,
                               reverse_complement = FALSE,
                               output_format = "target_right",
                               ambiguous_nuc = "zeros",
                               use_quality_score = FALSE,
                               proportion_per_seq = NULL,
                               padding = TRUE,
                               added_label_path = NULL,
                               add_input_as_seq = NULL,
                               skip_amb_nuc = NULL,
                               max_samples = NULL,
                               concat_seq = NULL,
                               target_len = 1,
                               file_filter = NULL,
                               use_coverage = NULL,
                               proportion_entries = NULL,
                               sample_by_file_size = FALSE,
                               n_gram = NULL,
                               n_gram_stride = 1,
                               add_noise = NULL,
                               return_int = FALSE) {
  
  
  ##TODO: add check for n-gram
  
  total_seq_len <- maxlen + target_len
  gen <- generator_fasta_label_folder(path_corpus = path_corpus,
                                      format = format,
                                      batch_size = batch_size,
                                      maxlen = total_seq_len,
                                      max_iter = max_iter,
                                      vocabulary = vocabulary,
                                      shuffle_file_order = shuffle_file_order,
                                      step = step,
                                      seed = seed,
                                      shuffle_input = shuffle_input,
                                      file_limit = file_limit,
                                      path_file_log = path_file_log,
                                      reverse_complement = reverse_complement,
                                      reverse_complement_encoding = FALSE,
                                      num_targets = 1,
                                      ones_column = 1,
                                      ambiguous_nuc = ambiguous_nuc,
                                      proportion_per_seq = proportion_per_seq,
                                      read_data = FALSE,
                                      use_quality_score = use_quality_score,
                                      padding = padding,
                                      added_label_path = added_label_path,
                                      add_input_as_seq = add_input_as_seq,
                                      skip_amb_nuc = skip_amb_nuc,
                                      max_samples = max_samples,
                                      concat_seq = concat_seq,
                                      file_filter = file_filter,
                                      use_coverage = use_coverage,
                                      proportion_entries = proportion_entries,
                                      sample_by_file_size = sample_by_file_size,
                                      n_gram = n_gram,
                                      n_gram_stride = n_gram_stride,
                                      masked_lm = NULL,
                                      add_noise = add_noise,
                                      return_int = return_int)
  
  function() {
    
    if (is.null(added_label_path)) {
      xy <- gen()[[1]]
    } else {
      z <- gen()[[1]]
      added_input <- z[1:(length(z)-1)]
      xy <- z[length(z)][[1]]
    }
    
    xy_list <- slice_tensor_lm(xy = xy,
                               output_format = output_format,
                               target_len = target_len,
                               n_gram = n_gram,
                               total_seq_len = total_seq_len,
                               return_int = return_int)
    
    # if (batch_size == 1) {
    #   xy_list$x <- add_dim(xy_list$x)
    #   xy_list$y <- add_dim(xy_list$y)
    # }
    
    if (is.null(added_label_path)) {
      return(xy_list)
    } else {
      return(list(append(added_input, list(xy_list$x)), xy_list$y))
    }
    
    # add dim for batch size 1
    
  }
}

#' Data generator for fasta/fastq files and label targets
#'
#' @description Iterates over folder containing fasta/fastq files and produces encoding of predictor sequences
#' and target variables. Targets will be read from fasta headers or a separate csv file.
#'
#' @inheritParams generator_fasta_lm
#' @inheritParams train_model
#' @param vocabulary_label Character vector of possible targets. Targets outside \code{vocabulary_label} will get discarded.
#' @param target_from_csv Path to csv file with target mapping. One column should be called "file" and other entries in row are the targets.
#' @param target_split If target gets read from csv file, list of names to divide target tensor into list of tensors.
#' Example: if csv file has header names `"file", "label_1", "label_2", "label_3"` and `target_split = list(c("label_1", "label_2"), "label_3")`,
#' this will divide target matrix to list of length 2, where the first element contains columns named `"label_1"` and `"label_2"` and the
#' second entry contains the column named `"label_3"`.
#' @param read_data If `TRUE` the first element of output is a list of length 2, each containing one part of paired read. Maxlen should be 2*length of one read.
#' @import data.table
#' @export
generator_fasta_label_header_csv <- function(path_corpus,
                                             format = "fasta",
                                             batch_size = 256,
                                             maxlen = 250,
                                             max_iter = 10000,
                                             vocabulary = c("a", "c", "g", "t"),
                                             verbose = FALSE,
                                             shuffle_file_order = FALSE,
                                             step = 1,
                                             seed = 1234,
                                             shuffle_input = FALSE,
                                             file_limit = NULL,
                                             path_file_log = NULL,
                                             vocabulary_label = c("x", "y", "z"),
                                             reverse_complement = TRUE,
                                             ambiguous_nuc = "zero",
                                             proportion_per_seq = NULL,
                                             read_data = FALSE,
                                             use_quality_score = FALSE,
                                             padding = TRUE,
                                             skip_amb_nuc = NULL,
                                             max_samples = NULL,
                                             concat_seq = NULL,
                                             added_label_path = NULL,
                                             add_input_as_seq = NULL,
                                             target_from_csv = NULL,
                                             target_split = NULL,
                                             file_filter = NULL,
                                             use_coverage = NULL,
                                             proportion_entries = NULL,
                                             sample_by_file_size = FALSE,
                                             reverse_complement_encoding = FALSE,
                                             n_gram = NULL,
                                             n_gram_stride = 1,
                                             add_noise = NULL,
                                             return_int = FALSE) {
  
  if (is.null(use_coverage)) {
    use_coverage <- FALSE
    cov_vector <- NULL
    max_cov <- NULL
  } else {
    max_cov <- use_coverage
    use_coverage <- TRUE
  }
  if (!is.null(concat_seq) && (!all(stringr::str_split(concat_seq,"")[[1]] %in% vocabulary))) {
    stop("Characters of separating sequence should be in vocabulary")
  }
  if (reverse_complement_encoding) {
    test_len <- length(vocabulary) != 4
    if (test_len || all(sort(stringr::str_to_lower(vocabulary)) != c("a", "c", "g", "t"))) {
      stop("reverse_complement_encoding only implemented for A,C,G,T vocabulary yet")
    }
  }
  
  discard_amb_nuc <- ifelse(ambiguous_nuc == "discard", TRUE, FALSE)
  vocabulary <- stringr::str_to_lower(vocabulary)
  vocabulary_label <- stringr::str_to_lower(vocabulary_label)
  start_index_list <- vector("list")
  file_index <- 1
  num_samples <- 0
  start_index <- 1
  iter <- 1
  concat <- !is.null(concat_seq)
  if (concat & is.null(target_from_csv)) {
    stop("Cannot concatenate fasta sequences when reading label from header")
  }
  additional_labels <- !is.null(added_label_path)
  seq_vector <- NULL
  
  # # adjust maxlen for n_gram
  # if (!is.null(n_gram)) {
  #   stop("n-gram encoding not implemented yet for classification")
  #   maxlen <- maxlen + n_gram - 1
  # }
  
  for (i in letters) {
    if (!(i %in% stringr::str_to_lower(vocabulary))) {
      amb_nuc_token <- i
      break
    }
  }
  tokenizer_pred <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), c(vocabulary, amb_nuc_token))
  tokenizer_target <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = FALSE, lower = TRUE, filters = "\t\n"),
                                                vocabulary_label)
  
  fasta.files <-  list_fasta_files(path_corpus = path_corpus,
                                   format = format,
                                   file_filter = file_filter)
  num_files <- length(fasta.files)
  
  if (sample_by_file_size) {
    shuffle_file_order <- FALSE
    file_prob <- file.info(fasta.files)$size/sum(file.info(fasta.files)$size)
  }
  
  set.seed(seed)
  if (shuffle_file_order) fasta.files <- sample(fasta.files, replace = FALSE)
  
  # target from csv
  if (!is.null(target_from_csv)) {
    .datatable.aware = TRUE
    output_label_csv <- read.csv2(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
    if (dim(output_label_csv)[2] == 1) {
      output_label_csv <- read.csv(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
    }
    output_label_csv <- data.table::as.data.table(output_label_csv)
    if ("file" %in% names(output_label_csv) & "header" %in% names(output_label_csv)) {
      stop('names in target_from_csv should contain "header" or "file" not both')
    } else if ("header" %in% names(output_label_csv)) {
      added_label_by_header_target <- TRUE
      #output_label_csv$file <- stringr::str_to_lower(as.character(output_label_csv$header))
      data.table::setkey(output_label_csv, header)
    } else if ("file" %in% names(output_label_csv)) {
      added_label_by_header_target <- FALSE
      #output_label_csv$file <- stringr::str_to_lower(as.character(output_label_csv$file))
      data.table::setkey(output_label_csv, file)
    } else {
      stop('file in target_from_csv must contain one column named "header" or "file"')
    }
    
    # remove files without target label
    if (!added_label_by_header_target) {
      fasta.files <- fasta.files[basename(fasta.files) %in% output_label_csv$file]
      if (length(fasta.files) == 0) {
        stop("No overlap between files and 'file' column in target_from_csv")
      }
    }
    
    # if (!added_label_by_header_target) {
    #   fasta.file$Header <- rep(basename(fasta.files[file_index]), nrow(fasta.file))
    # }
    col_name <- ifelse(added_label_by_header_target, "header", "file")
    # header_vector <- fasta.file$Header
    vocabulary_label <- names(output_label_csv)
    vocabulary_label <- vocabulary_label[vocabulary_label != "header" & vocabulary_label != "file"]
    if (!is.null(target_split)) {
      check_header_names(target_split = target_split, vocabulary_label = vocabulary_label)
    }
    
    if (any(duplicated(output_label_csv$file))) {
      stop("csv file with label contains duplicate file names in 'file' column")
    }
  }
  
  # regular expression for chars outside vocabulary
  pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
  
  while (length(seq_vector) == 0) {
    
    # pre-load the first file
    fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                   shuffle_input = shuffle_input, proportion_entries = proportion_entries,
                                   reverse_complement = reverse_complement, fasta.files = fasta.files,
                                   vocabulary_label = vocabulary_label, filter_header = TRUE, target_from_csv = target_from_csv)
    
    if (concat) {
      if (use_coverage) {
        cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
      } 
      fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                               stringsAsFactors = FALSE)
    }
    
    # skip file that can't produce one sample
    if (!padding) {
      if (read_data) {
        seq_too_short <- all(nchar(as.character(fasta.file$Sequence)) < (maxlen/2))
      } else {
        seq_too_short <- all(nchar(as.character(fasta.file$Sequence)) < maxlen)
      }
      while((nrow(fasta.file) == 0) || seq_too_short) {
        file_index <- file_index + 1
        iter <- iter + 1
        if (file_index > length(fasta.files) || iter > max_iter) {
          stop("Can not extract enough samples, try reducing maxlen parameter")
        }
        
        fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                       shuffle_input = shuffle_input, proportion_entries = proportion_entries,
                                       reverse_complement = reverse_complement, fasta.files = fasta.files,
                                       vocabulary_label = vocabulary_label, filter_header = TRUE, target_from_csv = target_from_csv)
        
        if (concat) {
          if (use_coverage) {
            cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
          } 
          fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                                   stringsAsFactors = FALSE)
        }
        
        if (read_data) {
          seq_too_short <- all(nchar(as.character(fasta.file$Sequence)) < (maxlen/2))
        } else {
          seq_too_short <- all(nchar(as.character(fasta.file$Sequence)) < maxlen)
        }
      }
    } else {
      while(nrow(fasta.file) == 0) {
        file_index <- file_index + 1
        iter <- iter + 1
        if (file_index > length(fasta.files) || iter > max_iter) {
          stop("Can not extract enough samples, try reducing maxlen parameter")
        }
        fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                       shuffle_input = shuffle_input, proportion_entries = proportion_entries,
                                       reverse_complement = reverse_complement, fasta.files = fasta.files,
                                       vocabulary_label = vocabulary_label, filter_header = TRUE, target_from_csv = target_from_csv)
        
        if (concat) {
          if (use_coverage) {
            cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
          } 
          fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                                   stringsAsFactors = FALSE)
        }
      }
    }
    
    if (use_coverage) {
      cov_vector <- get_coverage(fasta.file)
    }
    
    # take random subset
    if (!is.null(proportion_per_seq)) {
      if (!read_data) {
        fasta_width <- nchar(fasta.file$Sequence)
        perc_length <- floor(fasta_width * proportion_per_seq)
        sample_range <- fasta_width - perc_length + 1
        start <- mapply(sample_range, FUN = sample, size = 1)
        stop <- start + perc_length - 1
        seq_vector <- mapply(fasta.file$Sequence, FUN = substr, start = start, stop = stop)
        if (use_quality_score) {
          quality_scores <- mapply(fasta.file$Quality, FUN = substr, start = start, stop = stop)
        }
      } else {
        if (use_quality_score) {
          quality_scores <- fasta.file$Quality
        }
      }
    } else {
      seq_vector <- fasta.file$Sequence
      if (use_quality_score) {
        quality_scores <- fasta.file$Quality
      }
    }
    
    seq_vector <- stringr::str_to_lower(seq_vector)
    seq_vector <- stringr::str_replace_all(string = seq_vector, pattern = pattern, amb_nuc_token)
    length_vector <- nchar(seq_vector)
    label_vector <- trimws(stringr::str_to_lower(fasta.file$Header))
    
    # label from csv
    if (!is.null(added_label_path)) {
      label_list <- list()
      # extra input from csv
      if (additional_labels) {
        if (length(added_label_path) != length(add_input_as_seq)) {
          stop("added_label_path and add_input_as_seq must have the same length")
        }
        added_label_list <- list()
        for (i in 1:length(added_label_path)) {
          added_label_list[[i]] <- input_from_csv(added_label_path[i])
        }
      }
      #added_label_by_header <- ifelse(added_label_list[[1]]$col_name == "header", TRUE, FALSE)
      added_label_by_header <- FALSE
    }
    
    # sequence vector collects strings until one batch can be created
    sequence_list <- vector("list")
    target_list <- vector("list")
    coverage_list <- vector("list")
    if (!use_quality_score) {
      quality_list <- NULL
    } else {
      quality_list <- vector("list")
    }
    
    if (!use_coverage) {
      coverage_list <- NULL
    } else {
      coverage_list <- vector("list")
    }
    
    if (!is.null(added_label_path)) {
      label_list <- vector("list")
    }
    
    if (!is.null(target_from_csv)) {
      output_label_list <- vector("list")
    }
    sequence_list_index <- 1
    
    # pad short sequences with zeros or discard
    short_seq_index <- which(length_vector < maxlen)
    if (padding) {
      for (i in short_seq_index) {
        seq_vector[i] <- paste0(paste(rep("0", maxlen - length_vector[i]), collapse = ""), seq_vector[i])
        if (use_quality_score) {
          quality_scores[i] <- paste0(paste(rep("!", maxlen - length_vector[i]), collapse = ""), quality_scores[i])
        }
        length_vector[i] <- maxlen
      }
    } else {
      if (length(short_seq_index) > 0) {
        seq_vector <- seq_vector[-short_seq_index]
        length_vector <- length_vector[-short_seq_index]
        label_vector <- label_vector[-short_seq_index]
        if (use_quality_score) {
          quality_scores <- quality_scores[-short_seq_index]
        }
        if (use_coverage) {
          cov_vector <- cov_vector[-short_seq_index]
        }
      }
    }
    
    nucSeq <- paste(seq_vector, collapse = "")
    
    if (use_quality_score) {
      quality_vector <- paste(quality_scores, collapse = "") %>% quality_to_probability()
    } else {
      quality_vector <- NULL
    }
    
    if (use_coverage) {
      cov_vector <- rep(cov_vector, times = nchar(seq_vector))
    } else {
      cov_vector <- NULL
    }
    
    if (length(seq_vector) == 0) {
      
      if(iter > max_iter) {
        stop('exceeded max_iter value, try reducing maxlen parameter')
        break
      }
      iter <- iter + 1
      
      file_index <- file_index + 1
      start_index <- 1
      
      if (file_index > length(fasta.files)) {
        if (shuffle_file_order) fasta.files <- sample(fasta.files, replace = FALSE)
        file_index <- 1
      }
    }
    
  }
  
  # vocabulary distribution
  nuc_dist_list <- vector("list")
  if (ambiguous_nuc == "empirical") {
    nuc_table <- table(stringr::str_split(nucSeq, ""))[vocabulary]
    nuc_dist <- vector("numeric")
    for (i in 1:length(vocabulary)) {
      nuc_dist[vocabulary[i]] <- nuc_table[vocabulary[i]]/sum(nuc_table)
    }
    nuc_dist[is.na(nuc_dist)] <- 0
    nuc_dist_list[[sequence_list_index]] <- nuc_dist
  } else {
    nuc_dist <- 0
  }
  
  startNewEntry <- cumsum(c(1, length_vector[-length(length_vector)]))
  if (!read_data) {
    start_indices <- get_start_ind(seq_vector = seq_vector, length_vector = length_vector, maxlen = maxlen, step = step,
                                   discard_amb_nuc = discard_amb_nuc, vocabulary = vocabulary)
  } else {
    start_indices <- startNewEntry
  }
  
  # limit samples per file
  if (!is.null(max_samples) && length(start_indices) > max_samples) {
    max_samples_subsample <- sample(1:(length(start_indices) - max_samples + 1), 1)
    start_indices <- start_indices[max_samples_subsample:(max_samples_subsample + max_samples - 1)]
  }
  
  nucSeq <- keras::texts_to_sequences(tokenizer_pred, nucSeq)[[1]] - 1
  
  # use subset of files
  if (!is.null(file_limit) && (file_limit < length(fasta.files))) {
    fasta.files <- fasta.files[1:file_limit]
    num_files <- length(fasta.files)
  }
  
  # log file
  if (!is.null(path_file_log)) {
    if (!endsWith(path_file_log, ".csv")) path_file_log <- paste0(path_file_log, ".csv")
    write.table(x = fasta.files[1], file = path_file_log, row.names = FALSE, col.names = FALSE)
  }
  
  if (verbose) message("Initializing ...")
  rngstate <- .GlobalEnv$.Random.seed
  
  function() {
    
    .GlobalEnv$.Random.seed <- rngstate
    on.exit(rngstate <<- .GlobalEnv$.Random.seed)
    iter <- 1
    # loop until enough samples collected
    while(num_samples < batch_size) {
      # loop through sub-sequences/files until sequence of suitable length is found
      while((start_index > length(start_indices)) | length(start_indices) == 0) {
        
        # go to next file
        if (sample_by_file_size) {
          file_index <<- sample(1:num_files, size = 1, prob = file_prob)
        } else {
          file_index <<- file_index + 1
        }
        start_index <<- 1
        
        if (file_index > length(fasta.files)) {
          if (shuffle_file_order) fasta.files <<- sample(fasta.files, replace = FALSE)
          file_index <<- 1
        }
        
        # skip empty files
        while(TRUE) {
          fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                         shuffle_input = shuffle_input, proportion_entries = proportion_entries,
                                         reverse_complement = reverse_complement, fasta.files = fasta.files,
                                         vocabulary_label = vocabulary_label, filter_header = TRUE, target_from_csv = target_from_csv)
          
          if (concat) {
            if (use_coverage) {
              cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
            } 
            fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                                     stringsAsFactors = FALSE)
          }
          
          if(iter > max_iter) {
            stop('exceeded max_iter value, try reducing maxlen or skip_amb_nuc parameter')
            break
          }
          iter <- iter + 1
          if (nrow(fasta.file) > 0) break
          file_index <<- file_index + 1
          if (file_index > length(fasta.files)) {
            if (shuffle_file_order) fasta.files <<- sample(fasta.files, replace = FALSE)
            file_index <<- 1
          }
        }
        
        if (use_coverage) {
          cov_vector <<- get_coverage(fasta.file)
        }
        
        # take random subset
        if (!is.null(proportion_per_seq)) {
          if (!read_data) {
            fasta_width <- nchar(fasta.file$Sequence)
            perc_length <- floor(fasta_width * proportion_per_seq)
            sample_range <- fasta_width - perc_length + 1
            start <- mapply(sample_range, FUN = sample, size = 1)
            stop <- start + perc_length - 1
            seq_vector <- mapply(fasta.file$Sequence, FUN = substr, start = start, stop = stop)
            if (use_quality_score) {
              quality_scores <- mapply(fasta.file$Quality, FUN = substr, start = start, stop = stop)
            }
          } else {
            if (use_quality_score) {
              quality_scores <- fasta.file$Quality
            }
          }
        } else {
          seq_vector <- fasta.file$Sequence
          if (use_quality_score) {
            quality_scores <- fasta.file$Quality
          }
        }
        
        seq_vector <- stringr::str_to_lower(seq_vector)
        seq_vector <- stringr::str_replace_all(string = seq_vector, pattern = pattern, amb_nuc_token)
        length_vector <- nchar(seq_vector)
        label_vector <- trimws(stringr::str_to_lower(fasta.file$Header))
        
        # log file
        if (!is.null(path_file_log)) {
          write.table(x = fasta.files[file_index], file = path_file_log, append = TRUE, col.names = FALSE, row.names = FALSE)
        }
        
        # pad short sequences with zeros or discard
        short_seq_index <<- which(length_vector < maxlen)
        if (padding) {
          for (i in short_seq_index) {
            seq_vector[i] <- paste0(paste(rep("0", maxlen - length_vector[i]), collapse = ""), seq_vector[i])
            if (use_quality_score) {
              quality_scores[i] <- paste0(paste(rep("!", maxlen - length_vector[i]), collapse = ""), quality_scores[i])
            }
            length_vector[i] <- maxlen
          }
        } else {
          if (length(short_seq_index) > 0) {
            seq_vector <- seq_vector[-short_seq_index]
            length_vector <- length_vector[-short_seq_index]
            label_vector <- label_vector[-short_seq_index]
            if (use_quality_score) {
              quality_scores <<- quality_scores[-short_seq_index]
            }
            if (use_coverage) {
              cov_vector <<- cov_vector[-short_seq_index]
            }
          }
        }
        
        label_vector <<- label_vector
        length_vector <<- length_vector
        
        # skip empty file
        if (length(seq_vector) == 0) {
          start_indices <<- NULL
          next
        }
        
        nucSeq <<- paste(seq_vector, collapse = "")
        
        if (use_quality_score) {
          quality_vector <<- paste(quality_scores, collapse = "") %>% quality_to_probability()
        } else {
          quality_vector <<- NULL
        }
        
        if (use_coverage) {
          cov_vector <<- rep(cov_vector, times = nchar(seq_vector))
        } else {
          cov_vector <<- NULL
        }
        
        # vocabulary distribution
        if (ambiguous_nuc == "empirical") {
          nuc_table <<- table(stringr::str_split(nucSeq, ""))[vocabulary]
          nuc_dist_temp <<- vector("numeric")
          for (i in 1:length(vocabulary)) {
            nuc_dist_temp[vocabulary[i]] <- nuc_table[vocabulary[i]]/sum(nuc_table)
          }
          nuc_dist_temp[is.na(nuc_dist)] <- 0
          nuc_dist <<- nuc_dist_temp
        }
        
        startNewEntry <<- cumsum(c(1, length_vector[-length(length_vector)]))
        if (!read_data) {
          start_indices <<- get_start_ind(seq_vector = seq_vector, length_vector = length_vector, maxlen = maxlen, step = step,
                                          discard_amb_nuc = discard_amb_nuc, vocabulary = vocabulary)
        } else {
          start_indices <<- startNewEntry
        }
        
        # limit samples per file
        if (!is.null(max_samples) && length(start_indices) > max_samples) {
          max_samples_subsample <- sample(1:(length(start_indices) - max_samples + 1), 1)
          start_indices <<- start_indices[max_samples_subsample:(max_samples_subsample + max_samples - 1)]
        }
        
        nucSeq <<- keras::texts_to_sequences(tokenizer_pred, nucSeq)[[1]] - 1
        
        if(iter > max_iter) {
          stop('exceeded max_iter value, try reducing maxlen parameter')
          break
        }
        iter <- iter + 1
      }
      
      # go as far as possible in sequence or stop when enough samples are collected
      remainingSamples <- batch_size - num_samples
      end_index <- min(length(start_indices), start_index + remainingSamples  - 1)
      
      subsetStartIndices <- start_indices[start_index:end_index]
      sequence_list[[sequence_list_index]] <- nucSeq[subsetStartIndices[1] : (subsetStartIndices[length(subsetStartIndices)] + maxlen - 1)]
      # collect targets
      if (is.null(target_from_csv)) {
        target_list[[sequence_list_index]] <- as.character(cut(subsetStartIndices, breaks = c(startNewEntry, length(nucSeq)),
                                                               labels = label_vector, include.lowest = TRUE, right = FALSE))
      }
      nuc_dist_list[[sequence_list_index]] <- nuc_dist
      
      if (!is.null(added_label_path)) {
        if (added_label_by_header) {
          label_list[[sequence_list_index]] <- as.character(cut(subsetStartIndices, breaks = c(startNewEntry, length(nucSeq)),
                                                                labels = header_vector, include.lowest = TRUE, right = FALSE))
        } else {
          label_list[[sequence_list_index]] <- basename(fasta.files[file_index])
        }
      }
      
      if (!is.null(target_from_csv)) {
        if (added_label_by_header_target) {
          output_label_list[[sequence_list_index]] <- as.character(cut(subsetStartIndices, breaks = c(startNewEntry, length(nucSeq)),
                                                                       labels = header_vector, include.lowest = TRUE, right = FALSE))
        } else {
          output_label_list[[sequence_list_index]] <- basename(fasta.files[file_index])
        }
      }
      
      if (use_quality_score) {
        quality_list[[sequence_list_index]] <- quality_vector[subsetStartIndices[1] : (subsetStartIndices[length(subsetStartIndices)] + maxlen - 1)]
      }
      
      if (use_coverage) {
        coverage_list[[sequence_list_index]] <- cov_vector[subsetStartIndices[1] : (subsetStartIndices[length(subsetStartIndices)] + maxlen - 1)]
      }
      
      start_index_list[[sequence_list_index]] <- subsetStartIndices
      sequence_list_index <<- sequence_list_index + 1
      num_new_samples <- end_index - start_index + 1
      num_samples <- num_samples + num_new_samples
      start_index <<- end_index + 1
    }
    
    # one hot encode strings collected in sequence_list and connect arrays
    array_x_list <- purrr::map(1:length(sequence_list), ~seq_encoding_label(sequence_list[[.x]], ambiguous_nuc = ambiguous_nuc, adjust_start_ind = TRUE,
                                                                            maxlen = maxlen, vocabulary = vocabulary, nuc_dist = nuc_dist_list[[.x]],
                                                                            start_ind =  start_index_list[[.x]], quality_vector = quality_list[[.x]],
                                                                            cov_vector = coverage_list[[.x]], return_int = return_int,
                                                                            use_coverage = use_coverage, max_cov = max_cov, n_gram = n_gram, 
                                                                            n_gram_stride = n_gram_stride)
    )
    
    # one hot encode targets
    if (is.null(target_from_csv)) {
      target_int <- unlist(keras::texts_to_sequences(tokenizer_target, unlist(target_list))) - 1
      y  <- keras::to_categorical(target_int, num_classes = length(vocabulary_label))
    }
    
    x <- array_x_list[[1]]
    
    if (length(array_x_list) > 1) {
      for (i in 2:length(array_x_list)) {
        x <- abind::abind(x, array_x_list[[i]], along = 1)
      }
    }
    
    if (additional_labels) {
      .datatable.aware = TRUE
      added_label_vector <- unlist(label_list) %>% stringr::str_to_lower()
      label_tensor_list <- list()
      for (i in 1:length(added_label_path)) {
        # added_label_by_header <- ifelse(added_label_list[[i]]$col_name == "header", TRUE, FALSE)
        label_tensor_list[[i]] <- csv_to_tensor(label_csv = added_label_list[[i]]$label_csv,
                                                added_label_vector = added_label_vector,
                                                added_label_by_header = added_label_by_header,
                                                batch_size = batch_size, start_index_list = start_index_list)
        if (add_input_as_seq[i]) {
          label_tensor_list[[i]] <- seq_encoding_label(as.vector(t(label_tensor_list[[i]])), nuc_dist = NULL, adjust_start_ind = TRUE,
                                                       maxlen = ncol(label_tensor_list[[i]]), vocabulary = vocabulary, ambiguous_nuc = ambiguous_nuc,
                                                       start_ind =  1 + ncol(label_tensor_list[[i]]) * (0:(nrow(label_tensor_list[[i]]) - 1)),
                                                       quality_vector = NULL)
        }
      }
    }
    
    if (!is.null(target_from_csv)) {
      .datatable.aware = TRUE
      output_label_vector <- unlist(output_label_list) # %>% stringr::str_to_lower()
      target_tensor <- matrix(0, ncol = ncol(output_label_csv) - 1, nrow = batch_size, byrow = TRUE)
      
      if (added_label_by_header_target) {
        header_unique <- unique(output_label_vector)
        for (i in header_unique) {
          output_label_from_csv <- output_label_csv[ .(i), -"header"]
          index__output_label_vector <- output_label_vector == i
          if (nrow(output_label_from_csv) > 0) {
            target_tensor[index__output_label_vector, ] <- matrix(as.matrix(output_label_from_csv[1, ]),
                                                                  nrow = sum(index_output_label_vector), ncol = ncol(target_tensor), byrow = TRUE)
          }
        }
      } else {
        row_index <- 1
        for (i in 1:length(output_label_vector)) {
          row_filter <- output_label_vector[i]
          output_label_from_csv <- output_label_csv[data.table(row_filter), -"file"]
          samples_per_file <- length(start_index_list[[i]])
          assign_rows <-  row_index:(row_index + samples_per_file - 1)
          
          if (nrow(na.omit(output_label_from_csv)) > 0) {
            target_tensor[assign_rows, ] <- matrix(as.matrix(output_label_from_csv[1, ]),
                                                   nrow = samples_per_file, ncol = ncol(target_tensor), byrow = TRUE)
          }
          row_index <- row_index + samples_per_file
        }
      }
      y <- target_tensor
    }
    
    # coerce y type to matrix
    if (dim(x)[1] == 1) {
      dim(y) <- c(1, length(vocabulary_label))
    }
    
    # empty sequence_list for next batch
    start_index_list <<- vector("list")
    sequence_list <<- vector("list")
    target_list <<- vector("list")
    nuc_dist_list <<- vector("list")
    if (use_quality_score) {
      quality_list <<- vector("list")
    }
    coverage_list <<- vector("list")
    sequence_list_index <<- 1
    num_samples <<- 0
    if (additional_labels) {
      if (length(x) == 2) {
        label_tensor_list[[length(label_tensor_list) + 1]] <- x[[1]]
        label_tensor_list[[length(label_tensor_list) + 1]] <- x[[2]]
        x <- label_tensor_list
      } else {
        label_tensor_list[[length(label_tensor_list) + 1]] <- x
        x <- label_tensor_list
      }
    }
    if (!is.null(target_split)) {
      colnames(y) <- vocabulary_label
      y <- slice_tensor(tensor = y, target_split = target_split)
    }
    
    if (!is.null(add_noise)) {
      noise_args <- c(add_noise, list(x = x))
      x <- do.call(add_noise_tensor, noise_args)
    }
    
    if (reverse_complement_encoding){
      x_1 <- x
      x_2 <- array(x_1[ , (dim(x)[2]):1, 4:1], dim = dim(x))
      x <- list(x_1, x_2)
    }
    
    return(list(X = x, Y = y))
  }
}

#' Data generator for fasta/fasta files
#'
#' @description Iterates over folder containing fasta/fastq files and produces encoding of predictor sequences
#' and target variables. Files in \code{path_corpus} should all belong to one class.  
#'
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_header_csv
#' @inheritParams train_model
#' @param num_targets Number of columns of target matrix.
#' @param ones_column Which column of target matrix contains ones.
#' @param read_data If `TRUE` the first element of output is a list of length 2, each containing one part of paired read. Maxlen should be 2*length of one read.
#' @param masked_lm If not `NULL`, input and target are equal except some parts of the input are masked or random.
#' Must be list with the following arguments: 
#' \itemize{
#' \item `mask_rate`: Rate of input to mask (rate of input to replace with mask token).
#' \item `random_rate`: Rate of input to set to random token.
#' \item `identity_rate`: Rate of input where sample weights are applied but input and output are identical. 
#' \item `include_sw`: Whether to include sample weigths.  
#' \item `block_len` (optional): Masked/random/itentity regions appear in blocks of size `block_len`.     
#' }
#' @export
generator_fasta_label_folder <- function(path_corpus,
                                         format = "fasta",
                                         batch_size = 256,
                                         maxlen = 250,
                                         max_iter = 10000,
                                         vocabulary = c("a", "c", "g", "t"),
                                         verbose = FALSE,
                                         shuffle_file_order = FALSE,
                                         step = 1,
                                         seed = 1234,
                                         shuffle_input = FALSE,
                                         file_limit = NULL,
                                         path_file_log = NULL,
                                         reverse_complement = TRUE,
                                         reverse_complement_encoding = FALSE,
                                         num_targets,
                                         ones_column,
                                         ambiguous_nuc = "zero",
                                         proportion_per_seq = NULL,
                                         read_data = FALSE,
                                         use_quality_score = FALSE,
                                         padding = TRUE,
                                         added_label_path = NULL,
                                         add_input_as_seq = NULL,
                                         skip_amb_nuc = NULL,
                                         max_samples = NULL,
                                         concat_seq = NULL,
                                         file_filter = NULL,
                                         use_coverage = NULL,
                                         proportion_entries = NULL,
                                         sample_by_file_size = FALSE,
                                         n_gram = NULL,
                                         n_gram_stride = 1,
                                         masked_lm = NULL,
                                         add_noise = NULL,
                                         return_int = FALSE) {
  
  #n_gram <- NULL
  if (is.null(use_coverage)) {
    use_coverage <- FALSE
    cov_vector <- NULL
    max_cov <- NULL
  } else {
    max_cov <- use_coverage
    use_coverage <- TRUE
  }
  if (!is.null(concat_seq) && (!all(stringr::str_split(concat_seq,"")[[1]] %in% vocabulary))) {
    stop("Characters of separating sequence should be in vocabulary")
  }
  if (reverse_complement_encoding) {
    test_len <- length(vocabulary) != 4
    if (test_len || all(sort(stringr::str_to_lower(vocabulary)) != c("a", "c", "g", "t"))) {
      stop("reverse_complement_encoding only implemented for A,C,G,T vocabulary yet")
    }
  }
  stopifnot(!(read_data & padding))
  stopifnot(ones_column <= num_targets)
  if (read_data & !is.null(skip_amb_nuc)) {
    stop("Using read data and skipping files at the same time not implemented yet")
  }
  additional_labels <- ifelse(is.null(added_label_path), FALSE, TRUE)
  
  # need to declare variables before nameless function() corpus for indexing
  #path_corpus <- path_corpus
  batch_size <- batch_size
  format <- format
  shuffle_file_order <- shuffle_file_order
  step <- step
  seed <- seed
  shuffle_input <- shuffle_input
  file_limit <- file_limit
  path_file_log <- path_file_log
  reverse_complement <- reverse_complement
  ambiguous_nuc <- ambiguous_nuc
  proportion_per_seq <- proportion_per_seq
  
  # # adjust maxlen for n_gram
  # if (!is.null(n_gram)) {
  #   stop("n-gram encoding not implemented yet for classification")
  #   maxlen <- maxlen + n_gram - 1
  # }
  
  discard_amb_nuc <- ifelse(ambiguous_nuc == "discard", TRUE, FALSE)
  vocabulary <- stringr::str_to_lower(vocabulary)
  start_index_list <- vector("list")
  file_index <- 1
  num_samples <- 0
  start_index <- 1
  iter <- 1
  concat <- !is.null(concat_seq)
  seq_vector <- NULL
  
  for (i in letters) {
    if (!(i %in% stringr::str_to_lower(vocabulary))) {
      amb_nuc_token <- i
      break
    }
  }
  tokenizer_pred <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), c(vocabulary, amb_nuc_token))
  
  # get fasta files
  fasta.files <- list_fasta_files(path_corpus = path_corpus,
                                  format = format,
                                  file_filter = file_filter)
  num_files <- length(fasta.files)
  
  if (sample_by_file_size) {
    shuffle_file_order <- FALSE
    file_prob <- file.info(fasta.files)$size/sum(file.info(fasta.files)$size)
  }
  
  set.seed(seed)
  if (shuffle_file_order) fasta.files <- sample(fasta.files, replace = FALSE)
  
  if (read_data) {
    contains_R1 <-  stringr::str_detect(fasta.files, "R1")
    fasta.files <- fasta.files[contains_R1]
  }
  
  # regular expression for chars outside vocabulary
  pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
  
  while (length(seq_vector) == 0) {
    
    fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                   shuffle_input = shuffle_input, proportion_entries = proportion_entries,
                                   reverse_complement = reverse_complement, fasta.files = fasta.files)
    
    if (concat) {
      if (use_coverage) {
        cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
      } 
      fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                               stringsAsFactors = FALSE)
    }
    
    # skip file that can't produce one sample
    if (!padding) {
      if (read_data) {
        seq_too_short <- all(nchar(as.character(fasta.file$Sequence)) < (maxlen/2))
      } else {
        seq_too_short <- all(nchar(as.character(fasta.file$Sequence)) < maxlen)
      }
      while((nrow(fasta.file) == 0) || seq_too_short) {
        file_index <- file_index + 1
        iter <- iter + 1
        if (file_index > length(fasta.files) || iter > max_iter) {
          stop("Can not extract enough samples, try reducing maxlen parameter")
        }
        fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                       shuffle_input = shuffle_input, proportion_entries = proportion_entries,
                                       reverse_complement = reverse_complement, fasta.files = fasta.files)
        
        if (concat) {
          if (use_coverage) {
            cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
          } 
          fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                                   stringsAsFactors = FALSE)
        }
        
        if (read_data) {
          seq_too_short <- all(nchar(as.character(fasta.file$Sequence)) < (maxlen/2))
        } else {
          seq_too_short <- all(nchar(as.character(fasta.file$Sequence)) < maxlen)
        }
      }
    } else {
      while(nrow(fasta.file) == 0) {
        file_index <- file_index + 1
        iter <- iter + 1
        if (file_index > length(fasta.files) || iter > max_iter) {
          stop("Can not extract enough samples, try reducing maxlen parameter")
        }
        fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                       shuffle_input = shuffle_input, proportion_entries = proportion_entries,
                                       reverse_complement = reverse_complement, fasta.files = fasta.files)
        
        if (concat) {
          if (use_coverage) {
            cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
          } 
          fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                                   stringsAsFactors = FALSE)
        }
      }
    }
    
    if (use_coverage) {
      cov_vector <- get_coverage(fasta.file)
    }
    
    # combine pairs to one string
    if (read_data) {
      second_read_path <- stringr::str_replace_all(fasta.files[file_index], pattern = "R1", replacement = "R2")
      if (format == "fasta") {
        fasta.file_2 <- microseq::readFasta(second_read_path)
      }
      if (format == "fastq") {
        fasta.file_2 <- microseq::readFastq(second_read_path)
      }
      df_1 <- as.data.frame(fasta.file)
      df_2 <- as.data.frame(fasta.file_2)
      fasta.file <- data.frame(Sequence = paste0(df_1$Sequence, df_2$Sequence), Quality = paste0(df_1$Quality, df_2$Quality))
    }
    
    # take random subset
    if (!is.null(proportion_per_seq)) {
      if (!read_data) {
        fasta_width <- nchar(fasta.file$Sequence)
        perc_length <- floor(fasta_width * proportion_per_seq)
        sample_range <- fasta_width - perc_length + 1
        start <- mapply(sample_range, FUN = sample, size = 1)
        stop <- start + perc_length - 1
        seq_vector <- mapply(fasta.file$Sequence, FUN = substr, start = start, stop = stop)
        if (use_quality_score) {
          quality_scores <- mapply(fasta.file$Quality, FUN = substr, start = start, stop = stop)
        }
      } else {
        # sample_index <- sample(nrow(fasta.file), ceiling(proportion_per_seq * nrow(fasta.file)))
        # fasta.file <- fasta.file[sample_index,]
        # seq_vector <- fasta.file$Sequence
        if (use_quality_score) {
          quality_scores <- fasta.file$Quality
        }
      }
    } else {
      seq_vector <- fasta.file$Sequence
      if (use_quality_score) {
        quality_scores <- fasta.file$Quality
      }
    }
    
    seq_vector <- stringr::str_to_lower(seq_vector)
    seq_vector <- stringr::str_replace_all(string = seq_vector, pattern = pattern, amb_nuc_token)
    length_vector <- nchar(seq_vector)
    
    # extra input from csv
    if (additional_labels) {
      label_list <- list()
      if (length(added_label_path) != length(add_input_as_seq)) {
        stop("added_label_path and add_input_as_seq must have the same length")
      }
      added_label_list <- list()
      for (i in 1:length(added_label_path)) {
        added_label_list[[i]] <- input_from_csv(added_label_path[i])
      }
      # added_label_by_header <- ifelse(added_label_list[[1]]$col_name == "header", TRUE, FALSE)
      added_label_by_header <- FALSE
    }
    
    # sequence vector collects strings until one batch can be created
    sequence_list <- vector("list")
    target_list <- vector("list")
    quality_list <- vector("list")
    coverage_list <- vector("list")
    
    if (!use_quality_score) {
      quality_list <- NULL
    }
    if (additional_labels) {
      label_list <- vector("list")
    }
    sequence_list_index <- 1
    
    # pad short sequences with zeros or discard
    short_seq_index <- which(length_vector < maxlen)
    if (padding) {
      for (i in short_seq_index) {
        seq_vector[i] <- paste0(paste(rep("0", maxlen - length_vector[i]), collapse = ""), seq_vector[i])
        if (use_quality_score) {
          quality_scores[i] <- paste0(paste(rep("!", maxlen - length_vector[i]), collapse = ""), quality_scores[i])
        }
        length_vector[i] <- maxlen
      }
    } else {
      if (length(short_seq_index) > 0) {
        seq_vector <- seq_vector[-short_seq_index]
        length_vector <- length_vector[-short_seq_index]
        if (use_quality_score) {
          quality_scores <- quality_scores[-short_seq_index]
        }
        if (use_coverage) {
          cov_vector <- cov_vector[-short_seq_index]
        }
        if (additional_labels) {
          header_vector <- header_vector[-short_seq_index]
        }
      }
    }
    
    if (length(seq_vector) == 0) {
      
      if(iter > max_iter) {
        stop('exceeded max_iter value, try reducing maxlen parameter')
        break
      }
      iter <- iter + 1
      
      file_index <- file_index + 1
      start_index <- 1
      
      if (file_index > length(fasta.files)) {
        if (shuffle_file_order) fasta.files <- sample(fasta.files, replace = FALSE)
        file_index <- 1
      }
    }
    
  }
  
  nucSeq <- paste(seq_vector, collapse = "")
  
  if (use_quality_score) {
    quality_vector <- paste(quality_scores, collapse = "") %>% quality_to_probability()
  } else {
    quality_vector <- NULL
  }
  
  if (use_coverage) {
    cov_vector <- rep(cov_vector, times = nchar(seq_vector))
  } else {
    cov_vector <- NULL
  }
  
  # vocabulary distribution
  nuc_dist_list <- vector("list")
  if (ambiguous_nuc == "empirical") {
    nuc_table <- table(stringr::str_split(nucSeq, ""))[vocabulary]
    nuc_dist <- vector("numeric")
    for (i in 1:length(vocabulary)) {
      nuc_dist[vocabulary[i]] <- nuc_table[vocabulary[i]]/sum(nuc_table)
    }
    nuc_dist[is.na(nuc_dist)] <- 0
    nuc_dist_list[[sequence_list_index]] <- nuc_dist
  } else {
    nuc_dist <- 0
  }
  
  startNewEntry <- cumsum(c(1, length_vector[-length(length_vector)]))
  if (!read_data) {
    start_indices <- get_start_ind(seq_vector = seq_vector, length_vector = length_vector, maxlen = maxlen, step = step,
                                   discard_amb_nuc = discard_amb_nuc, vocabulary = vocabulary)
  } else {
    start_indices <- startNewEntry
  }
  
  # limit samples per file
  if (!is.null(max_samples) && length(start_indices) > max_samples) {
    max_samples_subsample <- sample(1:(length(start_indices) - max_samples + 1), 1)
    start_indices <- start_indices[max_samples_subsample:(max_samples_subsample + max_samples - 1)]
  }
  
  nucSeq <- keras::texts_to_sequences(tokenizer_pred, nucSeq)[[1]] - 1
  
  # use subset of files
  if (!is.null(file_limit) && (file_limit < length(fasta.files))) {
    fasta.files <- fasta.files[1:file_limit]
    num_files <- length(fasta.files)
  }
  
  # log file
  if (!is.null(path_file_log)) {
    if (!endsWith(path_file_log, ".csv")) path_file_log <- paste0(path_file_log, ".csv")
    append <- file.exists(path_file_log)
    write.table(x = fasta.files[file_index], file = path_file_log, col.names = FALSE, row.names = FALSE, append = append)
  }
  
  if (verbose) message("Initializing ...")
  
  rngstate <- .GlobalEnv$.Random.seed
  
  function() {
    
    .GlobalEnv$.Random.seed <- rngstate
    on.exit(rngstate <<- .GlobalEnv$.Random.seed)
    
    iter <- 1
    
    # loop until enough samples collected
    while(num_samples < batch_size) {
      
      # loop through sub-sequences/files until sequence of suitable length is found
      while((start_index > length(start_indices)) | length(start_indices) == 0) {
        # go to next file
        if (sample_by_file_size) {
          file_index <<- sample(1:num_files, size = 1, prob = file_prob)
        } else {
          file_index <<- file_index + 1
        }
        start_index <<- 1
        
        if (file_index > length(fasta.files)) {
          if (shuffle_file_order) fasta.files <<- sample(fasta.files, replace = FALSE)
          file_index <<- 1
        }
        
        # skip empty files
        while(TRUE) {
          fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                         shuffle_input = shuffle_input, proportion_entries = proportion_entries,
                                         reverse_complement = reverse_complement, fasta.files = fasta.files)
          
          if (concat) {
            if (use_coverage) {
              cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
            } 
            fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                                     stringsAsFactors = FALSE)
          }
          
          if(iter > max_iter) {
            stop('exceeded max_iter value, try reducing maxlen or skip_amb_nuc parameter')
            break
          }
          iter <- iter + 1
          if (nrow(fasta.file) > 0) break
          file_index <<- file_index + 1
          if (file_index > length(fasta.files)) {
            if (shuffle_file_order) fasta.files <<- sample(fasta.files, replace = FALSE)
            file_index <<- 1
          }
        }
        
        if (use_coverage) {
          cov_vector <<- get_coverage(fasta.file)
        }
        
        # combine pairs to one string
        if (read_data) {
          second_read_path <- stringr::str_replace_all(fasta.files[file_index], pattern = "R1", replacement = "R2")
          if (format == "fasta") {
            fasta.file_2 <- microseq::readFasta(second_read_path)
          }
          if (format == "fastq") {
            fasta.file_2 <- microseq::readFastq(second_read_path )
          }
          df_1 <- as.data.frame(fasta.file)
          df_2 <- as.data.frame(fasta.file_2)
          fasta.file <- data.frame(Sequence = paste0(df_1$Sequence, df_2$Sequence), Quality = paste0(df_1$Quality, df_2$Quality))
        }
        
        # take random subset
        if (!is.null(proportion_per_seq)) {
          if (!read_data) {
            fasta_width <- nchar(fasta.file$Sequence)
            perc_length <- floor(fasta_width * proportion_per_seq)
            sample_range <- fasta_width - perc_length + 1
            start <- mapply(sample_range, FUN = sample, size = 1)
            stop <- start + perc_length - 1
            seq_vector <- mapply(fasta.file$Sequence, FUN = substr, start = start, stop = stop)
            if (use_quality_score) {
              quality_scores <- mapply(fasta.file$Quality, FUN = substr, start = start, stop = stop)
            }
          } else {
            # sample_index <- sample(nrow(fasta.file), ceiling(proportion_per_seq * nrow(fasta.file)))
            # fasta.file <- fasta.file[sample_index,]
            # seq_vector <- fasta.file$Sequence
            if (use_quality_score) {
              quality_scores <- fasta.file$Quality
            }
          }
        } else {
          seq_vector <- fasta.file$Sequence
          if (use_quality_score) {
            quality_scores <- fasta.file$Quality
          }
        }
        
        header_vector <<- fasta.file$Header
        seq_vector <- stringr::str_to_lower(seq_vector)
        seq_vector <- stringr::str_replace_all(string = seq_vector, pattern = pattern, amb_nuc_token)
        length_vector <- nchar(seq_vector)
        
        # log file
        if (!is.null(path_file_log)) {
          write.table(x = fasta.files[file_index], file = path_file_log, append = TRUE, col.names = FALSE, row.names = FALSE)
        }
        
        # pad short sequences with zeros or discard
        short_seq_index <<- which(length_vector < maxlen)
        if (padding) {
          for (i in short_seq_index) {
            seq_vector[i] <- paste0(paste(rep("0", maxlen - length_vector[i]), collapse = ""), seq_vector[i])
            if (use_quality_score) {
              quality_scores[i] <- paste0(paste(rep("!", maxlen - length_vector[i]), collapse = ""), quality_scores[i])
            }
            length_vector[i] <- maxlen
          }
        } else {
          if (length(short_seq_index) > 0) {
            seq_vector <- seq_vector[-short_seq_index]
            length_vector <- length_vector[-short_seq_index]
            if (use_quality_score) {
              quality_scores <- quality_scores[-short_seq_index]
            }
            if (use_coverage) {
              cov_vector <<- cov_vector[-short_seq_index]
            }
            if (additional_labels) {
              header_vector <- header_vector[-short_seq_index]
              header_vector <<- header_vector
            }
          }
        }
        
        # skip empty file
        if (length(seq_vector) == 0) {
          start_indices <<- NULL
          next
        }
        
        nucSeq <<- paste(seq_vector, collapse = "")
        if (use_quality_score) {
          quality_vector <<- paste(quality_scores, collapse = "") %>% quality_to_probability()
        } else {
          quality_vector <<- NULL
        }
        if (use_coverage) {
          cov_vector <<- rep(cov_vector, times = nchar(seq_vector))
        } else {
          cov_vector <<- NULL
        }
        
        # vocabulary distribution
        if (ambiguous_nuc == "empirical") {
          nuc_table <<- table(stringr::str_split(nucSeq, ""))[vocabulary]
          nuc_dist_temp <<- vector("numeric")
          for (i in 1:length(vocabulary)) {
            nuc_dist_temp[vocabulary[i]] <- nuc_table[vocabulary[i]]/sum(nuc_table)
          }
          nuc_dist_temp[is.na(nuc_dist)] <- 0
          nuc_dist <<- nuc_dist_temp
        }
        
        startNewEntry <<- cumsum(c(1, length_vector[-length(length_vector)]))
        if (!read_data) {
          start_indices <<- get_start_ind(seq_vector = seq_vector, length_vector = length_vector, maxlen = maxlen, step = step,
                                          discard_amb_nuc = discard_amb_nuc, vocabulary = vocabulary)
        } else {
          start_indices <<- startNewEntry
        }
        
        # limit samples per file
        if (!is.null(max_samples) && length(start_indices) > max_samples) {
          max_samples_subsample <- sample(1:(length(start_indices) - max_samples + 1), 1)
          start_indices <<- start_indices[max_samples_subsample:(max_samples_subsample + max_samples - 1)]
        }
        
        nucSeq <<- keras::texts_to_sequences(tokenizer_pred, nucSeq)[[1]] - 1
        
        if(iter > max_iter) {
          stop('exceeded max_iter value, try reducing maxlen parameter')
          break
        }
        iter <- iter + 1
      }
      
      # go as far as possible in sequence or stop when enough samples are collected
      remainingSamples <- batch_size - num_samples
      end_index <- min(length(start_indices), start_index + remainingSamples  - 1)
      
      subsetStartIndices <- start_indices[start_index:end_index]
      sequence_list[[sequence_list_index]] <- nucSeq[subsetStartIndices[1] : (subsetStartIndices[length(subsetStartIndices)] + maxlen - 1)]
      if (use_quality_score) {
        quality_list[[sequence_list_index]] <- quality_vector[subsetStartIndices[1] : (subsetStartIndices[length(subsetStartIndices)] + maxlen - 1)]
      }
      if (use_coverage) {
        coverage_list[[sequence_list_index]] <- cov_vector[subsetStartIndices[1] : (subsetStartIndices[length(subsetStartIndices)] + maxlen - 1)]
      }
      nuc_dist_list[[sequence_list_index]] <- nuc_dist
      start_index_list[[sequence_list_index]] <- subsetStartIndices
      if (additional_labels) {
        if (added_label_by_header) {
          label_list[[sequence_list_index]] <- as.character(cut(subsetStartIndices, breaks = c(startNewEntry, length(nucSeq)),
                                                                labels = header_vector, include.lowest = TRUE, right = FALSE))
        } else {
          label_list[[sequence_list_index]] <- basename(fasta.files[file_index])
        }
      }
      sequence_list_index <<- sequence_list_index + 1
      
      num_new_samples <- end_index - start_index + 1
      num_samples <- num_samples + num_new_samples
      start_index <<- end_index + 1
    }
    
    # one hot encode strings collected in sequence_list and connect arrays
    if (is.null(masked_lm)) {
      
      array_x_list <- purrr::map(1:length(sequence_list), ~seq_encoding_label(sequence = sequence_list[[.x]], ambiguous_nuc = ambiguous_nuc,
                                                                              maxlen = maxlen, vocabulary = vocabulary, nuc_dist = nuc_dist_list[[.x]],
                                                                              start_ind =  start_index_list[[.x]], return_int = return_int,
                                                                              quality_vector = quality_list[[.x]], cov_vector = coverage_list[[.x]],
                                                                              max_cov = max_cov, use_coverage = use_coverage, n_gram = n_gram,
                                                                              n_gram_stride = n_gram_stride, adjust_start_ind = TRUE)
      )
      
      x <- array_x_list[[1]]
      
      if (length(array_x_list) > 1) {
        for (i in 2:length(array_x_list)) {
          x <- abind::abind(x, array_x_list[[i]], along = 1)
        }
      }
      
      # one hot encode targets
      y  <- matrix(0, ncol = num_targets, nrow = dim(x)[1])
      y[ , ones_column] <- 1
      
      # coerce y type to matrix
      if (dim(x)[1] == 1) {
        dim(y) <- c(1, num_targets)
      }
      
    } else {
      
      if (is.null(masked_lm$include_sw)) {
        include_sw <- FALSE
      } else {
        include_sw <- masked_lm$include_sw
      }
      array_lists <- purrr::map(1:length(sequence_list), ~seq_encoding_label(sequence = sequence_list[[.x]], ambiguous_nuc = ambiguous_nuc,
                                                                             maxlen = maxlen, vocabulary = vocabulary, nuc_dist = nuc_dist_list[[.x]],
                                                                             start_ind =  start_index_list[[.x]], masked_lm = masked_lm,
                                                                             quality_vector = quality_list[[.x]], cov_vector = coverage_list[[.x]],
                                                                             max_cov = max_cov, use_coverage = use_coverage, n_gram = n_gram,
                                                                             n_gram_stride = n_gram_stride, adjust_start_ind = TRUE, return_int = return_int)
      )
      
      x_y_sw <- reorder_masked_lm_lists(array_lists, include_sw = include_sw)
      x <- x_y_sw$x
      y <- x_y_sw$y
      if (include_sw) {
        sample_weight <- x_y_sw$sw
      } else {
        sample_weight <- NULL
      }
      rm(x_y_sw)
    }
    
    if (reverse_complement_encoding){
      x_1 <- x
      x_2 <- array(x_1[ , (dim(x)[2]):1, 4:1], dim = dim(x))
      x <- list(x_1, x_2)
    }
    
    if (read_data) {
      input_dim <- dim(x)/c(1,2,1)
      x_1 <- array(x[ , 1:(maxlen/2), ], dim = input_dim)
      x_2 <- array(x[ , ((maxlen/2) + 1) : maxlen, ], dim = input_dim)
      x <- list(x_1, x_2)
    }
    
    if (additional_labels) {
      .datatable.aware = TRUE
      added_label_vector <- unlist(label_list) %>% stringr::str_to_lower()
      label_tensor_list <- list()
      for (i in 1:length(added_label_path)) {
        # added_label_by_header <- ifelse(added_label_list[[i]]$col_name == "header", TRUE, FALSE)
        label_tensor_list[[i]] <- csv_to_tensor(label_csv = added_label_list[[i]]$label_csv,
                                                added_label_vector = added_label_vector,
                                                added_label_by_header = added_label_by_header,
                                                batch_size = batch_size, start_index_list = start_index_list)
        if (add_input_as_seq[i]) {
          label_tensor_list[[i]] <- seq_encoding_label(as.vector(t(label_tensor_list[[i]])), nuc_dist = NULL, adjust_start_ind = TRUE,
                                                       maxlen = ncol(label_tensor_list[[i]]), vocabulary = vocabulary, ambiguous_nuc = ambiguous_nuc,
                                                       start_ind =  1 + ncol(label_tensor_list[[i]]) * (0:(nrow(label_tensor_list[[i]]) - 1)), 
                                                       quality_vector = NULL)
        }
      }
    }
    
    # empty lists for next batch
    start_index_list <<- vector("list")
    sequence_list <<- vector("list")
    target_list <<- vector("list")
    if (use_quality_score) {
      quality_list <<- vector("list")
    }
    nuc_dist_list <<- vector("list")
    coverage_list <<- vector("list")
    sequence_list_index <<- 1
    num_samples <<- 0
    if (reverse_complement_encoding) {
      return(list(X = x, Y = y))
    }
    
    if (additional_labels) {
      if (length(x) == 2) {
        label_tensor_list[[length(label_tensor_list) + 1]] <- x[[1]]
        label_tensor_list[[length(label_tensor_list) + 1]] <- x[[2]]
        x <- label_tensor_list
      } else {
        label_tensor_list[[length(label_tensor_list) + 1]] <- x
        x <- label_tensor_list
      }
    }
    
    if (!is.null(add_noise)) {
      noise_args <- c(add_noise, list(x = x))
      x <- do.call(add_noise_tensor, noise_args)
    }
    
    if (!is.null(masked_lm) && include_sw) return(list(x, y, sample_weight))
    return(list(X = x, Y = y))
    
  }
}

#' Initializes generators defined by \code{generator_fasta_label_folder} function
#'
#' Initializes generators defined by \code{\link{generator_fasta_label_folder}} function. Targets get encoded in order of directories.
#' Number of classes is given by length of \code{directories}.
#'
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_header_csv
#' @inheritParams generator_fasta_label_folder
#' @inheritParams train_model
#' @param directories Vector of paths to folder containing fasta files. Files in one folder should belong to one class.
#' @param val Logical, call initialized generator "genY" or "genValY" where Y is an integer between 1 and length of directories.
#' @param target_middle Split input sequence into two sequences while removing nucleotide in middle. If input is x_1,..., x_(n+1), input gets split into
#' input_1 = x_1,..., x_m and input_2 = x_(n+1),..., x_(m+2) where m = ceiling((n+1)/2) and n = maxlen. Note that x_(m+1) is not used.
#' @param read_data If true the first element of output is a list of length 2, each containing one part of paired read.
#' @export
generator_initialize <- function(directories,
                                 format = "fasta",
                                 batch_size = 256,
                                 maxlen = 250,
                                 max_iter = 10000,
                                 vocabulary = c("a", "c", "g", "t"),
                                 verbose = FALSE,
                                 shuffle_file_order = FALSE,
                                 step = 1,
                                 seed = 1234,
                                 shuffle_input = FALSE,
                                 file_limit = NULL,
                                 path_file_log = NULL,
                                 reverse_complement = FALSE,
                                 reverse_complement_encoding = FALSE,
                                 val = FALSE,
                                 ambiguous_nuc = "zero",
                                 proportion_per_seq = NULL,
                                 target_middle = FALSE,
                                 read_data = FALSE,
                                 use_quality_score = FALSE,
                                 padding = TRUE,
                                 added_label_path = NULL,
                                 add_input_as_seq = NULL,
                                 skip_amb_nuc = NULL,
                                 max_samples = NULL,
                                 file_filter = NULL,
                                 concat_seq = NULL,
                                 use_coverage = NULL,
                                 set_learning = NULL,
                                 proportion_entries = NULL,
                                 sample_by_file_size = FALSE,
                                 n_gram = NULL,
                                 n_gram_stride = 1,
                                 add_noise = NULL,
                                 return_int = FALSE) {
  
  # adjust batch_size
  if (is.null(set_learning) && (length(batch_size) == 1) && (batch_size %% length(directories) != 0)) {
    batch_size <- ceiling(batch_size/length(directories)) * length(directories)
    if (!val) {
      message(paste("Batch size needs to be multiple of number of targets. Setting batch_size to", batch_size))
    }
  }
  
  num_targets <- length(directories)
  
  if (!is.null(set_learning)) {
    reshape_mode <- set_learning$reshape_mode
    samples_per_target <- set_learning$samples_per_target
    buffer_len <- set_learning$buffer_len
    maxlen <- set_learning$maxlen
    concat_maxlen <- NULL
    
    if (reshape_mode == "concat") {
      
      if (sum(batch_size) %% length(directories) != 0) {
        stop_text <- paste("batch_size is", batch_size, "but needs to be multiple of number of classes (",
                           length(directories), ") for set learning with 'concat'")
        stop(stop_text)
      }
      buffer_len <- ifelse(is.null(set_learning$buffer_len), 0, set_learning$buffer_len)
      concat_maxlen <- (maxlen * samples_per_target) + (buffer_len * (samples_per_target - 1))
      if (any(c("z", "Z") %in% vocabulary) & !is.null(set_learning$buffer_len)) {
        stop("'Z' is used as token for separating sequences and can not be in vocabulary.")
      }
      if (!is.null(set_learning$buffer_len)) {
        vocabulary <- c(vocabulary, "Z")
      }
    }
    
    if (any(batch_size[1] != batch_size)) {
      stop("Set learning only implemented for uniform batch_size for all classes.")
    }
    new_batch_size <- batch_size
    batch_size <- samples_per_target * batch_size
  }
  
  
  if (length(batch_size) == 1) {
    batch_size <- rep(batch_size/num_targets, num_targets)
  }
  
  argg <- c(as.list(environment()))
  # variables with just one entry
  argg["directories"] <- NULL
  argg["file_filter"] <- NULL
  argg["val"] <- NULL
  argg["vocabulary"] <- NULL
  argg["num_targets"] <- NULL
  argg["verbose"] <- NULL
  argg["maxlen"] <- NULL
  argg["max_iter"] <- NULL
  argg["read_data"] <- NULL
  argg["use_quality_score"] <- NULL
  argg["added_label_path"] <- NULL
  argg["add_input_as_seq"] <- NULL
  argg["skip_amb_nuc"] <- NULL
  argg["concat_seq"] <- NULL
  argg["reverse_complement_encoding"] <- NULL
  argg["use_coverage"] <- NULL
  argg["set_learning"] <- NULL
  argg["proportion_entries"] <- NULL
  argg["sample_by_file_size"] <- NULL
  argg["n_gram"] <- NULL
  argg["n_gram_stride"] <- NULL
  argg[["add_noise"]] <- NULL
  argg[["return_int"]] <- NULL
  
  for (i in 1:length(argg)) {
    if (length(argg[[i]]) == 1) {
      assign(names(argg)[i], rep(argg[[i]], num_targets))
    }
    if ((length(argg[[i]]) != 1) & (length(argg[[i]]) != num_targets) & !(is.null(argg[[i]]))) {
      stop_message <- paste("Incorrect argument length,", names(argg[i]), "argument vector must have length 1 or", num_targets)
      stop(stop_message)
    }
  }
  
  if (!val) {
    # create generator for every folder
    for (i in 1:length(directories)) {
      numberedGen <- paste0("gen", as.character(i))
      genAsText <- paste(numberedGen, "<<- generator_fasta_label_folder(path_corpus = directories[[i]],
                                         format = format[i],
                                         batch_size = batch_size[i],
                                         maxlen = maxlen,
                                         max_iter = max_iter,
                                         vocabulary = vocabulary,
                                         verbose = verbose,
                                         shuffle_file_order = shuffle_file_order[i],
                                         step = step[i],
                                         seed = seed[i],
                                         shuffle_input = shuffle_input[i],
                                         file_limit = file_limit[i],
                                         path_file_log = path_file_log[i],
                                         reverse_complement = reverse_complement[i],
                                         reverse_complement_encoding = reverse_complement_encoding,
                                         num_targets = num_targets,
                                         ones_column = i,
                                         ambiguous_nuc = ambiguous_nuc[i],
                                         proportion_per_seq = proportion_per_seq[i],
                                         read_data = read_data,
                                         use_quality_score = use_quality_score,
                                         padding = padding[i],
                                         file_filter = file_filter,
                                         added_label_path = added_label_path,
                                         add_input_as_seq = add_input_as_seq,
                                         skip_amb_nuc = skip_amb_nuc,
                                         max_samples = max_samples[i],
                                         concat_seq = concat_seq,
                                         use_coverage = use_coverage,
                                         proportion_entries = proportion_entries,
                                         sample_by_file_size = sample_by_file_size,
                                         n_gram = n_gram,
                                         masked_lm = NULL,
                                         n_gram_stride = n_gram_stride,
                                         add_noise = add_noise,
                                         return_int = return_int
    )"
      )
      eval(parse(text = genAsText))
    }
  } else {
    # create generator for every folder
    for (i in 1:length(directories)) {
      # different names for validation generators
      numberedGenVal <- paste0("genVal", as.character(i))
      genAsTextVal <- paste(numberedGenVal, "<<- generator_fasta_label_folder(path_corpus = directories[[i]],
                                         format = format[i],
                                         batch_size = batch_size[i],
                                         maxlen = maxlen,
                                         max_iter = max_iter,
                                         vocabulary = vocabulary,
                                         verbose = verbose,
                                         shuffle_file_order = shuffle_file_order[i],
                                         step = step[i],
                                         seed = seed[i],
                                         shuffle_input = shuffle_input[i],
                                         file_limit = file_limit[i],
                                         path_file_log = path_file_log[i],
                                         reverse_complement = reverse_complement[i],
                                         reverse_complement_encoding = reverse_complement_encoding,
                                         num_targets = num_targets,
                                         ones_column = i,
                                         ambiguous_nuc = ambiguous_nuc[i],
                                         proportion_per_seq = proportion_per_seq[i],
                                         read_data = read_data,
                                         use_quality_score = use_quality_score,
                                         padding = padding[i],
                                         file_filter = file_filter,
                                         added_label_path = added_label_path,
                                         add_input_as_seq = add_input_as_seq,
                                         skip_amb_nuc = skip_amb_nuc,
                                         max_samples = max_samples[i],
                                         concat_seq = concat_seq,
                                         use_coverage = use_coverage,
                                         proportion_entries = proportion_entries,
                                         sample_by_file_size = sample_by_file_size,
                                         masked_lm = NULL,
                                         n_gram = n_gram,
                                         n_gram_stride = n_gram_stride,
                                         add_noise = add_noise,
                                         return_int = return_int
    )"
      )
      eval(parse(text = genAsTextVal))
    }
  }
}

#' Generator wrapper
#'
#' Combines generators created by \code{\link{generator_initialize}} into a single generator.
#'
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_folder
#' @inheritParams train_model
#' @inheritParams reshape_tensor
#' @param val Train or validation generator.
#' @param path Path to input files.
#' @param new_batch_size Only applied if \code{set_learning} is not NULL
#' @param voc_len Length of vocabulary.
#' @export
generator_fasta_label_folder_wrapper <- function(val, new_batch_size = NULL,
                                                 batch_size = NULL,
                                                 path = NULL, voc_len = NULL, 
                                                 maxlen = NULL,
                                                 set_learning = NULL) {
  
  if (is.null(set_learning)) {
    samples_per_target <- NULL
    new_batch_size <- NULL
    reshape_mode <- NULL
    buffer_len <- NULL
  } else {
    reshape_mode <- set_learning$reshape_mode
    samples_per_target <- set_learning$samples_per_target
    buffer_len <- set_learning$buffer_len
    maxlen <- set_learning$maxlen
    concat_maxlen <- NULL
    
    if (reshape_mode == "concat") {
      
      if (sum(batch_size) %% length(path) != 0) {
        stop_text <- paste("batch_size is", batch_size, "but needs to be multiple of number of classes (",
                           length(path), ") for set learning with 'concat'")
        stop(stop_text)
      }
      buffer_len <- ifelse(is.null(set_learning$buffer_len), 0, set_learning$buffer_len)
      concat_maxlen <- (maxlen * samples_per_target) + (buffer_len * (samples_per_target - 1))
      if (!is.null(set_learning$buffer_len)) {
        vocabulary <- c(vocabulary, "Z")
      }
      
    }
    
    if (any(batch_size[1] != batch_size)) {
      stop("Set learning only implemented for uniform batch_size for all classes.")
    }
    new_batch_size <- batch_size
    batch_size <- samples_per_target * batch_size
  }
  
  if (!val) {
    function() {
      directories <- path
      # combine generators to create one batch
      subBatchTrain <<- eval(parse(text = paste0("gen", as.character("1"), "()")))
      if (is.list(subBatchTrain[[1]])) {
        num_inputs <- length(subBatchTrain[[1]])
      } else {
        num_inputs <- 1
      }
      xTrain <- subBatchTrain[[1]]
      yTrain <- subBatchTrain[[2]]
      if (num_inputs > 1) {
        x_train_list <- list()
        for (i in 1:num_inputs) {
          x_train_list[[i]] <- xTrain[[i]]
        }
      }
      
      if (length(directories) > 1) {
        for (i in 2:length(directories)) {
          subBatchTrain <<- eval(parse(text = paste0("gen", as.character(i), "()")))
          yTrain <- rbind(yTrain, subBatchTrain[[2]])
          
          if (num_inputs == 1) {
            xTrain <- abind::abind(xTrain, subBatchTrain[[1]], along = 1)
          } else {
            for (j in 1:num_inputs) {
              x_train_list[[j]] <- abind::abind(x_train_list[[j]], subBatchTrain[[1]][[j]], along = 1)
            }
          }
        }
      }
      if (num_inputs > 1) {
        xTrain <- x_train_list
      }
      
      if (!is.null(samples_per_target)) {
        l <- reshape_tensor(x = xTrain, y = yTrain,
                            new_batch_size = new_batch_size,
                            samples_per_target = samples_per_target,
                            buffer_len = buffer_len,
                            reshape_mode = reshape_mode)
        return(l)
      } else {
        return(list(X = xTrain, Y = yTrain))
      }
      
    }
  } else {
    function() {
      directories <- path
      # combine generators to create one batch
      subBatchVal <<- eval(parse(text = paste0("genVal", as.character("1"), "()")))
      if (is.list(subBatchVal[[1]])) {
        num_inputs <- length(subBatchVal[[1]])
      } else {
        num_inputs <- 1
      }
      xVal <- subBatchVal[[1]]
      yVal <- subBatchVal[[2]]
      if (num_inputs > 1) {
        x_val_list <- list()
        for (i in 1:num_inputs) {
          x_val_list[[i]] <- xVal[[i]]
        }
      }
      
      if (length(directories) > 1) {
        for (i in 2:length(directories)) {
          subBatchVal <<- eval(parse(text = paste0("genVal", as.character(i), "()")))
          yVal <- rbind(yVal, subBatchVal[[2]])
          
          if (num_inputs == 1) {
            xVal <- abind::abind(xVal, subBatchVal[[1]], along = 1)
          } else {
            for (j in 1:num_inputs) {
              x_val_list[[j]] <- abind::abind(x_val_list[[j]], subBatchVal[[1]][[j]], along = 1)
            }
          }
        }
      }
      if (num_inputs > 1) {
        xVal <- x_val_list
      }
      if (!is.null(samples_per_target)) {
        l <- reshape_tensor(x = xVal, y = yVal,
                            new_batch_size = new_batch_size,
                            samples_per_target = samples_per_target,
                            #batch_size = batch_size,
                            #path = path,
                            #voc_len = voc_len,
                            #maxlen = maxlen,
                            #concat_maxlen = concat_maxlen,
                            buffer_len = buffer_len,
                            reshape_mode = reshape_mode)
        return(l)
      } else {
        return(list(X = xVal, Y = yVal))
      }
    }
  }
}

#' Random data generator
#' 
#' Creates a random input/target list once and repeatedly returns list. 
#'
#' @inheritParams generator_fasta_lm
#' @param model A keras model.
#' @param sparse_loss Whether to adjust target to sparse loss.
#' @examples 
#' model <- create_model_lstm_cnn(
#'   maxlen = 10,
#'   layer_lstm = c(4),
#'   layer_dense = c(5))
#' gen <- generator_dummy(model, 12)
#' z <- gen()
#' x <- z[[1]]
#' y <- z[[2]]
#' dim(x)
#' dim(y)
#' @export
generator_dummy <- function(model, batch_size) {
  
  # sparse loss
  if (!is.null(model$loss) && 
      stringr::str_detect(stringr::str_to_lower(model$loss), "sparse")) {
    sparse_loss <- TRUE
  } else {
    sparse_loss <- FALSE
  }

  # stateful model 
  if (!is.null(model$input_shape[[1]][[1]])) {
    batch_size <- model$input_shape[[1]]
  }
  
  num_input_layers <- ifelse(is.list(model$input), length(model$inputs), 1)
  x <- list()
  if (num_input_layers == 1) {
    input_dim <- batch_size
    for (j in 2:length(model$input_shape)) {
      input_dim  <- c(input_dim, model$input_shape[[j]])
    }
    x <- array(runif(n = prod(input_dim)), dim = input_dim)
  } else {
    
    for (i in 1:num_input_layers) {
      input_dim <- batch_size
      input_list <- model$input_shape[[i]]
      for (j in 2:length(input_list)) {
        input_dim  <- c(input_dim, input_list[[j]])
      }
      x[[i]] <- array(runif(n = prod(input_dim)), dim = input_dim)
    }
  }
  
  num_output_layers <- ifelse(is.list(model$output), length(model$outputs), 1)
  y <- list()
  if (num_output_layers == 1) {
    output_dim <- batch_size
    for (j in 2:length(model$output_shape)) {
      output_dim  <- c(output_dim, model$output_shape[[j]])
    }
    if (sparse_loss) output_dim <- output_dim[-length(output_dim)]
    y <- array(runif(n = prod(output_dim)), dim = output_dim)
  } else {
    for (i in 1:num_output_layers) {
      output_dim <- batch_size
      output_list <- model$output_shape[[i]]
      for (j in 2:length(output_list)) {
        output_dim  <- c(output_dim, output_list[[j]])
      }
      if (sparse_loss) output_dim <- output_dim[-length(output_dim)]
      y[[i]] <- array(runif(n = prod(output_dim)), dim = output_dim)
    }
  }
  
  function() {
    return(list(x,y))
  }
}

#' Rds data generator
#' 
#' Creates training batches from rds files. Rds files must contain a  
#' list of length 2 (input/target) or of length 1 (for language model).
#' If \code{target_len} is not NULL will take the last \code{target_len} entries of 
#' the first list element as targets and the rest as input.    
#' 
#' @inheritParams generator_fasta_label_header_csv
#' @param rds_folder Path to input files.
#' @param target_len Number of target nucleotides for language model.
#' @param delete_used_files Whether to delete file once used. Only applies for rds files. 
#' @examples 
#' # create 3 rds files
#' rds_folder <- tempfile()
#' dir.create(rds_folder)
#' batch_size <- 7
#' maxlen <- 11
#' voc_len <- 4
#' for (i in 1:3) {
#'   x <- sample(0:(voc_len-1), maxlen*batch_size, replace = TRUE)
#'   x <- keras::to_categorical(x, num_classes = voc_len)
#'   x <- array(x, dim = c(batch_size, maxlen, voc_len))
#'   y <- sample(0:2, batch_size ,replace = TRUE)
#'   y <- keras::to_categorical(y, num_classes = 3)
#'   xy_list <- list(x, y)
#'   file_name <- paste0(rds_folder, "/file_", i, ".rds")
#'   saveRDS(xy_list, file_name) 
#' }
#' 
#' # create generator
#' gen <- generator_rds(rds_folder, batch_size = 2)
#' z <- gen()
#' x <- z[[1]]
#' y <- z[[2]]
#' x[1, , ]
#' y[1, ]
#' @export
generator_rds <- function(rds_folder, batch_size, path_file_log = NULL,
                          max_samples = NULL,
                          proportion_per_seq = NULL,
                          target_len = NULL,
                          seed = NULL, delete_used_files = FALSE,
                          reverse_complement = FALSE,
                          sample_by_file_size = FALSE,
                          n_gram = NULL, n_gram_stride = 1,
                          reverse_complement_encoding = FALSE,
                          add_noise = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  is_lm <- !is.null(target_len)
  
  if (!is.null(n_gram) & is_lm && (target_len < n_gram)) {
    stop("target_len needs to be at least as big as n_gram.")
  }
  
  rds_files <- list_fasta_files(rds_folder, format = "rds", file_filter = NULL)
  num_files <- length(rds_files)
  
  read_success <- FALSE
  while (!read_success) {
    tryCatch(
      expr = {
        if (sample_by_file_size) {
          file_prob <- file.info(rds_files)$size/sum(file.info(rds_files)$size)
          file_index <- sample(1:num_files, size = 1, prob = file_prob)
        } else {
          file_prob <- NULL
          rds_files <- sample(rds_files)
          file_index <- 1
        }
        rds_file <- readRDS(rds_files[file_index])
        read_success <- TRUE
      },
      error = function(e){ 
        
      }
    )
  }
  
  if (is.list(rds_file)) {
    x_complete <- rds_file[[1]]
  } else {
    x_complete <- rds_file
  }
  
  if (!is_lm) y_complete <- rds_file[[2]]
  # TODO: adjust for different input fomrat (input mix of 3D and 1D etc.)
  multi_input <- ifelse(is.list(x_complete), TRUE, FALSE)
  multi_output <- ifelse(length(rds_file) > 1 && is.list(rds_file[[2]]), TRUE, FALSE)
  
  if (multi_input) {
    x_dim_list <- list()
    size_splits_in <- list()
    for (i in 1:length(x_complete)) {
      x_dim_list[[i]] <- dim(x_complete[[i]])
      size_splits_in[[i]] <- x_dim_list[[i]][length(x_dim_list[[i]])] 
      if (i > 1) {
        if (length(x_dim_list[[i]]) != length(x_dim_list[[i-1]])) {
          stop("rds generator only works if separate inputs have same dimension size")
        }
      }
    }
    x_complete <- tensorflow::tf$concat(x_complete, 
                                        axis = as.integer(length(x_dim_list[[1]]) - 1)) %>% as.array()
  } 
  x_dim_start <- dim(x_complete)
  
  if (!is_lm) {
    if (multi_output) {
      y_dim_list <- list()
      size_splits_out <- list()
      for (i in 1:length(y_complete)) {
        y_dim_list[[i]] <- dim(y_complete[[i]])
        size_splits_out[[i]] <- y_dim_list[[i]][length(y_dim_list[[i]])] 
        if (i > 1) {
          if (length(y_dim_list[[i]]) != length(y_dim_list[[i-1]])) {
            stop("rds generator only works if separate outputs have same dimension size")
          }
        }
      }
      y_complete <- tensorflow::tf$concat(y_complete,
                                          axis = as.integer(length(y_dim_list[[1]]) - 1)) %>% as.array()
    } 
    y_dim_start <- dim(y_complete)
    if (is.null(y_dim_start)) y_dim_start <- length(y_complete)
    if (x_dim_start[1] != y_dim_start[1]) {
      stop("Different number of samples for input and target")
    }
  }
  sample_index <- 1:x_dim_start[1]
  
  if (!is.null(proportion_per_seq)) {
    sample_index <- sample(sample_index, min(length(sample_index), length(sample_index) * proportion_per_seq))
  }
  
  if (!is.null(max_samples)) {
    sample_index <- sample(sample_index, min(length(sample_index), max_samples))
  }
  
  if (!is.null(path_file_log)) {
    write.table(x = rds_files[1], file = path_file_log, row.names = FALSE, col.names = FALSE)
  }
  
  x_dim <- x_dim_start
  if (!is_lm) y_dim <- y_dim_start
  
  function() {
    
    # TODO: adjust for multi input/output
    x_index <- 1
    x <- array(0, c(batch_size, x_dim[2:3]))
    if (is_lm) {
      y <- vector("list", target_len)
    } else {
      y <- array(0, c(batch_size, y_dim[2]))
    }
    
    while (x_index <= batch_size) {
      
      # go to next file if sample index empty
      if (length(sample_index) == 0) {
        if (num_files == 1) {
          sample_index <<- 1:x_dim[1]
        } else {
          
          if (delete_used_files) file.remove(rds_files[file_index])
          
          read_success <- FALSE
          while (!read_success) {
            tryCatch(
              expr = {
                if (sample_by_file_size) {
                  file_index <<- sample(1:num_files, size = 1, prob = file_prob)
                } else {
                  file_index <<- file_index + 1
                  if (file_index > num_files) {
                    file_index <<- 1
                    rds_files <<- sample(rds_files)
                  }
                }
                rds_file <<- readRDS(rds_files[file_index])
                read_success <- TRUE
              },
              error = function(e){ 
                
              }
            )
          }
          
          if (multi_input) {
            # combine inputs in one tensor
            x_complete <<- tensorflow::tf$concat(rds_file[[1]], 
                                                 axis = as.integer(length(x_dim_list[[1]]) - 1)) %>% as.array()
          } else {
            x_complete <<- rds_file[[1]]
          } 
          x_dim <<- dim(x_complete)
          
          if (!is_lm) {
            if (multi_output) {
              y_complete <- tensorflow::tf$concat(rds_file[[2]], 
                                                  axis = as.integer(length(y_dim_list[[1]]) - 1)) %>% as.array()
            } else {
              y_complete <<- rds_file[[2]]
            }
            y_dim <<- dim(y_complete)
            #if (is.null(y_dim)) y_dim <<- length(y_complete)
          }
          
          if (!is_lm && (x_dim[1] != y_dim[1])) {
            print(x_dim)
            print(y_dim)
            stop("Different number of samples for input and target")
          }
          
          sample_index <<- 1:x_dim[1]
          if (!is.null(proportion_per_seq)) {
            if (length(sample_index) > 1) {
              sample_index <<- sample(sample_index, min(length(sample_index), max(1, floor(length(sample_index) * proportion_per_seq))))
            }
          }
          
          if (!is.null(max_samples)) {
            if (length(sample_index) > 1) {
              sample_index <<- sample(sample_index, min(length(sample_index), max_samples))
            }
          }
        }
        
        # log file
        if (!is.null(path_file_log)) {
          write.table(x = rds_files[file_index], file = path_file_log, append = TRUE, col.names = FALSE, row.names = FALSE)
        }
      }
      
      if (length(sample_index) == 0) next
      if (length(sample_index) == 1) {
        index <- sample_index
      } else {
        index <- sample(sample_index, min(batch_size - x_index + 1, length(sample_index)))
      }
      
      # subsetting
      x[x_index:(x_index + length(index) - 1), , ] <- x_complete[index, , ]
      
      if (!is_lm) {
        y[x_index:(x_index + length(index) - 1), ] <- y_complete[index, ]
      }
      
      x_index <- x_index + length(index)
      sample_index <<- setdiff(sample_index, index)
      
    }
    
    if (is_lm) {
      for (m in 1:target_len) {
        if (batch_size == 1) {
          y[[m]] <- matrix(x[ , x_dim[2] - target_len + m, ], nrow = 1)
        } else {
          y[[m]] <- x[ , x_dim[2] - target_len + m, ]
        }
      }
      
      x <- x[ , 1:(x_dim[2] - target_len), ]
      if (batch_size == 1) {
        dim(x) <- c(1, dim(x))
      }
    }
    
    if (!is.null(n_gram) & is_lm) {
      y <- do.call(rbind, y)
      y_list <- list()
      for (i in 1:batch_size) {
        index <- (i-1)  + (1 + (0:(target_len-1)) * batch_size)
        n_gram_matrix <- n_gram_of_matrix(input_matrix = y[index, ], n = n_gram)
        y_list[[i]] <- n_gram_matrix
      }
      y_tensor <- keras::k_stack(y_list, axis = 1L) %>% keras::k_eval()
      y <- vector("list", dim(y_tensor)[2])
      for (i in 1:dim(y_tensor)[2]) {
        y_subset <- y_tensor[ , i, ]
        if (batch_size == 1) y_subset <- matrix(y_subset, nrow = 1)
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
    
    if (!is.null(add_noise)) {
      noise_args <- c(add_noise, list(x = x))
      x <- do.call(add_noise_tensor, noise_args)
    }
    
    if (reverse_complement_encoding){
      x_1 <- x
      x_2 <- array(x_1[ , (dim(x)[2]):1, 4:1], dim = dim(x))
      x <- list(x_1, x_2)
    }
    
    if (multi_input) {
      x <- tensorflow::tf$split(x, num_or_size_splits = size_splits_in, axis = as.integer(length(x_dim)-1))
    }
    
    if (multi_output) {
      y <- tensorflow::tf$split(y, num_or_size_splits = size_splits_out, axis = as.integer(length(y_dim)-1))
    }
    
    return(list(x, y))
  }
}

#' Randomly select samples from fasta files
#' 
#' Generator \code{\link{generator_fasta_lm}}, \code{\link{generator_fasta_label_header_csv}}
#' or \code{\link{generator_fasta_label_folder}} will randomly choose a consecutive sequence of samples when
#' a \code{max_samples} argument is supplied. \code{generator_random} will choose samples at random.
#'
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_header_csv
#' @inheritParams train_model
#' @param number_target_nt Number of target nucleotides for language model.
#' @export
generator_random <- function(
  train_type = "label_folder",
  output_format = NULL,
  seed = 123,
  format = "fasta",
  reverse_complement = TRUE,
  path = NULL,
  batch_size = c(100),
  maxlen = 4,
  ambiguous_nuc = "equal",
  padding = FALSE,
  vocabulary = c("a", "c", "g", "t"),
  number_target_nt = 1,
  n_gram = NULL,
  n_gram_stride = NULL,
  sample_by_file_size = TRUE,
  max_samples = 1,
  skip_amb_nuc = NULL,
  vocabulary_label = NULL,
  target_from_csv = NULL,
  target_split = NULL,
  max_iter = 1000,
  verbose = TRUE,
  set_learning = NULL,
  shuffle_input = TRUE,
  reverse_complement_encoding = FALSE,
  proportion_entries = NULL,
  masked_lm = NULL,
  concat_seq = NULL,
  return_int = FALSE) {
  
  path_len <- ifelse(train_type != "label_folder", 1, length(path))
  vocabulary <- stringr::str_to_lower(vocabulary)
  vocabulary_label <- stringr::str_to_lower(vocabulary_label)
  label_from_header <- ifelse(train_type == "label_header", TRUE, FALSE)
  label_from_csv <- ifelse(train_type == "label_csv", TRUE, FALSE)
  
  if (ambiguous_nuc == "empirical") {
    stop("Empirical option not implemented for random sampling, only 'zero', 'equal' and 'discard'")
  }
  
  if (reverse_complement_encoding) {
    test_len <- length(vocabulary) != 4
    if (test_len || all(sort(stringr::str_to_lower(vocabulary)) != c("a", "c", "g", "t"))) {
      stop("reverse_complement_encoding only implemented for A,C,G,T vocabulary yet")
    }
  }
  
  if (length(batch_size) == 1 & (train_type == "label_folder")) {
    if ((batch_size %% path_len != 0)) {
      batch_size <- ceiling(batch_size/path_len) * path_len
      if (verbose) {
        message(paste("Batch size needs to be multiple of number of targets. Setting batch_size to", batch_size))
      }
    }
    batch_size <- rep(batch_size/path_len, path_len)
  }
  
  # set learning
  if (is.null(set_learning)) {
    samples_per_target <- NULL
    new_batch_size <- NULL
    reshape_mode <- NULL
    buffer_len <- NULL
  } else {
    if (train_type != "label_folder") {
      stop("train_type must be 'label_folder' when using set learning")
    }
    reshape_mode <- set_learning$reshape_mode
    samples_per_target <- set_learning$samples_per_target
    buffer_len <- set_learning$buffer_len
    maxlen <- set_learning$maxlen
    concat_maxlen <- NULL
    
    if (reshape_mode == "concat") {
      if (sum(batch_size) %% length(path) != 0) {
        stop_text <- paste("batch_size is", batch_size, "but needs to be multiple of number of classes (",
                           length(path), ") for set learning with 'concat'")
        stop(stop_text)
      }
      buffer_len <- ifelse(is.null(set_learning$buffer_len), 0, set_learning$buffer_len)
      concat_maxlen <- (maxlen * samples_per_target) + (buffer_len * (samples_per_target - 1))
      if (any(c("z", "Z") %in% vocabulary) & !is.null(set_learning$buffer_len)) {
        stop("'Z' is used as token for separating sequences and can not be in vocabulary.")
      }
      if (!is.null(set_learning$buffer_len)) {
        vocabulary <- c(vocabulary, "Z")
      }
      
    }
    
    if (any(batch_size[1] != batch_size)) {
      stop("Set learning only implemented for uniform batch_size for all classes.")
    }
    
    new_batch_size <- batch_size
    batch_size <- samples_per_target * batch_size
  }
  
  if (is.null(max_samples)) {
    stop("Please specify a max_samples argument when using random sampling")
  }
  
  if (length(max_samples) == 1 & (train_type == "label_folder")) {
    max_samples <- rep(max_samples, path_len)
  }
  
  if (any(max_samples > batch_size)) {
    message(paste("max_samples should not be bigger than batch_size when using random sampling, since generator opens a new file every batch"))
  }
  
  set.seed(seed)
  if (train_type == "label_folder") stopifnot(length(vocabulary_label) == path_len)
  target_len <- number_target_nt
  if (train_type != "lm") {
    number_target_nt <- NULL
    if (!is.null(n_gram)) stop("n-gram encoding not implemented yet for classification")
    seq_len_total <- maxlen
  } else {
    stopifnot(output_format %in% c("target_right", "target_middle_lstm", "target_middle_cnn", "wavenet"))
    if (train_type == "lm") stopifnot(length(batch_size) == 1)
    if (is.null(n_gram_stride)) n_gram_stride <- 1
    stopifnot(number_target_nt %% n_gram_stride == 0)
    seq_len_total <- maxlen + number_target_nt
    target_list <- list()
    if (is.null(n_gram) || n_gram == 1) {
      num_targets <- target_len/n_gram_stride
    } else {
      num_targets <- seq(1, (target_len - n_gram + 1), by = n_gram_stride) %>% length()
    }
    for (k in 1:num_targets) {
      target_list[[k]] <- list()
    }
  }
  
  count <- 1
  batch_number <- 1
  stop_gen <- FALSE
  
  # token for ambiguous nucleotides
  for (i in letters) {
    if (!(i %in% stringr::str_to_lower(vocabulary))) {
      amb_nuc_token <- i
      break
    }
  }
  pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
  tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), c(vocabulary, amb_nuc_token))
  if (label_from_header) {
    tokenizer_target <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = FALSE, lower = TRUE, filters = "\t\n"),
                                                  vocabulary_label)
  }
  #seq_list <- vector("list", sum(batch_size))
  seq_list <- vector("list", path_len)
  
  fasta_files <- list()
  num_files <- list()
  file_prob <- list()
  start_ind <- list()
  
  # target from csv
  if (label_from_csv) {
    .datatable.aware = TRUE
    output_label_csv <- utils::read.csv2(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
    if (dim(output_label_csv)[2] == 1) {
      output_label_csv <- read.csv(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
    }
    output_label_csv <- data.table::as.data.table(output_label_csv)
    stopifnot("file" %in% names(output_label_csv))
    data.table::setkey(output_label_csv, file)
    
    vocabulary_label <- names(output_label_csv)
    vocabulary_label <- vocabulary_label[vocabulary_label != "header" & vocabulary_label != "file"]
    if (!is.null(target_split)) {
      check_header_names(target_split = target_split, vocabulary_label = vocabulary_label)
    }
  }
  
  for (i in 1:path_len) {
    
    if (train_type == "label_folder") {
      fasta_files[[i]] <- list_fasta_files(path[[i]], format, file_filter = NULL)
    } else {
      fasta_files[[i]] <- list_fasta_files(path, format, file_filter = NULL)
    }
    
    # remove files without target label
    if (train_type == "label_csv") {
      if (any(duplicated(output_label_csv$file))) {
        stop("csv file with label contains duplicate file names in 'file' column")
      }
      fasta_files[[i]] <- fasta_files[[i]][basename(fasta_files[[i]]) %in% unique(output_label_csv$file)]
    }
    
    num_files[[i]] <- length(fasta_files[[i]])
    start_ind[[i]] <- (0:(batch_size[[i]] - 1) * seq_len_total) + 1
    if (sample_by_file_size) {
      file_prob[[i]] <- file.info(fasta_files[[i]])$size/sum(file.info(fasta_files[[i]])$size)
    } else {
      file_prob <- NULL
    }
  }
  
  # if (!sample_by_file_size & any(unlist(num_files) > 1)) {
  #   message("It is adviced to use sample_by_file_size when using random sampling strategy.")
  # }
  
  function() {
    
    for (p in 1:path_len) {
      remaining_samples <- batch_size[p]
      nuc_vector <- vector("character", remaining_samples)
      target_list <- list()
      target_count <- 1
      
      while (remaining_samples > 0) {
        
        nuc_seq <- ""
        length_vector <- 0
        iter <- 1
        while (all(length_vector < seq_len_total)) {
          file_index <- sample(1:num_files[[p]], size = 1, prob = file_prob[[p]])
          fasta_file <- read_fasta_fastq(format = format, skip_amb_nuc = skip_amb_nuc, file_index = file_index, pattern = pattern,
                                         shuffle_input = shuffle_input, proportion_entries = proportion_entries,
                                         vocabulary_label = vocabulary_label,
                                         filter_header = ifelse(train_type == "label_header", TRUE, FALSE),
                                         reverse_complement = reverse_complement, fasta.files = fasta_files[[p]])
          
          if (nrow(fasta_file) == 0) next
          if (!is.null(concat_seq)) {
            fasta_file <- data.frame(Header = "header", Sequence = paste(fasta_file$Sequence, collapse = stringr::str_to_lower(concat_seq)),
                                     stringsAsFactors = FALSE)
          }
          length_vector <- nchar(fasta_file$Sequence)
          if (!padding) {
            fasta_file <- fasta_file[length_vector >= seq_len_total, ]
          } else {
            short_seq_index <- which(length_vector < seq_len_total)
            for (i in short_seq_index) {
              fasta_file$Sequence[i] <- paste0(paste(rep("0", seq_len_total - length_vector[i]), collapse = ""), fasta_file$Sequence[i])
            }
          }
          nuc_seq <- fasta_file$Sequence
          if (!is.null(concat_seq)) {
            paste(nuc_seq, collapse = concat_seq)
            length_vector <- nchar(nuc_seq)
          } else {
            length_vector <- nchar(fasta_file$Sequence)
          }
          if (iter >= max_iter) {
            stop(paste("Could not extract sample for", iter, "iterations. Either maxlen is too big or sequences in fasta files too short."))
          }
          iter <- iter + 1
        }
        
        sample_start <- get_start_ind(seq_vector = nuc_seq,
                                      length_vector = length_vector,
                                      maxlen = seq_len_total,
                                      step = 1,
                                      train_mode = "label",
                                      discard_amb_nuc = ifelse(ambiguous_nuc == "discard", TRUE, FALSE),
                                      vocabulary = vocabulary)
        
        if (length(sample_start) == 0) next
        nuc_seq <- paste(nuc_seq, collapse = "")
        sample_start <- sample(sample_start, size = min(remaining_samples, max_samples[p], length(sample_start)))
        sample_end <- sample_start + seq_len_total - 1
        nuc_sample <- vector("list", length(sample_start))
        nuc_vector_start_index <- (batch_size[p] - remaining_samples + 1)
        nuc_vector_end_index <- nuc_vector_start_index + length(sample_start) - 1
        nuc_vector[nuc_vector_start_index : nuc_vector_end_index] <- unlist(purrr::map(1:length(sample_start),
                                                                                       ~substr(nuc_seq, sample_start[.x], sample_end[.x])))
        remaining_samples <- remaining_samples - length(sample_start)
        
        if (label_from_csv) {
          file_row_name <- fasta_files[[p]][file_index] %>% basename
          label_row <- output_label_csv[.(file_row_name)][ , -"file"] %>% as.matrix()
          label_matrix <- t(replicate(length(sample_start), label_row, simplify = TRUE))
          target_list[[target_count]] <- label_matrix %>% as.matrix()
          target_count <- target_count + 1
        }
        
        if (label_from_header) {
          start_new_entry <- c(1, cumsum(length_vector))
          start_new_entry <- start_new_entry[-length(start_new_entry)]
          label_row <- as.character(cut(sample_start, breaks = c(start_new_entry, sum(length_vector)),
                                        labels = fasta_file$Header, include.lowest = TRUE, right = FALSE))
          target_list[[target_count]] <- label_row
          target_count <- target_count + 1
        }
        
      }
      
      nuc_vector <- paste(nuc_vector, collapse = "")
      nuc_vector <- stringr::str_to_lower(nuc_vector)
      nuc_vector <- stringr::str_replace_all(string = nuc_vector, pattern = pattern, amb_nuc_token)
      nuc_vector <- keras::texts_to_sequences(tokenizer, nuc_vector)[[1]] - 1
      
      if (train_type != "lm") {
        one_hot_sample <- seq_encoding_label(sequence = nuc_vector,
                                             maxlen = maxlen,
                                             adjust_start_ind = TRUE,
                                             vocabulary = vocabulary,
                                             masked_lm = masked_lm, return_int = return_int,
                                             start_ind = start_ind[[p]],
                                             ambiguous_nuc = ambiguous_nuc, nuc_dist = NULL,
                                             quality_vector = NULL, use_coverage = FALSE,
                                             max_cov = NULL, n_gram_stride = n_gram_stride,
                                             cov_vector = NULL, n_gram = n_gram)
      } else {
        one_hot_sample <- seq_encoding_label(sequence = nuc_vector,
                                             maxlen = maxlen + target_len,
                                             vocabulary = vocabulary, 
                                             adjust_start_ind = TRUE,
                                             start_ind = start_ind[[p]],
                                             ambiguous_nuc = ambiguous_nuc,
                                             n_gram = n_gram, 
                                             n_gram_stride = n_gram_stride,
                                             output_format = output_format)
      }
      
      if (train_type != "lm") {
        seq_list[[p]] <- one_hot_sample
      } else {
        xy_list <- slice_tensor_lm(xy = xy,
                                   output_format = output_format,
                                   target_len = target_len,
                                   n_gram = n_gram,
                                   total_seq_len = total_seq_len,
                                   return_int = return_int)
      }
    }
    
    batch_number <<- batch_number + 1
    
    if (train_type == "lm") return(xy_list)
    
    if (train_type == "label_folder") {
      x <- abind::abind(seq_list, along = 1)
      y_list <- list()
      for (i in 1:length(batch_size)) {
        m <- matrix(0, ncol = length(batch_size), nrow = batch_size[i])
        m[ , i] <- 1
        y_list[[i]] <- m
      }
      y <- do.call(rbind, y_list)
    }
    
    if (train_type == "label_csv") {
      x <- one_hot_sample
      y <- do.call(rbind, target_list) %>% as.matrix()
      colnames(y) <- NULL
    }
    
    if (train_type == "label_header") {
      x <- one_hot_sample
      target_int <- unlist(keras::texts_to_sequences(tokenizer_target, unlist(target_list))) - 1
      y  <- keras::to_categorical(target_int, num_classes = length(vocabulary_label))
    }
    
    if (train_type == "masked_lm") {
      return(one_hot_sample)
    }
    
    if (reverse_complement_encoding){
      x_1 <- x
      x_2 <- array(x_1[ , (dim(x)[2]):1, 4:1], dim = dim(x))
      x <- list(x_1, x_2)
    }
    
    if (!is.null(set_learning)) {
      l <- reshape_tensor(x = x, y = y,
                          new_batch_size = sum(new_batch_size),
                          samples_per_target = samples_per_target,
                          buffer_len = buffer_len,
                          reshape_mode = reshape_mode)
      return(l)
    } else {
      return(list(X = x, Y = y))
    }
  }
}


#' Wrapper for generator functions
#' 
#' For a detailed description see the data generator [tutorial](https://deepg.de/articles/data_generator.html).
#' Will choose one of the generators from \code{\link{generator_fasta_lm}}, 
#' \code{\link{generator_fasta_label_folder}}, \code{\link{generator_fasta_label_header_csv}}, 
#' \code{\link{generator_rds}}, \code{\link{generator_random}}, \code{\link{generator_dummy}} or 
#' \code{\link{generator_fasta_lm}} according to the \code{train_type} and \code{random_sampling}
#' arguments.
#'
#' @inheritParams train_model
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_folder
#' @inheritParams generator_fasta_label_header_csv
#' @inheritParams generator_rds
#' @inheritParams generator_random
#' @inheritParams generator_initialize
#' @param path_file_logVal Path to csv file logging used validation files.
#' @param new_batch_size Batch size after reshaping data for set learning.
#' @examples 
#' # create dummy fasta files
#' fasta_path <- tempfile()
#' dir.create(fasta_path)
#' create_dummy_data(file_path = fasta_path,
#'                   num_files = 3,
#'                   seq_length = 10,
#'                   num_seq = 5,
#'                   vocabulary = c("a", "c", "g", "t"))
#' 
#' gen <- get_generator(path = fasta_path,
#'                      maxlen = 5, train_type = "lm",
#'                      output_format = "target_right",
#'                      step = 3, batch_size = 7)
#' z <- gen()
#' x <- z[[1]]
#' y <- z[[2]]
#' dim(x)
#' dim(y)
#' 
#' @export
get_generator <- function(path = NULL,
                          train_type,
                          batch_size,
                          maxlen = NULL,
                          step = NULL,
                          shuffle_file_order = FALSE,
                          vocabulary = c("A", "C", "G", "T"),
                          seed = 1,
                          proportion_entries = NULL,
                          shuffle_input = FALSE,
                          format = "fasta",
                          path_file_log = NULL,
                          reverse_complement = FALSE,
                          n_gram = NULL,
                          n_gram_stride = NULL,
                          output_format = "target_right",
                          ambiguous_nuc = "zero",
                          proportion_per_seq = NULL,
                          skip_amb_nuc = NULL,
                          use_quality_score = FALSE,
                          padding = FALSE,
                          added_label_path = NULL,
                          target_from_csv = NULL,
                          add_input_as_seq = NULL,
                          max_samples = NULL,
                          concat_seq = NULL,
                          target_len = 1,
                          file_filter = NULL,
                          use_coverage = NULL,
                          sample_by_file_size = FALSE,
                          add_noise = NULL,
                          random_sampling = FALSE,
                          set_learning = NULL,
                          file_limit = NULL,
                          reverse_complement_encoding = FALSE,
                          read_data = FALSE,
                          target_split = NULL,
                          path_file_logVal = NULL,
                          model = NULL,
                          vocabulary_label = NULL,
                          new_batch_size = NULL,
                          masked_lm = NULL,
                          val = FALSE,
                          return_int = FALSE,
                          delete_used_files = FALSE) {
  
  if (random_sampling) {
    if (use_quality_score) stop("use_quality_score not implemented for random sampling")
    if (read_data) stop("read_data not implemented for random sampling")
    if (!is.null(use_coverage)) stop("use_coverage not implemented for random sampling")
    if (!is.null(add_noise)) stop("add_noise not implemented for random sampling")
  }
  
  # adjust batch size
  if ((length(batch_size) == 1) && (batch_size %% length(path) != 0) & train_type == "label_folder") {
    batch_size <- ceiling(batch_size/length(path)) * length(path)
    if (!val) {
      message(paste("Batch size needs to be multiple of number of targets. Setting batch_size to", batch_size))
    }
  }
  
  if (is.null(step)) step <- maxlen
  
  if (train_type == "dummy_gen") {
    gen <- generator_dummy(model, ifelse(is.null(set_learning), batch_size, new_batch_size))
    removeLog <- FALSE
  }
  
  if (!is.null(added_label_path) & is.null(add_input_as_seq)) {
    add_input_as_seq <- rep(FALSE, length(added_label_path))
  }
  
  # language model
  if (train_type == "lm" & random_sampling) {
    
    gen <- generator_random(
      train_type = "lm",
      output_format = output_format,
      seed = seed[1],
      format = format,
      reverse_complement = reverse_complement,
      reverse_complement_encoding = reverse_complement_encoding,
      path = path,
      batch_size = batch_size,
      maxlen = maxlen,
      ambiguous_nuc = ambiguous_nuc,
      padding = padding,
      vocabulary = vocabulary,
      number_target_nt = target_len,
      target_split = target_split,
      target_from_csv = target_from_csv,
      n_gram = n_gram,
      n_gram_stride = n_gram_stride,
      sample_by_file_size = sample_by_file_size,
      max_samples = max_samples,
      skip_amb_nuc = skip_amb_nuc,
      vocabulary_label = vocabulary_label,
      shuffle_input = shuffle_input,
      proportion_entries = proportion_entries,
      return_int = return_int,
      concat_seq = concat_seq)
  } 
  
  if (train_type == "lm" & !random_sampling) {
    
    gen <- generator_fasta_lm(path_corpus = path, batch_size = batch_size,
                              maxlen = maxlen, step = step, shuffle_file_order = shuffle_file_order,
                              vocabulary = vocabulary, seed = seed[1], proportion_entries = proportion_entries,
                              shuffle_input = shuffle_input, format = format, n_gram_stride = n_gram_stride,
                              path_file_log = path_file_log, reverse_complement = reverse_complement, 
                              output_format = output_format, ambiguous_nuc = ambiguous_nuc,
                              proportion_per_seq = proportion_per_seq, skip_amb_nuc = skip_amb_nuc,
                              use_quality_score = use_quality_score, padding = padding, n_gram = n_gram,
                              added_label_path = added_label_path, add_input_as_seq = add_input_as_seq,
                              max_samples = max_samples, concat_seq = concat_seq, target_len = target_len,
                              file_filter = file_filter, use_coverage = use_coverage, 
                              sample_by_file_size = sample_by_file_size, add_noise = add_noise)
  }
  
  # label by folder
  if (train_type %in% c("label_folder", "masked_lm") & random_sampling) {
    
    gen <- generator_random(
      train_type = "label_folder",
      seed = seed[1],
      format = format,
      reverse_complement = reverse_complement,
      path = path,
      batch_size = batch_size,
      maxlen = maxlen,
      ambiguous_nuc = ambiguous_nuc,
      padding = padding,
      vocabulary = vocabulary,
      number_target_nt = NULL,
      n_gram = n_gram,
      n_gram_stride = n_gram_stride,
      sample_by_file_size = sample_by_file_size,
      max_samples = max_samples,
      skip_amb_nuc = skip_amb_nuc,
      shuffle_input = shuffle_input,
      set_learning = set_learning,
      reverse_complement_encoding = reverse_complement_encoding,
      vocabulary_label = vocabulary_label,
      proportion_entries = proportion_entries,
      masked_lm = masked_lm,
      return_int = return_int,
      concat_seq = concat_seq)
  } 
  
  if (train_type == "label_folder" & !random_sampling) {
    
    generator_initialize(directories = path, format = format, batch_size = batch_size, maxlen = maxlen, vocabulary = vocabulary,
                         verbose = FALSE, shuffle_file_order = shuffle_file_order, step = step, seed = seed[1],
                         shuffle_input = shuffle_input, file_limit = file_limit, skip_amb_nuc = skip_amb_nuc,
                         path_file_log = path_file_log, reverse_complement = reverse_complement,
                         reverse_complement_encoding = reverse_complement_encoding, return_int = return_int,
                         ambiguous_nuc = ambiguous_nuc, proportion_per_seq = proportion_per_seq,
                         read_data = read_data, use_quality_score = use_quality_score, val = val,
                         padding = padding, max_samples = max_samples, concat_seq = concat_seq,
                         added_label_path = added_label_path, add_input_as_seq = add_input_as_seq, use_coverage = use_coverage,
                         set_learning = set_learning, proportion_entries = proportion_entries,
                         sample_by_file_size = sample_by_file_size, n_gram = n_gram, n_gram_stride = n_gram_stride,
                         add_noise = add_noise)
    
    gen <- generator_fasta_label_folder_wrapper(val = val, path = path,  new_batch_size = new_batch_size,
                                                batch_size = batch_size, voc_len = length(vocabulary),
                                                maxlen = maxlen, set_learning = set_learning)
    
  }
  
  if (train_type == "masked_lm" & !random_sampling) {
    
    stopifnot(!is.null(masked_lm))
    
    gen <- generator_fasta_label_folder(path_corpus = unlist(path),
                                        format = format,
                                        batch_size = batch_size,
                                        maxlen = maxlen,
                                        vocabulary = vocabulary,
                                        shuffle_file_order = shuffle_file_order,
                                        step = step,
                                        seed = seed,
                                        shuffle_input = shuffle_input,
                                        file_limit = file_limit,
                                        path_file_log = path_file_log,
                                        reverse_complement = reverse_complement,
                                        reverse_complement_encoding = reverse_complement_encoding,
                                        num_targets = 1,
                                        ones_column = 1,
                                        ambiguous_nuc = ambiguous_nuc,
                                        proportion_per_seq = proportion_per_seq,
                                        read_data = read_data,
                                        use_quality_score = use_quality_score,
                                        padding = padding,
                                        added_label_path = added_label_path,
                                        add_input_as_seq = add_input_as_seq,
                                        skip_amb_nuc = skip_amb_nuc,
                                        max_samples = max_samples,
                                        concat_seq = concat_seq,
                                        file_filter = NULL,
                                        return_int = return_int,
                                        use_coverage = use_coverage,
                                        proportion_entries = proportion_entries,
                                        sample_by_file_size = sample_by_file_size,
                                        n_gram = n_gram,
                                        n_gram_stride = n_gram_stride,
                                        masked_lm = masked_lm,
                                        add_noise = add_noise) 
  }
  
  
  if ((train_type == "label_csv" | train_type == "label_header") & !random_sampling) {
    
    gen <- generator_fasta_label_header_csv(path_corpus = path, format = format, batch_size = batch_size, maxlen = maxlen,
                                            vocabulary = vocabulary, verbose = FALSE, shuffle_file_order = shuffle_file_order, step = step,
                                            seed = seed[1], shuffle_input = shuffle_input, return_int = return_int,
                                            path_file_log = path_file_log, vocabulary_label = vocabulary_label, reverse_complement = reverse_complement,
                                            ambiguous_nuc = ambiguous_nuc, proportion_per_seq = proportion_per_seq,
                                            read_data = read_data, use_quality_score = use_quality_score, padding = padding,
                                            added_label_path = added_label_path, add_input_as_seq = add_input_as_seq,
                                            skip_amb_nuc = skip_amb_nuc, max_samples = max_samples, concat_seq = concat_seq,
                                            target_from_csv = target_from_csv, target_split = target_split, file_filter = file_filter,
                                            use_coverage = use_coverage, proportion_entries = proportion_entries,
                                            sample_by_file_size = sample_by_file_size, n_gram = n_gram, n_gram_stride = n_gram_stride,
                                            add_noise = add_noise, reverse_complement_encoding = reverse_complement_encoding)
  }
  
  if ((train_type == "label_csv" | train_type == "label_header") & random_sampling) {
    
    gen <- generator_random(
      train_type = train_type, 
      output_format = output_format,
      seed = seed[1],
      format = format,
      reverse_complement = reverse_complement,
      reverse_complement_encoding = reverse_complement_encoding,
      path = path,
      batch_size = batch_size,
      maxlen = maxlen,
      ambiguous_nuc = ambiguous_nuc,
      padding = padding,
      vocabulary = vocabulary,
      number_target_nt = NULL,
      n_gram = n_gram,
      n_gram_stride = n_gram_stride,
      sample_by_file_size = sample_by_file_size,
      max_samples = max_samples,
      skip_amb_nuc = skip_amb_nuc,
      vocabulary_label = vocabulary_label,
      target_from_csv = target_from_csv,
      target_split = target_split,
      verbose = verbose,
      shuffle_input = shuffle_input,
      proportion_entries = proportion_entries,
      return_int = return_int,
      concat_seq = concat_seq)
  }
  
  if (train_type %in% c("label_rds", "lm_rds")) {
    reverse_complement <- FALSE
    step <- 1
    if (train_type == "label_rds") target_len <- NULL
    gen <- generator_rds(rds_folder = path, batch_size = batch_size, path_file_log = path_file_log,
                         max_samples = max_samples, proportion_per_seq = proportion_per_seq,
                         sample_by_file_size = sample_by_file_size, add_noise = add_noise,
                         reverse_complement_encoding = reverse_complement_encoding, seed = seed[1],
                         target_len = target_len, n_gram = n_gram, n_gram_stride = n_gram_stride,
                         delete_used_files = delete_used_files)
    
  }
  
  return(gen)
  
}


#' Collect samples from generator and store in rds or pickle file.
#'
#' Repeatedly generate samples with data generator and store output. Creates a separate rds or pickle file in \code{output_path} for each 
#' batch.
#'
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_header_csv
#' @inheritParams generator_initialize
#' @inheritParams generator_fasta_label_folder_wrapper
#' @inheritParams get_generator
#' @inheritParams train_model
#' @param iterations Number of batches (output files) to create.
#' @param output_path Output directory. Output files will be named `output_path` + `file_name_start` + x + ".rds" or ".pickle", where x is an index (from 1 to 
#' \code{iterations}) and file ending depends on \code{as_numpy_array} argument.
#' @param as_numpy_array Store output as list of numpy arrays if `TRUE` (otherwise as R array). 
#' @param shuffle Whether to shuffle samples within each batch.
#' @param file_name_start Start of output file names.
#' @param store_format Either "rds" or "pickle". 
#' @param ... further generator options. See \code{\link{get_generator}}.
#' @examples 
#' # create dummy fasta files
#' temp_dir <- tempfile()
#' dir.create(temp_dir)
#' create_dummy_data(file_path = temp_dir,
#'                   num_files = 3,
#'                   seq_length = 8, 
#'                   num_seq = 2)
#' 
#' # extract samples
#' out_dir <- tempfile()
#' dir.create(out_dir)
#' dataset_from_gen(output_path = out_dir,
#'                  iterations = 10,
#'                  train_type = "lm",
#'                  output_format = "target_right",
#'                  path_corpus = temp_dir, 
#'                  batch_size = 32,
#'                  maxlen = 5,
#'                  step = 1,
#'                  file_name_start = "batch_")
#' 
#' list.files(out_dir)
#' @export
dataset_from_gen <- function(output_path,
                             iterations = 10,
                             train_type = "lm",
                             output_format = "target_right",
                             path_corpus, 
                             batch_size = 32,
                             maxlen = 250,
                             step = NULL,
                             vocabulary = c("a", "c", "g", "t"),
                             shuffle = FALSE,
                             set_learning = NULL,
                             seed = NULL,
                             random_sampling,
                             store_format = "rds",
                             file_name_start = "batch_",
                             masked_lm = NULL,
                             ...) {
  
  stopifnot(train_type %in% c("lm", "label_header", "label_folder", "label_csv", "masked_lm", "dummy_gen"))
  stopifnot(store_format %in% c("pickle", "rds"))
  include_sample_weights <- !is.null(masked_lm) && masked_lm$include_sw
  
  if (is.null(step)) step <- maxlen
  if (is.null(seed)) seed <- sample(1:100000, 1)
  
  gen <- get_generator(path = path_corpus,
                       val = FALSE,
                       batch_size = batch_size,
                       maxlen = maxlen,
                       step = step,
                       model = NULL,
                       vocabulary = vocabulary,
                       file_filter = NULL,
                       train_type = train_type,
                       set_learning = set_learning,
                       path_file_logVal = NULL,
                       seed = seed,
                       masked_lm = masked_lm,
                       ...)
  
  for (batch_number in 1:iterations) {
    
    z <- gen()
    x <- z[[1]]
    y <- z[[2]]
    if (include_sample_weights) sw <- z[[3]]
    
    if (shuffle) {
      shuffle_index <- sample(dim(x)[1])
      x <- shuffle_batches(x, shuffle_index)
      y <- shuffle_batches(y, shuffle_index)
      if (include_sample_weights) sw <- shuffle_batches(sw, shuffle_index)
    }
    
    base_path <- file.path(output_path, file_name_start)
    filename <- paste0(base_path, batch_number, ".", store_format)
    if (include_sample_weights) {
      out_list <- list(x, y, sw)
    } else {
      out_list <- list(x, y)
    } 
    if (store_format == "pickle") {
      reticulate::py_save_object(object = reticulate::r_to_py(out_list), filename = filename)
    }
    if (store_format == "rds") {
      saveRDS(out_list, file = filename)
    }
    # if (store_format == "hdf5") {
    #   saveRDS(out_list, file = filename)
    # }
    
  }  
  
}
