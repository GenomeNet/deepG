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
#' @param read_data If `TRUE` the first element of input is a list of length 2, each containing one part of paired read. Maxlen should be 2*length of one read.
#' @rawNamespace import(data.table, except = c(first, last, between))
#' @examples
#' path_input <- tempfile()
#' dir.create(path_input)
#' # create 2 fasta files called 'file_1.fasta', 'file_2.fasta'
#' create_dummy_data(file_path = path_input, 
#'                   num_files = 2,
#'                   seq_length = 5,
#'                   num_seq = 1,
#'                   vocabulary = c("a", "c", "g", "t"))
#' dummy_labels <- data.frame(file = c('file_1.fasta', 'file_2.fasta'), # dummy labels
#'                            label1 = c(0, 1),
#'                            label2 = c(1, 0))
#' target_from_csv <- tempfile(fileext = '.csv')
#' write.csv(dummy_labels, target_from_csv, row.names = FALSE)
#' gen <- generator_fasta_label_header_csv(path_corpus = path_input, batch_size = 2, 
#'                                         maxlen = 5, target_from_csv = target_from_csv)
#' z <- gen()
#' dim(z[[1]])
#' z[[2]]
#' 
#' @returns A generator function.  
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
                                             return_int = FALSE,
                                             reshape_xy = NULL) {
  
  if (!is.null(reshape_xy)) {
    reshape_xy_bool <- TRUE
    reshape_x_bool <- ifelse(is.null(reshape_xy$x), FALSE, TRUE)
    if (reshape_x_bool && !all(c('x', 'y') %in% names(formals(reshape_xy$x)))) {
      stop("function reshape_xy$x needs to have arguments named x and y")
    }
    reshape_y_bool <- ifelse(is.null(reshape_xy$y), FALSE, TRUE)
    if (reshape_y_bool && !all(c('x', 'y') %in% names(formals(reshape_xy$y)))) {
      stop("function reshape_xy$y needs to have arguments named x and y")
    }
  } else {
    reshape_xy_bool <- FALSE
  }
  
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
  
  use_basename <- TRUE
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
  
  # target from csv
  if (!is.null(target_from_csv)) {
    .datatable.aware = TRUE
    if (!is.data.frame(target_from_csv)) {
      output_label_csv <- utils::read.csv2(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
      if (dim(output_label_csv)[2] == 1) {
        output_label_csv <- utils::read.csv(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
      }
    } else {
      output_label_csv <- target_from_csv
    }
    output_label_csv <- data.table::as.data.table(output_label_csv)
    if ("file" %in% names(output_label_csv) & "header" %in% names(output_label_csv)) {
      stop('names in target_from_csv should contain "header" or "file" not both')
    } else if ("header" %in% names(output_label_csv)) {
      # added_label_by_header_target <- TRUE
      # data.table::setkey(output_label_csv, header)
    } else if ("file" %in% names(output_label_csv)) {
      added_label_by_header_target <- FALSE
      #output_label_csv$file <- stringr::str_to_lower(as.character(output_label_csv$file))
      data.table::setkey(output_label_csv, file)
    } else {
      stop('file in target_from_csv must contain one column named "header" or "file"')
    }
    
    # remove files without target label
    if (!added_label_by_header_target) {
      # relative path or absolute path
      if (dirname(output_label_csv$file[1]) == ".") {
        index_basename <- basename(fasta.files) %in% output_label_csv$file
      } else {
        use_basename <- FALSE
        index_basename <- fasta.files %in% output_label_csv$file
      }
      index_abs_path <- fasta.files %in% output_label_csv$file
      index <- index_basename | index_abs_path
      fasta.files <- fasta.files[index]
      if (length(fasta.files) == 0) {
        stop("No overlap between files and 'file' column in target_from_csv")
      }
    }
    
    # if (!added_label_by_header_target) {
    #   fasta.file$Header <- rep(basename(fasta.files[file_index]), nrow(fasta.file))
    # }
    col_name <- ifelse(added_label_by_header_target, "header", "file")
    #header_vector <- NULL
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
    utils::write.table(x = fasta.files[1], file = path_file_log, row.names = FALSE, col.names = FALSE)
  }
  
  rngstate <- .GlobalEnv$.Random.seed
  
  function() {
    
    .GlobalEnv$.Random.seed <- rngstate
    on.exit(rngstate <<- .GlobalEnv$.Random.seed)
    # loop until enough samples collected
    while(num_samples < batch_size) {
      iter <- 1
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
          utils::write.table(x = fasta.files[file_index], file = path_file_log, append = TRUE, col.names = FALSE, row.names = FALSE)
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
          ff <- fasta.files[file_index]
          if (use_basename) ff <- basename(ff)
          label_list[[sequence_list_index]] <- ff
        }
      }
      
      if (!is.null(target_from_csv)) {
        if (added_label_by_header_target) {
          output_label_list[[sequence_list_index]] <- as.character(cut(subsetStartIndices, breaks = c(startNewEntry, length(nucSeq)),
                                                                       labels = header_vector, include.lowest = TRUE, right = FALSE))
        } else {
          ff <- fasta.files[file_index]
          if (use_basename) ff <- basename(ff)
          output_label_list[[sequence_list_index]] <- ff
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
          index_output_label_vector <- output_label_vector == i
          if (nrow(output_label_from_csv) > 0) {
            target_tensor[index_output_label_vector, ] <- matrix(as.matrix(output_label_from_csv[1, ]),
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
          
          if (nrow(stats::na.omit(output_label_from_csv)) > 0) {
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
    
    if (reshape_xy_bool) {
      l <- f_reshape(x = x, y = y,
                     reshape_xy = reshape_xy,
                     reshape_x_bool = reshape_x_bool,
                     reshape_y_bool = reshape_y_bool,
                     reshape_sw_bool = FALSE, sw = NULL)
      return(l)
    }
    
    return(list(X = x, Y = y))
  }
}
