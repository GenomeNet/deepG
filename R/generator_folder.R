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
#' \item `block_len` (optional): Masked/random/identity regions appear in blocks of size `block_len`.     
#' }
#' @examples
#' # create dummy fasta files
#' path_input_1 <- tempfile()
#' dir.create(path_input_1)
#' create_dummy_data(file_path = path_input_1, 
#'                   num_files = 2,
#'                   seq_length = 7,
#'                   num_seq = 1,
#'                   vocabulary = c("a", "c", "g", "t"))
#' 
#' gen <- generator_fasta_label_folder(path_corpus = path_input_1, batch_size = 2,
#'                                     num_targets = 3, ones_column = 2, maxlen = 7)
#' z <- gen()
#' dim(z[[1]])
#' z[[2]]
#' 
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
                                         return_int = FALSE,
                                         reshape_xy = NULL) {
  
  if (!is.null(reshape_xy)) {
    reshape_xy_bool <- TRUE
    reshape_x_bool <- ifelse(is.null(reshape_xy$x), FALSE, TRUE)
    reshape_y_bool <- ifelse(is.null(reshape_xy$y), FALSE, TRUE)
  } else {
    reshape_xy_bool <- FALSE
  }
  
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
    utils::write.table(x = fasta.files[file_index], file = path_file_log, col.names = FALSE, row.names = FALSE, append = append)
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
      if (reshape_xy_bool) {
        if (reshape_x_bool) x <- reshape_xy$x(x)
        if (reshape_y_bool) y <- reshape_xy$y(y)
      }
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
    
    if (reshape_xy_bool) {
      if (reshape_x_bool) x <- reshape_xy$x(x)
      if (reshape_y_bool) y <- reshape_xy$y(y)
    }
    
    if (!is.null(masked_lm) && include_sw) return(list(x, y, sample_weight))
    return(list(X = x, Y = y))
    
  }
}
