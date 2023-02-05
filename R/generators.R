#' Custom generator for fasta/fastq files
#'
#' @description
#' \code{fastaFileGenerator} Iterates over folder containing .fasta/.fastq files and produces one-hot-encoding of predictor sequences 
#' and target variables.
#' 
#' @param corpus.dir Input directory where .fasta files are located or path to single file ending with .fasta or .fastq 
#' (as specified in format argument). Can also be a list of directories. 
#' @param format File format, either fasta or fastq.
#' @param batch.size Number of batches.    
#' @param maxlen Length of predictor sequence.  
#' @param max_iter Stop after max_iter number of iterations failed to produce a new batch. 
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as 0-vector.
#' @param randomFiles Logical, whether to go through files randomly or sequential. 
#' @param step How often to take a sample.
#' @param showWarnings Logical, give warning if character outside vocabulary appears   
#' @param seed Sets seed for set.seed function, for reproducible results when using \code{randomFiles} or \code{shuffleFastaEntries}  
#' @param shuffleFastaEntries Logical, shuffle entries in every fasta file before connecting them to sequence.
#' @param verbose Whether to show message. 
#' @param numberOfFiles Use only specified number of files, ignored if greater than number of files in corpus.dir. 
#' @param fileLog Write name of files to csv file if path is specified.
#' @param reverseComplements Logical, for every new file decide randomly to use original data or its reverse complement.
#' @param output_format Determines shape of output tensor for language model.
#' Either "target_right", "target_middle_lstm", "target_middle_cnn" or "wavenet".
#' Assume a sequence "AACCGTA". Output correspond as follows
#' "target_right": X = "AACCGT", Y = "A"
#' "target_middle_lstm": X = (X_1 = "AAC", X_2 = "ATG"), Y = "C" (note reversed order of X_2)
#' "target_middle_cnn": X = "AACGTA", Y = "C" 
#' "wavenet": X = "AACCGT", Y = "ACCGTA"
#' @param ambiguous_nuc How to handle nucleotides outside vocabulary, either "zero", "discard", "empirical" or "equal". If "zero", input gets encoded as zero vector; 
#' if "equal" input is 1/length(vocabulary) x length(vocabulary). If "discard" samples containing nucleotides outside vocabulary get discarded. 
#' If "empirical" use nucleotide distribution of current file.   
#' @param proportion_per_file Numerical value between 0 and 1. Proportion of possible samples to take from one file. Takes samples from random subsequence.   
#' @param use_quality_score Whether to use fastq qualitiy scores. If TRUE input is not one-hot-encoding but corresponds to probabilities.
#' For example (0.97, 0.01, 0.01, 0.01) instead of (1, 0, 0, 0).   
#' @param padding Whether to pad sequences too short for one sample with zeros. 
#' @param added_label_path Path to folder with additional labels. Should be a csv file with one column named "file". Other columns should correspond to one label.
#' @param add_input_as_seq Boolean vector specifying for each entry in \code{added_label_path} if rows from csv should be encoded as a sequence or used directly.
#' If a row in your csv file is a sequence this should be true. For example you may want to add another sequence, say ACCGT. Then this would correspond to 1,2,2,3,4 in
#' csv file (if vocabulary = c("A", "C", "G", "T")).  If \code{add_input_as_seq} is TRUE, 12234 gets one-hot encoded, so added input is a 3D tensor.  If \code{add_input_as_seq} is 
#' FALSE this will feed network just raw data (a 2D tensor).
#' @param skip_amb_nuc Threshold of ambiguous nucleotides to accept in fasta entry. Complete entry will get discarded otherwise.  
#' @param max_samples Maximum number of samples to use from one file. If not NULL and file has more than \code{max_samples} samples, will randomly choose a 
#' subset of \code{max_samples} samples. 
#' @param concat_seq Character string or NULL. If not NULL all entries from file get concatenated to one sequence with concat_seq string between them.
#' Example: If 1.entry AACC, 2. entry TTTG and concat_seq = "ZZZ" this becomes AACCZZZTTTG.
#' @param target_len Number of nucleotides to predict at once. 
#' @param file_filter Vector of file names to use from corpus.dir.
#' @param use_coverage Integer or NULL. If not NULL, use coverage as encoding rather than one-hot encoding and normalize . 
#' @param proportion_entries Proportion of fasta entries to keep. For example, if fasta file has 50 entries and proportion_entries = 0.1, 
#' will randomly select 5 entries.
#' @param sample_by_file_size Sample new file weighted by file size (possible to repeatedly sample the same file). 
#' @param n_gram Integer, encode target not nucleotide wise but combine n nucleotides at once. For example for n=2, "AA" ->  (1, 0,..., 0),
#' "AC" ->  (0, 1, 0,..., 0), "TT" -> (0,..., 0, 1), where the one-hot vector have length length(vocabulary)^n.   
#' @import data.table 
#' @export
fastaFileGenerator <- function(corpus.dir,
                               format = "fasta",
                               batch.size = 256,
                               maxlen = 250,
                               max_iter = 10000,
                               vocabulary = c("a", "c", "g", "t"),
                               verbose = FALSE,
                               randomFiles = FALSE,
                               step = 1, 
                               showWarnings = FALSE,
                               seed = 1234,
                               shuffleFastaEntries = FALSE,
                               numberOfFiles = NULL,
                               fileLog = NULL,
                               reverseComplements = FALSE,
                               output_format = "target_right",
                               ambiguous_nuc = "zeros",
                               use_quality_score = FALSE,    
                               proportion_per_file = NULL,
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
                               n_gram_stride = 1) {
  
  if (!is.null(concat_seq) && (!all(stringr::str_split(concat_seq,"")[[1]] %in% vocabulary))) {
    stop("Characters of separating sequence should be in vocabulary")
  }
  if (is.null(use_coverage)) {
    use_coverage <- FALSE
    max_cov <- NULL
  } else {
    max_cov <- use_coverage
    use_coverage <- TRUE
  }
  stopifnot(output_format %in% c("target_right", "target_middle_lstm", "target_middle_cnn", "wavenet"))
  wavenet_format <- FALSE; target_middle <- FALSE; cnn_format <- FALSE
  if (output_format == "target_middle_lstm") target_middle <- TRUE 
  if (output_format == "target_middle_cnn") cnn_format <- TRUE 
  if (output_format == "wavenet") wavenet_format <- TRUE
  additional_labels <- ifelse(is.null(added_label_path), FALSE, TRUE)
  
  # token for ambiguous nucleotides
  for (i in letters) {
    if (!(i %in% stringr::str_to_lower(vocabulary))) {
      amb_nuc_token <- i
      break
    }
  }
  
  # adjust maxlen
  original_maxlen <- maxlen
  #if (is.null(n_gram)) {
  if (target_len > 1) {
    maxlen <- maxlen + target_len - 1
  }
  #} else {
  #  target_len <- target_len + n_gram - 1
  #}
  
  set.seed(seed)
  vocabulary <- stringr::str_to_lower(vocabulary)
  discard_amb_nuc <- ifelse(ambiguous_nuc == "discard", TRUE, FALSE)
  tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), c(vocabulary, amb_nuc_token))
  start_index_list <- vector("list")
  file_index <- 1
  num_samples <- 0
  start_index <- 1
  iter <- 1
  concat <- !is.null(concat_seq)
  seq_vector <- NULL
  
  fasta.files <- list_fasta_files(corpus.dir = corpus.dir,
                                  format = format,
                                  file_filter = file_filter)
  num_files <- length(fasta.files)
  
  if (sample_by_file_size) {
    randomFiles <- FALSE
    file_prob <- file.info(fasta.files)$size/sum(file.info(fasta.files)$size)
  }
  
  if (randomFiles) fasta.files <- sample(fasta.files, replace = FALSE)
  
  # regular expression for chars outside vocabulary
  pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
  
  # sequence vector collects strings until one batch can be created   
  sequence_list <- vector("list")
  quality_list <- vector("list") 
  coverage_list <- vector("list")
  
  if (!use_quality_score) {
    quality_list <- NULL
  }
  sequence_list_index <- 1
  
  while (length(seq_vector) == 0) {
    
    fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                   shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                   reverseComplements = reverseComplements, fasta.files = fasta.files)
    
    if (concat) {
      cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
      fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                               stringsAsFactors = FALSE)
    }
    
    # skip file that can't produce one sample
    if (!padding) {
      while((nrow(fasta.file) ==0) || all(nchar(fasta.file$Sequence) < (maxlen + 1))) {
        file_index <- file_index + 1
        iter <- iter + 1
        if (file_index > length(fasta.files) || iter > max_iter) {
          stop("Can not extract enough samples, try reducing maxlen parameter")
        }
        fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                       shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                       reverseComplements = reverseComplements, fasta.files = fasta.files)
        
        if (concat) {
          cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
          fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                                   stringsAsFactors = FALSE)
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
                                       shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                       reverseComplements = reverseComplements, fasta.files = fasta.files)
        
        if (concat) {
          cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
          fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                                   stringsAsFactors = FALSE)
        }
      }
    }
    
    if (use_coverage) {
      cov_vector <- get_coverage(fasta.file)
    }
    
    # take random subset
    if (!is.null(proportion_per_file)) {
      fasta_width <- nchar(fasta.file$Sequence)
      sample_range <- floor(fasta_width - (proportion_per_file * fasta_width))
      start <- mapply(sample_range, FUN = sample, size = 1)
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
      #added_label_by_header <- ifelse(added_label_list[[1]]$col_name == "header", TRUE, FALSE)
      added_label_by_header <- FALSE
    }
    
    # pad short sequences with zeros or discard
    short_seq_index <- which(length_vector < (maxlen + 1))
    if (padding) {
      for (i in short_seq_index) {
        seq_vector[i] <- paste0(paste(rep("0", (maxlen + 1) - length_vector[i]), collapse = ""), seq_vector[i])
        if (use_quality_score) {
          quality_scores[i] <- paste0(paste(rep("!", (maxlen + 1) - length_vector[i]), collapse = ""), quality_scores[i])
        }
        length_vector[i] <- maxlen + 1
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
        if (randomFiles) fasta.files <- sample(fasta.files, replace = FALSE)
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
  } else {
    nuc_dist <- 0
  }
  
  # start positions of samples
  start_indices <- getStartInd(seq_vector = seq_vector, length_vector = length_vector, maxlen = maxlen, step = step, train_mode = "lm",
                               discard_amb_nuc = discard_amb_nuc, vocabulary = vocabulary)
  
  # limit samples per file 
  if (!is.null(max_samples) && length(start_indices) > max_samples) {
    max_samples_subsample <- sample(1:(length(start_indices) - max_samples + 1), 1)
    start_indices <- start_indices[max_samples_subsample:(max_samples_subsample + max_samples - 1)]
  }
  
  startNewEntry <- cumsum(c(1, length_vector[-length(length_vector)]))
  
  nucSeq <- stringr::str_to_lower(nucSeq) 
  nucSeq <- keras::texts_to_sequences(tokenizer, nucSeq)[[1]] - 1
  
  # use subset of files
  if (!is.null(numberOfFiles) && (numberOfFiles < length(fasta.files))) {
    fasta.files <- fasta.files[1:numberOfFiles]
    num_files <- length(fasta.files)
  }
  
  # log file
  if (!is.null(fileLog)) {
    if (!endsWith(fileLog, ".csv")) fileLog <- paste0(fileLog, ".csv")
    # if (file.exists(fileLog)) {
    #   write.table(x = filePath, file = fileLog, append = TRUE, col.names = FALSE, row.names = FALSE)
    # } else {
    write.table(x = fasta.files[1], file = fileLog, row.names = FALSE, col.names = FALSE)
    # }
  }
  
  # test for chars outside vocabulary
  if (showWarnings) {
    charsOutsideVoc <- stringr::str_detect(stringr::str_to_lower(nucSeq), pattern)  
    if (charsOutsideVoc) warning("file ", filePath, " contains characters outside vocabulary")
  }
  
  if (verbose) message("Initializing ...")
  
  rngstate <- .GlobalEnv$.Random.seed
  
  function() {
    
    .GlobalEnv$.Random.seed <- rngstate
    on.exit(rngstate <<- .GlobalEnv$.Random.seed)
    iter <- 1
    # loop until enough samples collected
    while(num_samples < batch.size) {  
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
          if (randomFiles) fasta.files <<- sample(fasta.files, replace = FALSE)
          file_index <<- 1
        }
        
        filePath <<- fasta.files[[file_index]]
        
        # skip empty files
        while(TRUE) {
          fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                         shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                         reverseComplements = reverseComplements, fasta.files = fasta.files)
          
          if (concat) {
            cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
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
            if (randomFiles) fasta.files <<- sample(fasta.files, replace = FALSE)
            file_index <<- 1
          }
        }
        
        if (use_coverage) {
          cov_vector <<- get_coverage(fasta.file)
        }
        
        # log file
        if (!is.null(fileLog)) {
          write.table(x = filePath, file = fileLog, append = TRUE, col.names = FALSE, row.names = FALSE)
        }
        
        # test for chars outside vocabulary
        if (showWarnings) {
          charsOutsideVoc <- stringr::str_detect(stringr::str_to_lower(nucSeq), pattern)  
          if (charsOutsideVoc) warning("file ", filePath, " contains characters outside vocabulary")
        }
        
        # take random subset
        if (!is.null(proportion_per_file)) {
          fasta_width <- nchar(fasta.file$Sequence)
          sample_range <- floor(fasta_width - (proportion_per_file * fasta_width))
          start <- mapply(sample_range, FUN = sample, size = 1)
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
        }
        
        seq_vector <- stringr::str_to_lower(seq_vector)
        seq_vector <- stringr::str_replace_all(string = seq_vector, pattern = pattern, amb_nuc_token)
        length_vector <- nchar(seq_vector)
        if (additional_labels) {
          if (!added_label_by_header) {
            fasta.file$Header <- rep(basename(fasta.files[file_index]), nrow(fasta.file))
          }
          header_vector <<- fasta.file$Header
        }
        
        # pad short sequences with zeros or discard
        short_seq_index <- which(length_vector < (maxlen + 1))
        if (padding) {
          for (i in short_seq_index) {
            seq_vector[i] <- paste0(paste(rep("0", (maxlen + 1) - length_vector[i]), collapse = ""), seq_vector[i])
            if (use_quality_score) {
              quality_scores[i] <- paste0(paste(rep("!", (maxlen + 1) - length_vector[i]), collapse = ""), quality_scores[i])
            }
            length_vector[i] <- maxlen + 1
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
              header_vector <<- header_vector[-short_seq_index]
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
        
        # start positions of samples
        start_indices <<- getStartInd(seq_vector = seq_vector, length_vector = length_vector, maxlen = maxlen, step = step, train_mode = "lm",
                                      discard_amb_nuc = discard_amb_nuc, vocabulary = vocabulary)
        
        # limit samples per file 
        if (!is.null(max_samples) && length(start_indices) > max_samples) {
          max_samples_subsample <- sample(1:(length(start_indices) - max_samples + 1), 1)
          start_indices <<- start_indices[max_samples_subsample:(max_samples_subsample + max_samples - 1)]
        }
        
        startNewEntry <<- cumsum(c(1, length_vector[-length(length_vector)]))
        nucSeq <<- stringr::str_to_lower(nucSeq) 
        nucSeq <<- keras::texts_to_sequences(tokenizer, nucSeq)[[1]] - 1
        
        if(iter > max_iter) {
          stop('exceeded max_iter value, try reducing maxlen parameter')
          break
        }
        iter <- iter + 1
      }
      
      # go as far as possible in sequence or stop when enough samples are collected 
      remainingSamples <- batch.size - num_samples
      end_index <- min(length(start_indices), start_index + remainingSamples  - 1)
      
      subsetStartIndices <- start_indices[start_index:end_index]
      sequence_list[[sequence_list_index]] <- nucSeq[subsetStartIndices[1]:(subsetStartIndices[length(subsetStartIndices)] + maxlen)]
      if (use_quality_score) {
        quality_list[[sequence_list_index]] <- quality_vector[subsetStartIndices[1] : (subsetStartIndices[length(subsetStartIndices)] + maxlen)]
      }
      if (use_coverage) {
        coverage_list[[sequence_list_index]] <- cov_vector[subsetStartIndices[1] : (subsetStartIndices[length(subsetStartIndices)] + maxlen)]
      }
      start_index_list[[sequence_list_index]] <- subsetStartIndices
      if (additional_labels) {
        if (added_label_by_header) {
          label_list[[sequence_list_index]] <- as.character(cut(subsetStartIndices, breaks = c(startNewEntry, length(nucSeq)),
                                                                labels = header_vector, include.lowest = TRUE, right = FALSE))
        } else { 
          label_list[[sequence_list_index]] <- basename(fasta.files[file_index])
        }
      }
      
      nuc_dist_list[[sequence_list_index]] <- nuc_dist
      sequence_list_index <<- sequence_list_index + 1
      num_new_samples <- end_index - start_index + 1   
      num_samples <- num_samples + num_new_samples 
      start_index <<- end_index + 1  
    }
    
    if (!wavenet_format) {
      
      array_list <- purrr::map(1:length(sequence_list),
                               ~sequenceToArray(sequence_list[[.x]], target_middle = target_middle, nuc_dist = nuc_dist_list[[.x]],
                                                maxlen = maxlen, vocabulary = vocabulary, ambiguous_nuc = ambiguous_nuc,
                                                startInd =  start_index_list[[.x]], use_quality = use_quality_score, 
                                                quality_vector = quality_list[[.x]], cnn_format = cnn_format, target_len = target_len,
                                                cov_vector = coverage_list[[.x]], use_coverage = use_coverage, max_cov = max_cov,
                                                n_gram = n_gram, n_gram_stride = n_gram_stride)
      )
      
      if (!is.list(array_list[[1]][[2]])) {
        if (!target_middle) {
          x <- array_list[[1]][[1]]
          y <- array_list[[1]][[2]]
          if (length(array_list) > 1) {
            for (i in 2:length(array_list)) {
              x <- abind::abind(x, array_list[[i]][[1]], along = 1)
              y <- rbind(y, array_list[[i]][[2]])
            }
          }
          
          # coerce y type to matrix
          if (dim(x)[1] == 1) {
            if (is.null(n_gram)) {
              dim(y) <-  c(1, length(vocabulary))
            } else {
              dim(y) <-  c(1, length(vocabulary)^n_gram)
            }
          }
        } else {
          x_1 <- array_list[[1]][[1]][[1]]
          x_2 <- array_list[[1]][[1]][[2]]
          y <- array_list[[1]][[2]]
          if (length(array_list) > 1) {
            for (i in 2:length(array_list)) {
              x_1 <- abind::abind(x_1, array_list[[i]][[1]][[1]], along = 1)
              x_2 <- abind::abind(x_2, array_list[[i]][[1]][[2]], along = 1)
              y <- rbind(y, array_list[[i]][[2]])
            }
          }
          x <- list(x_1, x_2)
          
          # coerce y type to matrix
          if (dim(x_1)[1] == 1) {
            if (is.null(n_gram)) {
              dim(y) <-  c(1, length(vocabulary))
            } else {
              dim(y) <-  c(1, length(vocabulary)^n_gram)
            }
          }
        }
      } else {
        if (!target_middle) {
          x <- array_list[[1]][[1]]
          y <- array_list[[1]][[2]]
          if (length(array_list) > 1) {
            for (i in 2:length(array_list)) {
              x <- abind::abind(x, array_list[[i]][[1]], along = 1)
              for (j in 1:length(y)) {
                y[[j]] <- rbind(y[[j]], array_list[[i]][[2]][[j]] )
              }
            }
          }
          
          # coerce y type to matrix
          if (dim(x)[1] == 1) {
            for (i in 1:length(y)) {
              if (is.null(n_gram)) {
                dim(y[[i]]) <-  c(1, length(vocabulary))
              } else {
                dim(y[[i]]) <-  c(1, length(vocabulary)^n_gram)
              }
            }
          }
        } else {
          x_1 <- array_list[[1]][[1]][[1]]
          x_2 <- array_list[[1]][[1]][[2]]
          y <- array_list[[1]][[2]]
          if (length(array_list) > 1) {
            for (i in 2:length(array_list)) {
              x_1 <- abind::abind(x_1, array_list[[i]][[1]][[1]], along = 1)
              x_2 <- abind::abind(x_2, array_list[[i]][[1]][[2]], along = 1)
              for (j in 1:length(y)) {
                y[[j]] <- rbind(y[[j]], array_list[[i]][[2]][[j]] )
              }
            }
          }
          x <- list(x_1, x_2)
          
          # coerce y type to matrix
          if (dim(x_1)[1] == 1) {
            for (i in 1:length(y)) {
              if (is.null(n_gram)) {
                dim(y[[i]]) <-  c(1, length(vocabulary))
              } else {
                dim(y[[i]]) <-  c(1, length(vocabulary)^n_gram)
              }
            }
          }
        }
      }  
      
      # wavenet format  
    } else {
      
      if (target_len > 1) {
        stop("target_len must be 1 when using wavenet_format")
      }
      
      # one hot encode strings collected in sequence_list and connect arrays
      array_list <- purrr::map(1:length(sequence_list), 
                               ~sequenceToArray(sequence_list[[.x]], ambiguous_nuc = ambiguous_nuc, 
                                                maxlen = maxlen, vocabulary = vocabulary, nuc_dist = nuc_dist_list[[.x]],
                                                startInd =  start_index_list[[.x]], wavenet_format = TRUE, use_quality = use_quality_score, 
                                                quality_vector = quality_list[[.x]], cnn_format = FALSE, n_gram = n_gram,
                                                cov_vector = coverage_list[[.x]], use_coverage = use_coverage, max_cov = max_cov)
      )
      
      x <- array_list[[1]][[1]]
      y <- array_list[[1]][[2]]
      if (length(array_list) > 1) {
        for (i in 2:length(array_list)) {
          x <- abind::abind(x, array_list[[i]][[1]], along = 1)
          y <- abind::abind(y, array_list[[i]][[2]], along = 1)
        }
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
                                                batch.size = batch.size, start_index_list = start_index_list)
        if (add_input_as_seq[i]) {
          label_tensor_list[[i]] <- sequenceToArrayLabel(as.vector(t(label_tensor_list[[i]])), nuc_dist = NULL,
                                                         maxlen = ncol(label_tensor_list[[i]]), vocabulary = vocabulary, ambiguous_nuc = ambiguous_nuc,
                                                         startInd =  1 + ncol(label_tensor_list[[i]]) * (0:(nrow(label_tensor_list[[i]]) - 1)), use_quality = FALSE, 
                                                         quality_vector = list())
        }
      }
    }
    
    
    # empty lists for next batch 
    start_index_list <<- vector("list")
    sequence_list <<- vector("list")
    nuc_dist_list <<- vector("list")
    quality_list <<- vector("list") 
    coverage_list <<- vector("list") 
    sequence_list_index <<- 1
    num_samples <<- 0
    if (additional_labels) {
      label_list <<- vector("list")
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
    return(list(X = x, Y = y))
  }
}

#' Custom generator for fasta files and label targets
#' 
#' @description
#' \code{fastaLabelGenerator} Iterates over folder containing .fasta files and produces one-hot-encoding of predictor sequences 
#' and target variables. Targets will be read from fasta headers.
#' 
#' @inheritParams fastaFileGenerator
#' @param format File format, either fasta or fastq.
#' @param batch.size Number of batches.    
#' @param maxlen Length of predictor sequence.  
#' @param max_iter Stop after max_iter number of iterations failed to produce a new batch. 
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as 0-vector.
#' @param randomFiles Logical, whether to go through files randomly or sequential. 
#' @param step How often to take a sample.
#' @param showWarnings Logical, give warning if character outside vocabulary appears.   
#' @param seed Sets seed for set.seed function, for reproducible results when using \code{randomFiles} or \code{shuffleFastaEntries}  
#' @param shuffleFastaEntries Logical, shuffle fasta entries.
#' @param verbose Whether to show message. 
#' @param numberOfFiles Use only specified number of files, ignored if greater than number of files in corpus.dir. 
#' @param fileLog Write name of files to csv file if path is specified.
#' @param labelVocabulary Character vector of possible targets. Targets outside \code{labelVocabulary} will get discarded.
#' @param reverseComplements Logical, half of batch contains sequences and other its reverse complements. Reverse complement 
#' is given by reversed order of sequence and switching A/T and C/G. \code{batch.size} argument has to be even, otherwise 1 will be added
#' to \code{batch.size}
#' @param ambiguous_nuc How to handle nucleotides outside vocabulary, either "zero", "discard", "empirical" or "equal". If "zero", input gets encoded as zero vector; 
#' if "equal" input is 1/length(vocabulary) x length(vocabulary). If "discard" samples containing nucleotides outside vocabulary get discarded. 
#' If "empirical" use nucleotide distribution of current file.   
#' @param proportion_per_file Numerical value between 0 and 1. Proportion of possible samples to take from one file. Takes samples from random subsequence.  
#' @param read_data If true the first element of output is a list of length 2, each containing one part of paired read. Maxlen should be 2*length of one read.
#' @param use_quality_score Whether to use fastq qualitiy scores. If TRUE input is not one-hot-encoding but corresponds to probabilities.
#' For example (0.97, 0.01, 0.01, 0.01) instead of (1, 0, 0, 0).   
#' @param padding Whether to pad sequences too short for one sample with zeros. 
#' @param added_label_path If not NULL, get output from csv file. The file should have one column named "file" and one column for every label. 
#' @param target_from_csv Path to csv file with target mapping. One column should be called "file" and other entries in row are the targets. 
#' @param skip_amb_nuc Threshold of ambiguous nucleotides to accept in fasta entry. Complete entry will get discarded otherwise.  
#' @param target_split If target gets read from csv file, list of names to devide target tensor into list of tensors.
#' Example: if csv file has header names "file", "label_1", "label_2", "label_3" und target_split = list(c("label_1", "label_2"), "label_3"),
#' this will devide target matrix to list of length 2, where the first element contains columns named "label_1" and "label_2" and the 
#' second entry contains the column named "label_3".
#' @return A list of length 2. First element is a 3-dimensional tensor with dimensions (batch.size, maxlen, length(vocabulary)), encoding 
#' the predictor sequences. Second element is a matrix with dimensions (batch.size, length(labelVocabulary)), encoding the targets. 
#' @import data.table 
#' @export
fastaLabelGenerator <- function(corpus.dir,
                                format = "fasta",
                                batch.size = 256,
                                maxlen = 250,
                                max_iter = 10000,
                                vocabulary = c("a", "c", "g", "t"),
                                verbose = FALSE,
                                randomFiles = FALSE,
                                step = 1, 
                                showWarnings = FALSE,
                                seed = 1234,
                                shuffleFastaEntries = FALSE,
                                numberOfFiles = NULL,
                                fileLog = NULL,
                                labelVocabulary = c("x", "y", "z"),
                                reverseComplements = TRUE,
                                ambiguous_nuc = "zero",
                                proportion_per_file = NULL,
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
                                n_gram = NULL) {
  
  n_gram <- NULL
  if (is.null(use_coverage)) {
    use_coverage <- FALSE
    max_cov <- NULL
  } else {
    max_cov <- use_coverage
    use_coverage <- TRUE
  }
  if (!is.null(concat_seq) && (!all(stringr::str_split(concat_seq,"")[[1]] %in% vocabulary))) {
    stop("Characters of separating sequence should be in vocabulary")
  }
  discard_amb_nuc <- ifelse(ambiguous_nuc == "discard", TRUE, FALSE)
  vocabulary <- stringr::str_to_lower(vocabulary)
  labelVocabulary <- stringr::str_to_lower(labelVocabulary)
  start_index_list <- vector("list")
  file_index <- 1
  num_samples <- 0
  start_index <- 1
  iter <- 1
  concat <- !is.null(concat_seq)
  additional_labels <- !is.null(added_label_path) 
  seq_vector <- NULL
  
  # adjust maxlen for n_gram
  if (!is.null(n_gram)) {
    maxlen <- maxlen + n_gram - 1
  }
  
  for (i in letters) {
    if (!(i %in% stringr::str_to_lower(vocabulary))) {
      amb_nuc_token <- i
      break
    }
  }
  tokenizer_pred <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), c(vocabulary, amb_nuc_token))
  tokenizer_target <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = FALSE, lower = TRUE, filters = "\t\n"),
                                                labelVocabulary) 
  
  fasta.files <-  list_fasta_files(corpus.dir = corpus.dir,
                                   format = format,
                                   file_filter = file_filter)
  num_files <- length(fasta.files)
  
  if (sample_by_file_size) {
    randomFiles <- FALSE
    file_prob <- file.info(fasta.files)$size/sum(file.info(fasta.files)$size)
  }
  
  set.seed(seed)
  if (randomFiles) fasta.files <- sample(fasta.files, replace = FALSE)
  
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
    }
    
    # if (!added_label_by_header_target) {
    #   fasta.file$Header <- rep(basename(fasta.files[file_index]), nrow(fasta.file))
    # }
    col_name <- ifelse(added_label_by_header_target, "header", "file")
    # header_vector <- fasta.file$Header
    labelVocabulary <- names(output_label_csv)
    labelVocabulary <- labelVocabulary[labelVocabulary != "header" & labelVocabulary != "file"]
    if (!is.null(target_split)) {
      check_header_names(target_split = target_split, labelVocabulary = labelVocabulary) 
    }
  }
  
  # regular expression for chars outside vocabulary
  pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
  
  while (length(seq_vector) == 0) {
    
    # pre-load the first file
    fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                   shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                   reverseComplements = reverseComplements, fasta.files = fasta.files,
                                   labelVocabulary = labelVocabulary, filter_header = TRUE, target_from_csv = target_from_csv)
    
    if (concat) {
      cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
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
                                       shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                       reverseComplements = reverseComplements, fasta.files = fasta.files,
                                       labelVocabulary = labelVocabulary, filter_header = TRUE, target_from_csv = target_from_csv)
        
        if (concat) {
          cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
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
                                       shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                       reverseComplements = reverseComplements, fasta.files = fasta.files,
                                       labelVocabulary = labelVocabulary, filter_header = TRUE, target_from_csv = target_from_csv)
        
        if (concat) {
          cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
          fasta.file <- data.frame(Header = basename(fasta.files[file_index]), Sequence = paste(fasta.file$Sequence, collapse = concat_seq),
                                   stringsAsFactors = FALSE)
        }
      }
    }
    
    if (use_coverage) {
      cov_vector <- get_coverage(fasta.file)
    }
    
    # take random subset
    if (!is.null(proportion_per_file)) {
      if (!read_data) {
        fasta_width <- nchar(fasta.file$Sequence)
        sample_range <- floor(fasta_width - (proportion_per_file * fasta_width))
        start <- mapply(sample_range, FUN = sample, size = 1)
        perc_length <- floor(fasta_width * proportion_per_file)
        stop <- start + perc_length
        seq_vector <- mapply(fasta.file$Sequence, FUN = substr, start = start, stop = stop)
        if (use_quality_score) {
          quality_scores <- mapply(fasta.file$Quality, FUN = substr, start = start, stop = stop)
        }
      } else {
        sample_index <- sample(nrow(fasta.file), ceiling(proportion_per_file * nrow(fasta.file)))
        fasta.file <- fasta.file[sample_index,]
        seq_vector <- fasta.file$Sequence
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
        if (randomFiles) fasta.files <- sample(fasta.files, replace = FALSE)
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
    start_indices <- getStartInd(seq_vector = seq_vector, length_vector = length_vector, maxlen = maxlen, step = step,
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
  if (!is.null(numberOfFiles) && (numberOfFiles < length(fasta.files))) {
    fasta.files <- fasta.files[1:numberOfFiles]
    num_files <- length(fasta.files)
  }
  
  # log file
  if (!is.null(fileLog)) {
    if (!endsWith(fileLog, ".csv")) fileLog <- paste0(fileLog, ".csv")
    # if (file.exists(fileLog)) {
    #   write.table(x = filePath, file = fileLog, append = TRUE, col.names = FALSE, row.names = FALSE)
    # } else {
    write.table(x = fasta.files[1], file = fileLog, row.names = FALSE, col.names = FALSE)
    # }
  }
  
  # test for chars outside vocabulary
  if (showWarnings) {
    charsOutsideVoc <- stringr::str_detect(stringr::str_to_lower(nucSeq), pattern)  
    if (charsOutsideVoc) warning("file ", filePath, " contains characters outside vocabulary")
  }
  
  if (verbose) message("Initializing ...")
  rngstate <- .GlobalEnv$.Random.seed
  
  function() {
    
    .GlobalEnv$.Random.seed <- rngstate
    on.exit(rngstate <<- .GlobalEnv$.Random.seed)
    iter <- 1
    # loop until enough samples collected
    while(num_samples < batch.size) {  
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
          if (randomFiles) fasta.files <<- sample(fasta.files, replace = FALSE)
          file_index <<- 1
        }
        
        # skip empty files
        while(TRUE) {
          fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                         shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                         reverseComplements = reverseComplements, fasta.files = fasta.files, 
                                         labelVocabulary = labelVocabulary, filter_header = TRUE, target_from_csv = target_from_csv)
          
          if (concat) {
            cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
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
            if (randomFiles) fasta.files <<- sample(fasta.files, replace = FALSE)
            file_index <<- 1
          }
        }
        
        if (use_coverage) {
          cov_vector <<- get_coverage(fasta.file)
        }
        
        # take random subset
        if (!is.null(proportion_per_file)) {
          if (!read_data) {
            fasta_width <- nchar(fasta.file$Sequence)
            sample_range <- floor(fasta_width - (proportion_per_file * fasta_width))
            start <- mapply(sample_range, FUN = sample, size = 1)
            perc_length <- floor(fasta_width * proportion_per_file)
            stop <- start + perc_length
            seq_vector <- mapply(fasta.file$Sequence, FUN = substr, start = start, stop = stop)
            if (use_quality_score) {
              quality_scores <- mapply(fasta.file$Quality, FUN = substr, start = start, stop = stop)
            }
          } else {
            sample_index <- sample(nrow(fasta.file), ceiling(proportion_per_file * nrow(fasta.file)))
            fasta.file <- fasta.file[sample_index,]
            seq_vector <- fasta.file$Sequence
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
        if (!is.null(fileLog)) {
          write.table(x = fasta.files[file_index], file = fileLog, append = TRUE, col.names = FALSE, row.names = FALSE)
        }
        # test for chars outside vocabulary
        if (showWarnings) {
          charsOutsideVoc <- stringr::str_detect(stringr::str_to_lower(nucSeq), pattern)  
          if (charsOutsideVoc) warning("file ", filePath, " contains characters outside vocabulary")
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
              quality_scores <- quality_scores[-short_seq_index]
            }
            if (use_coverage) {
              cov_vector <- cov_vector[-short_seq_index]
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
          start_indices <<- getStartInd(seq_vector = seq_vector, length_vector = length_vector, maxlen = maxlen, step = step,
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
      remainingSamples <- batch.size - num_samples
      end_index <- min(length(start_indices), start_index + remainingSamples  - 1)
      
      subsetStartIndices <- start_indices[start_index:end_index]
      sequence_list[[sequence_list_index]] <- nucSeq[subsetStartIndices[1] : (subsetStartIndices[length(subsetStartIndices)] + maxlen - 1)]
      # collect targets 
      target_list[[sequence_list_index]] <- as.character(cut(subsetStartIndices, breaks = c(startNewEntry, length(nucSeq)),
                                                             labels = label_vector, include.lowest = TRUE, right = FALSE))
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
    array_x_list <- purrr::map(1:length(sequence_list), ~sequenceToArrayLabel(sequence_list[[.x]], ambiguous_nuc = ambiguous_nuc,
                                                                              maxlen = maxlen, vocabulary = vocabulary, nuc_dist = nuc_dist_list[[.x]],
                                                                              startInd =  start_index_list[[.x]], quality_vector = quality_list[[.x]],
                                                                              use_quality = use_quality_score, cov_vector = coverage_list[[.x]],
                                                                              use_coverage = use_coverage, max_cov = max_cov, n_gram = n_gram)
    )
    
    # one hot encode targets
    if (is.null(added_label_path)) {
      target_int <- unlist(keras::texts_to_sequences(tokenizer_target, unlist(target_list))) - 1
      y  <- keras::to_categorical(target_int, num_classes = length(labelVocabulary))
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
                                                batch.size = batch.size, start_index_list = start_index_list)
        if (add_input_as_seq[i]) {
          label_tensor_list[[i]] <- sequenceToArrayLabel(as.vector(t(label_tensor_list[[i]])), nuc_dist = NULL,
                                                         maxlen = ncol(label_tensor_list[[i]]), vocabulary = vocabulary, ambiguous_nuc = ambiguous_nuc,
                                                         startInd =  1 + ncol(label_tensor_list[[i]]) * (0:(nrow(label_tensor_list[[i]]) - 1)), use_quality = FALSE, 
                                                         quality_vector = list())
        }
      }
    }
    
    if (!is.null(target_from_csv)) {
      .datatable.aware = TRUE
      output_label_vector <- unlist(output_label_list) # %>% stringr::str_to_lower()
      target_tensor <- matrix(0, ncol = ncol(output_label_csv) - 1, nrow = batch.size, byrow = TRUE)
      
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
      dim(y) <- c(1, length(labelVocabulary))
    }
    
    # empty sequence_list for next batch 
    start_index_list <<- vector("list")
    sequence_list <<- vector("list")
    target_list <<- vector("list")
    nuc_dist_list <<- vector("list")
    quality_list <<- vector("list") 
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
      colnames(y) <- labelVocabulary
      y <- slice_tensor(tensor = y, target_split = target_split)
    }
    return(list(X = x, Y = y))
  }
}

#' Custom generator for fasta files
#'
#' @description
#' \code{labelByFolderGenerator} Iterates over folder containing .fasta files and produces one-hot-encoding of predictor sequences 
#' and target variables. Files in \code{corpus.dir} should all belong to one class.
#' 
#' @inheritParams fastaFileGenerator
#' @param corpus.dir Input directory where .fasta files are located or path to single file ending with .fasta or .fastq 
#' (as specified in format argument).
#' @param format File format, either fasta or fastq.
#' @param batch.size Number of batches.    
#' @param maxlen Length of predictor sequence.  
#' @param max_iter Stop after max_iter number of iterations failed to produce a new batch. 
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as 0-vector.
#' @param randomFiles Logical, whether to go through files randomly or sequential. 
#' @param step How often to take a sample.
#' @param showWarnings Logical, give warning if character outside vocabulary appears   
#' @param seed Sets seed for set.seed function, for reproducible results when using \code{randomFiles} or \code{shuffleFastaEntries}  
#' @param shuffleFastaEntries Logical, shuffle fasta entries.
#' @param verbose Whether to show message. 
#' @param numberOfFiles Use only specified number of files, ignored if greater than number of files in corpus.dir. 
#' @param fileLog Write name of files to csv file if path is specified.
#' @param reverseComplements Logical, for every new file decide randomly to use original data or its reverse complement.
#' @param reverseComplementEncoding Logical, use both original sequence and reverse.complement as two input sequences.
#' @param numTargets Number of columns of target matrix.  
#' @param onesColumn Which column of target matrix contains ones
#' @param ambiguous_nuc How to handle nucleotides outside vocabulary, either "zero", "discard", "empirical" or "equal". If "zero", input gets encoded as zero vector; 
#' if "equal" input is 1/length(vocabulary) x length(vocabulary). If "discard" samples containing nucleotides outside vocabulary get discarded. 
#' If "empirical" use nucleotide distribution of current file.    
#' @param proportion_per_file Numerical value between 0 and 1. Proportion of possible samples to take from one file. Takes samples from random subsequence. 
#' @param read_data If true the first element of output is a list of length 2, each containing one part of paired read. Maxlen should be 2*length of one read.
#' @param use_quality_score Whether to use fastq qualitiy scores. If TRUE input is not one-hot-encoding but corresponds to probabilities.
#' For example (0.97, 0.01, 0.01, 0.01) instead of (1, 0, 0, 0).   
#' @param padding Whether to pad sequences too short for one sample with zeros. 
#' @return A list of length 2. First element is a 3-dimensional tensor with dimensions (batch.size, maxlen, length(vocabulary)), encoding 
#' the predictor sequences. Second element is a matrix with dimensions (batch.size, numTargets), encoding the targets. If 
#' \code{read_data = TRUE} first element is a list of two 3-dimensional tensor with dimensions (batch.size, maxlen, length(vocabulary)) each.
#' @param skip_amb_nuc Threshold of ambiguous nucleotides to accept in fasta entry. Complete entry will get discarded otherwise.  
#' @param split_seq Split input sequence into two sequences while removing nucleotide in middle. If input is x_1,..., x_(n+1), input gets split into 
#' input_1 = x_1,..., x_m and input_2 = x_(n+1),..., x_(m+2) where m = ceiling((n+1)/2) and n = maxlen. Note that x_(m+1) is not used. Can be used for transfer learning,
#' when switching from language model trained with target in middle to label classification. 
#' @export
labelByFolderGenerator <- function(corpus.dir,
                                   format = "fasta",
                                   batch.size = 256,
                                   maxlen = 250,
                                   max_iter = 10000,
                                   vocabulary = c("a", "c", "g", "t"),
                                   verbose = FALSE,
                                   randomFiles = FALSE,
                                   step = 1, 
                                   showWarnings = FALSE,
                                   seed = 1234,
                                   shuffleFastaEntries = FALSE,
                                   numberOfFiles = NULL,
                                   fileLog = NULL,
                                   reverseComplements = TRUE,
                                   reverseComplementEncoding = FALSE,
                                   numTargets,
                                   onesColumn,
                                   ambiguous_nuc = "zero",
                                   proportion_per_file = NULL,
                                   read_data = FALSE,
                                   use_quality_score = FALSE,
                                   padding = TRUE,
                                   added_label_path = NULL,
                                   add_input_as_seq = NULL,
                                   skip_amb_nuc = NULL,
                                   max_samples = NULL,
                                   split_seq = FALSE,
                                   concat_seq = NULL,
                                   use_coverage = NULL,
                                   proportion_entries = NULL,
                                   sample_by_file_size = FALSE,
                                   n_gram = NULL) {
  
  n_gram <- NULL
  if (is.null(use_coverage)) {
    use_coverage <- FALSE
    max_cov <- NULL
  } else {
    max_cov <- use_coverage
    use_coverage <- TRUE
  }
  if (!is.null(concat_seq) && (!all(stringr::str_split(concat_seq,"")[[1]] %in% vocabulary))) {
    stop("Characters of separating sequence should be in vocabulary")
  }
  if (reverseComplementEncoding) {
    test_len <- length(vocabulary) != 4 
    if (test_len || all(sort(stringr::str_to_lower(vocabulary)) != c("a", "c", "g", "t"))) {
      stop("reverseComplementEncoding only implemented for A,C,G,T vocabulary yet")
    }
  } 
  stopifnot(!split_seq | !read_data)
  stopifnot(!split_seq | !reverseComplementEncoding)
  stopifnot(!(read_data & padding))
  stopifnot(onesColumn <= numTargets)
  if (read_data & !is.null(skip_amb_nuc)) {
    stop("Using read data and skipping files at the same time not implemented yet")
  }
  additional_labels <- ifelse(is.null(added_label_path), FALSE, TRUE)
  
  # adjust maxlen for n_gram
  if (!is.null(n_gram)) {
    maxlen <- maxlen + n_gram - 1
  }
  
  if (split_seq) {
    maxlen <- maxlen + 1
    len_input_1 <- floor(maxlen/2)
  }  
  
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
  
  if (is.list(corpus.dir)) {
    fasta.files <- list()
    for (i in 1:length(corpus.dir)) {
      fasta.files[[i]] <- list.files(
        path = xfun::normalize_path(corpus.dir[[i]]),
        pattern = paste0("\\.", format, "$"),
        full.names = TRUE)
    }
    fasta.files <- unlist(fasta.files)
    num_files <- length(fasta.files)
  } else {
    
    # single file
    if (endsWith(corpus.dir, paste0(".", format))) {
      num_files <- 1 
      fasta.files <- corpus.dir   
    } else {
      
      fasta.files <- list.files(
        path = xfun::normalize_path(corpus.dir),
        pattern = paste0("\\.", format, "$"),
        full.names = TRUE)
      num_files <- length(fasta.files)
    }
  }
  
  if (sample_by_file_size) {
    randomFiles <- FALSE
    file_prob <- file.info(fasta.files)$size/sum(file.info(fasta.files)$size)
  }
  
  set.seed(seed)
  if (randomFiles) fasta.files <- sample(fasta.files, replace = FALSE)
  
  if (read_data) {
    contains_R1 <-  stringr::str_detect(fasta.files, "R1")
    fasta.files <- fasta.files[contains_R1]
  }
  
  # regular expression for chars outside vocabulary
  pattern <- paste0("[^", paste0(vocabulary, collapse = ""), "]")
  
  while (length(seq_vector) == 0) {
    
    fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                   shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                   reverseComplements = reverseComplements, fasta.files = fasta.files)
    
    if (concat) {
      cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
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
                                       shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                       reverseComplements = reverseComplements, fasta.files = fasta.files)
        
        if (concat) {
          cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
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
                                       shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                       reverseComplements = reverseComplements, fasta.files = fasta.files)
        
        if (concat) {
          cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
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
    if (!is.null(proportion_per_file)) {
      if (!read_data) {
        fasta_width <- nchar(fasta.file$Sequence)
        sample_range <- floor(fasta_width - (proportion_per_file * fasta_width))
        start <- mapply(sample_range, FUN = sample, size = 1)
        perc_length <- floor(fasta_width * proportion_per_file)
        stop <- start + perc_length
        seq_vector <- mapply(fasta.file$Sequence, FUN = substr, start = start, stop = stop)
        if (use_quality_score) {
          quality_scores <- mapply(fasta.file$Quality, FUN = substr, start = start, stop = stop)
        }
      } else {
        sample_index <- sample(nrow(fasta.file), ceiling(proportion_per_file * nrow(fasta.file)))
        fasta.file <- fasta.file[sample_index,]
        seq_vector <- fasta.file$Sequence
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
        if (randomFiles) fasta.files <- sample(fasta.files, replace = FALSE)
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
    start_indices <- getStartInd(seq_vector = seq_vector, length_vector = length_vector, maxlen = maxlen, step = step,
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
  if (!is.null(numberOfFiles) && (numberOfFiles < length(fasta.files))) {
    fasta.files <- fasta.files[1:numberOfFiles]
    num_files <- length(fasta.files)
  }
  
  # log file
  if (!is.null(fileLog)) {
    if (!endsWith(fileLog, ".csv")) fileLog <- paste0(fileLog, ".csv")
    append <- file.exists(fileLog)
    write.table(x = fasta.files[file_index], file = fileLog, col.names = FALSE, row.names = FALSE, append = append)
  }
  
  # test for chars outside vocabulary
  if (showWarnings) {
    charsOutsideVoc <- stringr::str_detect(stringr::str_to_lower(nucSeq), pattern)  
    if (charsOutsideVoc) warning("file ", filePath, " contains characters outside vocabulary")
  }
  
  if (verbose) message("Initializing ...")
  
  rngstate <- .GlobalEnv$.Random.seed
  
  function() {
    
    .GlobalEnv$.Random.seed <- rngstate
    on.exit(rngstate <<- .GlobalEnv$.Random.seed)
    
    iter <- 1
    
    # loop until enough samples collected
    while(num_samples < batch.size) {  
      
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
          if (randomFiles) fasta.files <<- sample(fasta.files, replace = FALSE)
          file_index <<- 1
        }
        
        # skip empty files
        while(TRUE) {
          fasta.file <- read_fasta_fastq(format = format, skip_amb_nuc =  skip_amb_nuc, file_index = file_index, pattern = pattern,
                                         shuffleFastaEntries = shuffleFastaEntries, proportion_entries = proportion_entries,
                                         reverseComplements = reverseComplements, fasta.files = fasta.files)
          
          if (concat) {
            cov_vector <- get_coverage_concat(fasta.file = fasta.file, concat_seq = concat_seq)
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
            if (randomFiles) fasta.files <<- sample(fasta.files, replace = FALSE)
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
        if (!is.null(proportion_per_file)) {
          if (!read_data) {
            fasta_width <- nchar(fasta.file$Sequence)
            sample_range <- floor(fasta_width - (proportion_per_file * fasta_width))
            start <- mapply(sample_range, FUN = sample, size = 1)
            perc_length <- floor(fasta_width * proportion_per_file)
            stop <- start + perc_length
            seq_vector <- mapply(fasta.file$Sequence, FUN = substr, start = start, stop = stop)
            if (use_quality_score) {
              quality_scores <- mapply(fasta.file$Quality, FUN = substr, start = start, stop = stop)
            }
          } else {
            sample_index <- sample(nrow(fasta.file), ceiling(proportion_per_file * nrow(fasta.file)))
            fasta.file <- fasta.file[sample_index,]
            seq_vector <- fasta.file$Sequence
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
        if (!is.null(fileLog)) {
          write.table(x = fasta.files[file_index], file = fileLog, append = TRUE, col.names = FALSE, row.names = FALSE)
        }
        # test for chars outside vocabulary
        if (showWarnings) {
          charsOutsideVoc <- stringr::str_detect(stringr::str_to_lower(nucSeq), pattern)  
          if (charsOutsideVoc) warning("file ", filePath, " contains characters outside vocabulary")
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
              cov_vector <- cov_vector[-short_seq_index]
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
          start_indices <<- getStartInd(seq_vector = seq_vector, length_vector = length_vector, maxlen = maxlen, step = step,
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
      remainingSamples <- batch.size - num_samples
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
    array_x_list <- purrr::map(1:length(sequence_list), ~sequenceToArrayLabel(sequence = sequence_list[[.x]], ambiguous_nuc = ambiguous_nuc,
                                                                              maxlen = maxlen, vocabulary = vocabulary, nuc_dist = nuc_dist_list[[.x]],
                                                                              startInd =  start_index_list[[.x]], use_quality = use_quality_score,
                                                                              quality_vector = quality_list[[.x]], cov_vector = coverage_list[[.x]],
                                                                              max_cov = max_cov, use_coverage = use_coverage, n_gram = n_gram)
    )
    
    x <- array_x_list[[1]]
    
    if (length(array_x_list) > 1) {
      for (i in 2:length(array_x_list)) {
        x <- abind::abind(x, array_x_list[[i]], along = 1)
      }
    }
    
    # one hot encode targets
    y  <- matrix(0, ncol = numTargets, nrow = dim(x)[1])
    y[ , onesColumn] <- 1
    
    # coerce y type to matrix
    if (dim(x)[1] == 1) {
      dim(y) <- c(1, numTargets)
    }
    
    if (reverseComplementEncoding){
      x_1 <- x
      x_2 <- array(x_1[ , (dim(x)[2]):1, 4:1], dim = c(dim(x)[1],dim(x)[2],dim(x)[3])) 
      x <- list(x_1, x_2)
    }
    
    if (split_seq) {
      x_1 <- array(x[ , 1:len_input_1, ], dim = c(batch.size, len_input_1, length(vocabulary)))
      x_2 <- array(x[ , maxlen : (len_input_1 + 2), ], dim = c(batch.size, (maxlen - len_input_1 - 1), length(vocabulary))) 
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
                                                batch.size = batch.size, start_index_list = start_index_list)
        if (add_input_as_seq[i]) {
          label_tensor_list[[i]] <- sequenceToArrayLabel(as.vector(t(label_tensor_list[[i]])), nuc_dist = NULL,
                                                         maxlen = ncol(label_tensor_list[[i]]), vocabulary = vocabulary, ambiguous_nuc = ambiguous_nuc,
                                                         startInd =  1 + ncol(label_tensor_list[[i]]) * (0:(nrow(label_tensor_list[[i]]) - 1)), use_quality = FALSE, 
                                                         quality_vector = list())
        }
      }
    }
    
    # empty lists for next batch 
    start_index_list <<- vector("list")
    sequence_list <<- vector("list")
    target_list <<- vector("list")
    quality_list <<- vector("list") 
    nuc_dist_list <<- vector("list")
    coverage_list <<- vector("list") 
    sequence_list_index <<- 1
    num_samples <<- 0
    if (reverseComplementEncoding){
      return(list(X = x, Y = y))
    }
    if (split_seq) {
      return(list(X = list(x_1, x_2), Y = y)) 
      if (additional_labels) {
        label_tensor_list[[length(label_tensor_list) + 1]] <- x_1
        label_tensor_list[[length(label_tensor_list) + 1]] <- x_2 
        x <- label_tensor_list
      } 
      return(list(X = x, Y = y)) 
      
    } else {
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
      return(list(X = x, Y = y))
    }
  }
}

#' Initializes generators defined by labelByFolderGenerator function
#'
#' \code{initializeGenerators} Initializes generators defined by \code{\link{{labelByFolderGenerator}} function. Targets get one-hot-encoded in order of directories.
#' Number of classes is given by length of directories.
#' 
#' @inheritParams fastaFileGenerator
#' @param directories Vector of paths to folder containing fasta files. Files in one folder should belong to one class. 
#' @param format File format.
#' @param batch.size Number of batches, will get rounded to be multiple of number of targets if necessary.
#' @param maxlen Length of predictor sequence.  
#' @param max_iter Stop after max_iter number of iterations failed to produce a new batch. 
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as 0-vector.
#' @param randomFiles Logical, whether to go through files randomly or sequential. 
#' @param step How often to take a sample.
#' @param showWarnings Logical, give warning if character outside vocabulary appears.   
#' @param seed Sets seed for set.seed function, for reproducible results when using \code{randomFiles} or \code{shuffleFastaEntries}  
#' @param shuffleFastaEntries Logical, shuffle fasta entries.
#' @param verbose Whether to show message. 
#' @param numberOfFiles Use only specified number of files, ignored if greater than number of files in \code{directories}. 
#' @param fileLog Write name of files to csv file if path is specified.
#' @param reverseComplements Logical, half of batch contains sequences and other its reverse complements. Reverse complement 
#' is given by reversed order of sequence and switching A/T and C/G. \code{batch.size} argument has to be even, otherwise 1 will be added
#' to \code{batch.size}
#' @param reverseComplementEncoding Logical, use both original sequence and reverse.complement as two input sequences.
#' @param val Logical, call initialized generarator "genY" or "genValY" where Y is an integer between 1 and length of directories.
#' @param ambiguous_nuc How to handle nucleotides outside vocabulary, either "zero", "discard" or "equal". If "zero", input gets encoded as zero vector; 
#' if "equal" input is 1/length(vocabulary) x length(vocabulary). If "discard" samples containing nucleotides outside vocabulary get discarded.     
#' @param proportion_per_file Numerical value between 0 and 1. Proportion of possible samples to take from one file. Takes samples from random subsequence.  
#' @param target_middle Split input sequence into two sequences while removing nucleotide in middle. If input is x_1,..., x_(n+1), input gets split into 
#' input_1 = x_1,..., x_m and input_2 = x_(n+1),..., x_(m+2) where m = ceiling((n+1)/2) and n = maxlen. Note that x_(m+1) is not used.     
#' @param read_data If true the first element of output is a list of length 2, each containing one part of paired read.
#' @param use_quality_score Whether to use fastq qualitiy scores. If TRUE input is not one-hot-encoding but corresponds to probabilities.
#' For example (0.97, 0.01, 0.01, 0.01) instead of (1, 0, 0, 0).   
#' @param padding Whether to pad sequences too short for one sample with zeros. 
#' @export
initializeGenerators <- function(directories,
                                 format = "fasta",
                                 batch.size = 256,
                                 maxlen = 250,
                                 max_iter = 10000,
                                 vocabulary = c("a", "c", "g", "t"),
                                 verbose = FALSE,
                                 randomFiles = FALSE,
                                 step = 1, 
                                 showWarnings = FALSE,
                                 seed = 1234,
                                 shuffleFastaEntries = FALSE,
                                 numberOfFiles = NULL,
                                 fileLog = NULL,
                                 reverseComplements = FALSE, 
                                 reverseComplementEncoding = FALSE,
                                 val = FALSE,
                                 ambiguous_nuc = "zero",
                                 proportion_per_file = NULL,
                                 target_middle = FALSE,
                                 read_data = FALSE,
                                 use_quality_score = FALSE,
                                 padding = TRUE,
                                 added_label_path = NULL,
                                 add_input_as_seq = NULL,
                                 skip_amb_nuc = NULL,
                                 max_samples = NULL,
                                 split_seq = FALSE,
                                 concat_seq = NULL,
                                 use_coverage = NULL,
                                 set_learning = NULL,
                                 proportion_entries = NULL,
                                 sample_by_file_size = FALSE,
                                 n_gram = NULL) {
  
  # adjust batch.size
  if (is.null(set_learning)) {
    if (length(batch.size) == 1) {
      if (batch.size %% length(directories) != 0) {
        batchType <- ifelse(val, "validation", "training")
        batch.size <- ceiling(batch.size/length(directories)) * length(directories)
        if (!reverseComplements) {
          message(paste("Batch size needs to be multiple of number of targets. Setting", batchType, "batch size to", batch.size))
        } 
      }
    } 
  }
  
  numTargets <- length(directories)
  
  if (length(batch.size) == 1) {
    batch.size <- rep(batch.size/numTargets, numTargets)
  }
  
  argg <- c(as.list(environment()))
  # variables with just one entry  
  argg["directories"] <- NULL
  argg["val"] <- NULL
  argg["vocabulary"] <- NULL
  argg["numTargets"] <- NULL
  argg["verbose"] <- NULL
  argg["maxlen"] <- NULL
  argg["showWarnings"] <- NULL
  argg["max_iter"] <- NULL
  argg["read_data"] <- NULL
  argg["use_quality_score"] <- NULL
  argg["added_label_path"] <- NULL
  argg["add_input_as_seq"] <- NULL
  argg["skip_amb_nuc"] <- NULL
  argg["max_samples"] <- NULL
  argg["split_seq"] <- NULL
  argg["concat_seq"] <- NULL
  argg["reverseComplementEncoding"] <- NULL
  argg["use_coverage"] <- NULL
  argg["set_learning"] <- NULL
  argg["proportion_entries"] <- NULL
  argg["sample_by_file_size"] <- NULL
  argg["n_gram"] <- NULL
  
  for (i in 1:length(argg)) {
    if (length(argg[[i]]) == 1) {
      assign(names(argg)[i], rep(argg[[i]], numTargets))
    } 
    if ((length(argg[[i]]) != 1) & (length(argg[[i]]) != numTargets) & !(is.null(argg[[i]]))) {
      stop_message <- paste("Incorrect argument length,", names(argg[i]), "argument vector must have length 1 or", numTargets)
      stop(stop_message)
    }
  }
  
  if (!val) {
    # create generator for every folder
    for (i in 1:length(directories)) {
      numberedGen <- paste0("gen", as.character(i))
      genAsText <- paste(numberedGen, "<<- labelByFolderGenerator(corpus.dir = directories[[i]],
                                       format = format[i],
                                       batch.size = batch.size[i],
                                       maxlen = maxlen,
                                       max_iter = max_iter,
                                       vocabulary = vocabulary,
                                       verbose = verbose,
                                       randomFiles = randomFiles[i],
                                       step = step[i], 
                                       showWarnings = showWarnings,
                                       seed = seed[i],
                                       shuffleFastaEntries = shuffleFastaEntries[i],
                                       numberOfFiles = numberOfFiles[i],
                                       fileLog = fileLog[i],
                                       reverseComplements = reverseComplements[i],
                                       reverseComplementEncoding = reverseComplementEncoding,
                                       numTargets = numTargets,
                                       onesColumn = i,
                                       ambiguous_nuc = ambiguous_nuc[i],
                                       proportion_per_file = proportion_per_file[i],
                                       read_data = read_data,
                                       use_quality_score = use_quality_score,
                                       padding = padding[i],
                                       added_label_path = added_label_path,
                                       add_input_as_seq = add_input_as_seq,
                                       skip_amb_nuc = skip_amb_nuc,
                                       max_samples = max_samples,
                                       split_seq = split_seq,
                                       concat_seq = concat_seq,
                                       use_coverage = use_coverage,
                                       proportion_entries = proportion_entries,
                                       sample_by_file_size = sample_by_file_size,
                                       n_gram = n_gram
  )"
      )  
      eval(parse(text = genAsText))  
    } 
  } else {
    # create generator for every folder
    for (i in 1:length(directories)) {
      # different names for validation generators
      numberedGenVal <- paste0("genVal", as.character(i))
      genAsTextVal <- paste(numberedGenVal, "<<- labelByFolderGenerator(corpus.dir = directories[[i]],
                                       format = format[i],
                                       batch.size = batch.size[i],
                                       maxlen = maxlen,
                                       max_iter = max_iter,
                                       vocabulary = vocabulary,
                                       verbose = verbose,
                                       randomFiles = randomFiles[i],
                                       step = step[i], 
                                       showWarnings = showWarnings,
                                       seed = seed[i],
                                       shuffleFastaEntries = shuffleFastaEntries[i],
                                       numberOfFiles = numberOfFiles[i],
                                       fileLog = fileLog[i],
                                       reverseComplements = reverseComplements[i],
                                       reverseComplementEncoding = reverseComplementEncoding,
                                       numTargets = numTargets,
                                       onesColumn = i,
                                       ambiguous_nuc = ambiguous_nuc[i],
                                       proportion_per_file = proportion_per_file[i],
                                       read_data = read_data,
                                       use_quality_score = use_quality_score,
                                       padding = padding[i],
                                       added_label_path = added_label_path,
                                       add_input_as_seq = add_input_as_seq,
                                       skip_amb_nuc = skip_amb_nuc,
                                       max_samples = max_samples,
                                       split_seq = split_seq,
                                       concat_seq = concat_seq,
                                       use_coverage = use_coverage,
                                       proportion_entries = proportion_entries,
                                       sample_by_file_size = sample_by_file_size,
                                       n_gram = n_gram
  )"
      )
      eval(parse(text = genAsTextVal))  
    }
  }
}

#' Generator wrapper
#'
#' Iterates over generators created by \code{\link{{initializeGenerators}}
#' 
#' @param val Train or validation generator.
#' @param path Path 
#' @export
labelByFolderGeneratorWrapper <- function(val, new_batch_size = NULL,
                                          samples_per_target = NULL, batch.size = NULL,
                                          path = NULL, voc_len = NULL, maxlen = NULL,
                                          reshape_mode = NULL,  buffer_len = NULL,
                                          concat_maxlen = NULL) {
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
                            batch.size = batch.size,
                            path = path,
                            voc_len = voc_len,
                            maxlen = maxlen,
                            concat_maxlen = concat_maxlen,
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
                            batch.size = batch.size,
                            path = path,
                            voc_len = voc_len,
                            maxlen = maxlen,
                            concat_maxlen = concat_maxlen,
                            buffer_len = buffer_len,
                            reshape_mode = reshape_mode)
        return(l)
      } else {
        return(list(X = xVal, Y = yVal))
      }
    }
  }
}

#' creates random data
#' 
#' @param model A keras model.
#' @export
dummy_gen <- function(model, batch_size) {
  
  num_input_layers <- ifelse(is.list(model$input), length(model$inputs), 1)
  x <- list()
  if (num_input_layers == 1) {
    input_dim <- batch_size
    for (j in 2:length(model$input_shape)) {
      input_dim  <- c(input_dim, model$input_shape[[j]])
    }
    x <- array(rnorm(n = prod(input_dim)), dim = input_dim) 
  } else {
    
    for (i in 1:num_input_layers) {
      input_dim <- batch_size
      input_list <- model$input_shape[[i]]
      for (j in 2:length(input_list)) {
        input_dim  <- c(input_dim, input_list[[j]])
      }
      x[[i]] <- array(rnorm(n = prod(input_dim)), dim = input_dim) 
    }
  }
  
  num_output_layers <- ifelse(is.list(model$output), length(model$outputs), 1)
  y <- list()
  if (num_output_layers == 1) {
    output_dim <- batch_size
    for (j in 2:length(model$output_shape)) {
      output_dim  <- c(output_dim, model$output_shape[[j]])
    }
    y <- array(rnorm(n = prod(output_dim)), dim = output_dim) 
  } else {
    
    
    for (i in 1:num_output_layers) {
      output_dim <- batch_size
      output_list <- model$output_shape[[i]]
      for (j in 2:length(output_list)) {
        output_dim  <- c(output_dim, output_list[[j]])
      }
      y[[i]] <- array(rnorm(n = prod(output_dim)), dim = output_dim) 
    }
  }
  
  function() {
    return(list(x,y))  
  }
}  

#' Read batches from rds arrays.
#' 
#' @param rds_folder Path to input files.
#' @export
gen_rds <- function(rds_folder, batch_size, fileLog = NULL, 
                    max_samples = NULL,
                    proportion_per_file = NULL,
                    target_len = NULL,
                    n_gram = NULL, n_gram_stride = 1) {
  
  is_lm <- !is.null(target_len)
  
  if (!is.null(n_gram) & is_lm && (target_len < n_gram)) {
    stop("target_len needs to be at least as big as n_gram.")
  }
  
  # original_target_len <- target_len
  #if (!is.null(n_gram)) {
  #   target_len <- target_len + n_gram - 1
  #}    
  
  if (endsWith(rds_folder, ".rds")) {
    rds_files <- rds_folder
  } else {
    rds_files <- list.files(unlist(rds_folder), pattern = ".rds", full.names = TRUE)
  }
  
  file_index <- 1
  num_files <- length(rds_files)
  rds_file <- readRDS(rds_files[file_index])
  x_complete <- rds_file[[1]]
  if (!is_lm) y_complete <- rds_file[[2]]
  # if (!is.null(n_gram)) {
  #   x_complete <- n_gram_of_3d_tensor(tensor_3d = x_complete, n = n_gram)
  # }
  
  x_dim_start <- dim(x_complete)
  if (!is_lm) { 
    y_dim_start <- dim(y_complete) 
    if (x_dim_start[1] != y_dim_start[1]) {
      stop("Different number of samples for input and target")
    }
  }  
  sample_index <- 1:x_dim_start[1]  
  
  if (!is.null(proportion_per_file)) {
    sample_index <- sample(sample_index, min(length(sample_index), length(sample_index) * proportion_per_file))
  }
  
  if (!is.null(max_samples)) {
    sample_index <- sample(sample_index, min(length(sample_index), max_samples))
  }
  
  if (!is.null(fileLog)) {
    write.table(x = rds_files[1], file = fileLog, row.names = FALSE, col.names = FALSE)
  }
  
  x_dim <- x_dim_start
  if (!is_lm) y_dim <- y_dim_start
  
  function() {
    
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
          file_index <<- file_index + 1
          if (file_index > num_files) {
            file_index <<- 1
            rds_files <<- sample(rds_files)
          } 
          rds_file <<- readRDS(rds_files[file_index])
          x_complete <<- rds_file[[1]]
          # if (!is.null(n_gram)) {
          #   x_complete <<- n_gram_of_3d_tensor(tensor_3d = x_complete, n = n_gram)
          # }
          
          x_dim <<- dim(x_complete)
          
          if (!is_lm) {
            y_complete <<- rds_file[[2]]
            y_dim <<- dim(y_complete)
          }
          
          if (!is_lm && (x_dim[1] != y_dim[1])) {
            stop("Different number of samples for input and target")
          }
          
          sample_index <<- 1:x_dim[1] 
          if (!is.null(proportion_per_file)) {
            if (length(sample_index) > 1) {
              sample_index <<- sample(sample_index, min(length(sample_index), max(1, floor(length(sample_index) * proportion_per_file))))
            }
          }
          if (!is.null(max_samples)) {
            if (length(sample_index) > 1) {
              sample_index <<- sample(sample_index, min(length(sample_index), max_samples))
            }
          }  
        }
        
        # log file
        if (!is.null(fileLog)) {
          write.table(x = rds_files[file_index], file = fileLog, append = TRUE, col.names = FALSE, row.names = FALSE)
        }
      }
      
      if (length(sample_index) == 0) next
      if (length(sample_index) == 1) {
        index <- sample_index
      } else {
        index <- sample(sample_index, min(batch_size - x_index + 1, length(sample_index)))
      }
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
    
    # if (!is.null(n_gram) & is_lm) {
    #    x <- x[ , 1:(dim(x)[2] - n_gram + 1), ]
    #    if (is.matrix(x)) {
    #      x <- array(x, dim = c(1, dim(x)))
    #    }
    #  }
    
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
    
    return(list(x, y))
  }
}

#' Collect samples from generator and store in rds file.
#' 
#' Repeatedly take one random sample from random file and combine samples. Files are sampled weighted by file size.  
#' 
#' @inheritParams fastaFileGenerator
#' @inheritParams fastaLabelGenerator
#' @inheritParams initializeGenerators
#' @inheritParams labelByFolderGeneratorWrapper
#' @param output_path Where to store output rds file.
dataset_from_gen <- function(output_path,
                             num_samples = 1000,
                             train_type = "lm",
                             output_format = "target_right",
                             path,
                             format = "fasta",
                             maxlen = 250,
                             vocabulary = c("a", "c", "g", "t"),
                             seed = 1234,
                             reverseComplements = FALSE,
                             ambiguous_nuc = "zeros",
                             use_quality_score = FALSE,    
                             padding = FALSE,
                             added_label_path = NULL,
                             add_input_as_seq = NULL,
                             skip_amb_nuc = NULL,
                             concat_seq = NULL,
                             target_len = 1,
                             use_coverage = FALSE,
                             reshape_mode = NULL,
                             set_learning = NULL) {
  
  stopifnot(train_type %in% c("lm", "label_header", "label_folder", "label_csv"))
  
  if (train_type == "lm") {
    
    gen <- fastaFileGenerator(corpus.dir = path,
                              format = format,
                              batch.size = num_samples,
                              maxlen = maxlen,
                              max_iter = 10000,
                              vocabulary = vocabulary,
                              verbose = FALSE,
                              randomFiles = TRUE,
                              step = 1, 
                              showWarnings = FALSE,
                              seed = seed,
                              shuffleFastaEntries = TRUE,
                              numberOfFiles = NULL,
                              fileLog = NULL,
                              reverseComplements = reverseComplements,
                              output_format = output_format,
                              ambiguous_nuc = ambiguous_nuc,
                              use_quality_score = use_quality_score,    
                              proportion_per_file = NULL,
                              padding = padding,
                              added_label_path = added_label_path,
                              add_input_as_seq = add_input_as_seq,
                              skip_amb_nuc = skip_amb_nuc,
                              max_samples = 1,
                              concat_seq = concat_seq,
                              target_len = target_len,
                              file_filter = NULL, 
                              use_coverage = FALSE,
                              proportion_entries = NULL,
                              sample_by_file_size = TRUE)
  }
  
  if (train_type == "label_header" | train_type == "label_csv") {
    
    gen <- fastaLabelGenerator(corpus.dir = path,
                               format = format,
                               batch.size = num_samples,
                               maxlen = maxlen,
                               max_iter = 10000,
                               vocabulary = vocabulary,
                               verbose = FALSE,
                               randomFiles = TRUE,
                               step = 1, 
                               showWarnings = FALSE,
                               seed = seed,
                               shuffleFastaEntries = TRUE,
                               numberOfFiles = NULL,
                               fileLog = NULL,
                               labelVocabulary = labelVocabulary,
                               reverseComplements = reverseComplements,
                               ambiguous_nuc = ambiguous_nuc,
                               proportion_per_file = NULL,
                               read_data = FALSE,
                               use_quality_score = use_quality_score,
                               padding = padding,
                               skip_amb_nuc = skip_amb_nuc,
                               max_samples = 1,
                               concat_seq = concat_seq,
                               added_label_path = added_label_path,
                               add_input_as_seq = add_input_as_seq,
                               target_from_csv = target_from_csv,
                               target_split = NULL,
                               file_filter = NULL,
                               use_coverage = use_coverage,
                               proportion_entries = NULL,
                               sample_by_file_size = TRUE)
  }
  
  if (train_type == "label_folder") {
    
    if (is.null(set_learning)) {
      samples_per_target <- NULL
      new_batch_size <- NULL
      reshape_mode <- NULL
    } else {
      reshape_mode <- set_learning$reshape_mode
      maxlen <- set_learning$maxlen
      samples_per_target <- set_learning$samples_per_target
      new_batch_size <- batch.size
      batch.size <- samples_per_target * batch.size
    }
    
    initializeGenerators(directories = path, format = format, batch.size = num_samples, maxlen = maxlen, vocabulary = vocabulary,
                         verbose = FALSE, randomFiles = TRUE, step = 1, showWarnings = FALSE, seed = seed,
                         shuffleFastaEntries = TRUE, numberOfFiles = NULL, skip_amb_nuc = skip_amb_nuc,
                         fileLog = NULL, reverseComplements = reverseComplements, reverseComplementEncoding = FALSE,
                         val = FALSE, ambiguous_nuc = ambiguous_nuc,
                         proportion_per_file = NULL, read_data = FALSE, use_quality_score = use_quality_score,
                         padding = padding, max_samples = 1, split_seq = FALSE, concat_seq = concat_seq,
                         added_label_path = added_label_path, add_input_as_seq = add_input_as_seq, use_coverage = use_coverage,
                         set_learning = set_learning, proportion_entries = NULL, sample_by_file_size = TRUE)
    
    gen <- labelByFolderGeneratorWrapper(val = FALSE, path = path, new_batch_size = new_batch_size,
                                         samples_per_target = samples_per_target, 
                                         batch.size = batch.size, voc_len = length(vocabulary),
                                         maxlen = maxlen, reshape_mode = reshape_mode)
  }
  
  tensor_list <- gen()
  saveRDS(tensor_list, file = output_path)
}
