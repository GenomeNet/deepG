#' Randomly select samples from fasta files
#' 
#' Generator \code{\link{generator_fasta_lm}}, \code{\link{generator_fasta_label_header_csv}}
#' or \code{\link{generator_fasta_label_folder}} will randomly choose a consecutive sequence of samples when
#' a \code{max_samples} argument is supplied. \code{generator_random} will choose samples at random.
#'
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_header_csv
#' @inheritParams generator_fasta_label_folder
#' @inheritParams train_model
#' @param number_target_nt Number of target nucleotides for language model.
#' @examplesIf reticulate::py_module_available("tensorflow")
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
#' gen <- generator_random(path = path_input, batch_size = 2,
#'                         vocabulary_label = c('label_a', 'label_b'),
#'                         train_type = 'label_csv',
#'                         maxlen = 5, target_from_csv = target_from_csv)
#' z <- gen()
#' dim(z[[1]])
#' z[[2]]
#' 
#' @returns A generator function.  
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
    if (!is.data.frame(target_from_csv)) {
      output_label_csv <- utils::read.csv2(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
      if (dim(output_label_csv)[2] == 1) {
        output_label_csv <- utils::read.csv(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
      }
    } else {
      output_label_csv <- target_from_csv
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
      
      if (dirname(output_label_csv$file[1]) == ".") {
        # relative path 
        use_basename <- TRUE
        index_basename <- basename(fasta_files[[i]]) %in% output_label_csv$file
      } else {
        # absolute path
        use_basename <- FALSE
        index_basename <- fasta_files[[i]] %in% output_label_csv$file
      }
      index_abs_path <- fasta_files[[i]] %in% output_label_csv$file
      index <- index_basename | index_abs_path
      fasta_files[[i]] <- fasta_files[[i]][index]
      if (length(fasta_files[[i]]) == 0) {
        stop("No overlap between files and 'file' column in target_from_csv")
      }
      
    }
    
    num_files[[i]] <- length(fasta_files[[i]])
    start_ind[[i]] <- (0:(batch_size[[i]] - 1) * seq_len_total) + 1
    if (sample_by_file_size) {
      file_prob[[i]] <- file.info(fasta_files[[i]])$size/sum(file.info(fasta_files[[i]])$size)
    } else {
      file_prob <- NULL
    }
  }
  
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
          if (use_basename) {
            file_row_name <- fasta_files[[p]][file_index] %>% basename 
          } else {
            file_row_name <- fasta_files[[p]][file_index]
          }
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
                                             n_gram_stride = n_gram_stride)
      }
      
      if (train_type != "lm") {
        seq_list[[p]] <- one_hot_sample
      } else {
        xy_list <- slice_tensor_lm(xy = one_hot_sample,
                                   output_format = output_format,
                                   target_len = target_len,
                                   n_gram = n_gram,
                                   n_gram_stride = n_gram_stride,
                                   total_seq_len = seq_len_total,
                                   return_int = return_int)
      }
    }
    
    batch_number <<- batch_number + 1
    
    if (train_type == "lm") {
      if (reshape_xy_bool) {
        xy_list <- f_reshape(x = xy_list$x, y = xy_list$y,
                             reshape_xy = reshape_xy,
                             reshape_x_bool = reshape_x_bool,
                             reshape_y_bool = reshape_y_bool,
                             reshape_sw_bool = FALSE, sw = NULL)
      } 
      
      return(xy_list)
    } 
    
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
      y <- do.call(rbind, target_list) %>% matrix(nrow = batch_size)
      if (!is.null(target_split)) {
        colnames(y) <- vocabulary_label
        y <- slice_tensor(tensor = y, target_split = target_split)
      }
      colnames(y) <- NULL
    }
    
    if (train_type == "label_header") {
      x <- one_hot_sample
      target_int <- unlist(keras::texts_to_sequences(tokenizer_target, unlist(target_list))) - 1
      y  <- keras::to_categorical(target_int, num_classes = length(vocabulary_label))
    }
    
    if (train_type == "masked_lm") {
      if (reshape_xy_bool) {
        if (reshape_x_bool) one_hot_sample$x <- reshape_xy$x(x = one_hot_sample$x, y = one_hot_sample$y)
        if (reshape_y_bool) one_hot_sample$y <- reshape_xy$y(x = one_hot_sample$x, y = one_hot_sample$y)
      }
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
}
