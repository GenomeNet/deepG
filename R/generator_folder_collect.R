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
#' @examplesIf reticulate::py_module_available("tensorflow")
#' # create two folders with dummy fasta files
#' path_input_1 <- tempfile()
#' dir.create(path_input_1)
#' create_dummy_data(file_path = path_input_1, num_files = 2, seq_length = 5,
#'                   num_seq = 2, vocabulary = c("a", "c", "g", "t"))
#' path_input_2 <- tempfile()
#' dir.create(path_input_2)
#' create_dummy_data(file_path = path_input_2, num_files = 3, seq_length = 7,
#'                   num_seq = 5, vocabulary = c("a", "c", "g", "t"))
#' 
#' gen_list <- generator_initialize(directories = c(path_input_1, path_input_1),
#'                                         batch_size = 4, maxlen = 5)
#' z1 <- gen_list[[1]]()
#' z1[[1]]
#' z1[[2]]
#' 
#' @returns List of generator function. 
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
                                 return_int = FALSE,
                                 reshape_xy = NULL) {
  
  num_class <- length(directories)
  
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
  argg[["num_class"]] <- NULL
  argg[["reshape_xy"]] <- NULL
  
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
    gen_train_list <- list()
    # create generator for every folder
    for (i in 1:num_class) {
      numberedGen <- paste0("gen", i)
      gen_train_list[[numberedGen]] <- generator_fasta_label_folder(path_corpus = directories[[i]],
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
                                                                    return_int = return_int,
                                                                    reshape_xy = reshape_xy)
      
    }
  } else {
    gen_val_list <- list()
    # create generator for every folder
    for (i in 1:num_class) {
      # different names for validation generators
      numberedGenVal <- paste0("genVal", as.character(i))
      gen_val_list[[numberedGenVal]] <- generator_fasta_label_folder(path_corpus = directories[[i]],
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
                                                                     return_int = return_int,
                                                                     reshape_xy = reshape_xy)
    }
  }
  
  if (!val) {
    gen_train_list
  } else {
    gen_val_list
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
#' @param voc_len Length of vocabulary.
#' @param gen_list List of generator functions.
#' @examplesIf reticulate::py_module_available("tensorflow")
#' # create two folders with dummy fasta files
#' path_input_1 <- tempfile()
#' dir.create(path_input_1)
#' create_dummy_data(file_path = path_input_1, num_files = 2, seq_length = 5,
#'                   num_seq = 2, vocabulary = c("a", "c", "g", "t"))
#' path_input_2 <- tempfile()
#' dir.create(path_input_2)
#' create_dummy_data(file_path = path_input_2, num_files = 3, seq_length = 7,
#'                   num_seq = 5, vocabulary = c("a", "c", "g", "t"))
#' 
#' maxlen <- 5
#' p <- c(path_input_1, path_input_1)
#' gen_list <- generator_initialize(directories = p,
#'                                  batch_size = 4, maxlen = maxlen)
#' gen <- generator_fasta_label_folder_wrapper(val = FALSE, batch_size = 8,
#'                                             path = p, voc_len = 4, 
#'                                             maxlen = maxlen,
#'                                             gen_list = gen_list)
#' z <- gen()
#' dim(z[[1]])
#' z[[2]]
#' 
#' @returns A generator function.  
#' @export
generator_fasta_label_folder_wrapper <- function(val, 
                                                 batch_size = NULL,
                                                 path = NULL, voc_len = NULL, 
                                                 maxlen = NULL,
                                                 gen_list = NULL, 
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
      subBatchTrain <- gen_list[["gen1"]]()
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
          subBatchTrain <- gen_list[[paste0("gen", i)]]()
          yTrain <- abind::abind(yTrain, subBatchTrain[[2]], along = 1)
          
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
      subBatchVal <- gen_list[["genVal1"]]()
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
          subBatchVal <- gen_list[[paste0("genVal", i)]]()
          yVal <- abind::abind(yVal, subBatchVal[[2]], along = 1)
          
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
                            buffer_len = buffer_len,
                            reshape_mode = reshape_mode)
        return(l)
      } else {
        return(list(X = xVal, Y = yVal))
      }
    }
  }
  
}
