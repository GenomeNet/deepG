#' Random data generator
#' 
#' Creates a random input/target list once and repeatedly returns list. 
#'
#' @inheritParams generator_fasta_lm
#' @param model A keras model.
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
#' 
#' @export
generator_dummy <- function(model, batch_size) {
  
  # sparse loss (TODO: adjust for multi output mixes with sparse/non-sparse loss)
  if (!is.null(model$loss) && 
      any(stringr::str_detect(stringr::str_to_lower(model$loss), "sparse"))) {
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
    x <- array(stats::runif(n = prod(input_dim)), dim = input_dim)
  } else {
    
    for (i in 1:num_input_layers) {
      input_dim <- batch_size
      input_list <- model$input_shape[[i]]
      for (j in 2:length(input_list)) {
        input_dim  <- c(input_dim, input_list[[j]])
      }
      x[[i]] <- array(stats::runif(n = prod(input_dim)), dim = input_dim)
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
    y <- array(stats::runif(n = prod(output_dim)), dim = output_dim)
  } else {
    for (i in 1:num_output_layers) {
      output_dim <- batch_size
      output_list <- model$output_shape[[i]]
      for (j in 2:length(output_list)) {
        output_dim  <- c(output_dim, output_list[[j]])
      }
      if (sparse_loss) output_dim <- output_dim[-length(output_dim)]
      y[[i]] <- array(stats::runif(n = prod(output_dim)), dim = output_dim)
    }
  }
  
  function() {
    return(list(x,y))
  }
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
#' \code{iterations}) and file ending depends on \code{store_format} argument.
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
#' 
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
                             random_sampling = FALSE,
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
                       random_sampling = random_sampling,
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
