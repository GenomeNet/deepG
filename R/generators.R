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
#' @examplesIf reticulate::py_module_available("tensorflow")
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
#' @returns A generator function.
#' @export
get_generator <- function(path = NULL,
                          train_type,
                          batch_size,
                          maxlen,
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
                          masked_lm = NULL,
                          val = FALSE,
                          return_int = FALSE,
                          verbose = TRUE,
                          delete_used_files = FALSE,
                          reshape_xy = NULL) {
  
  if (random_sampling) {
    if (use_quality_score) stop("use_quality_score not implemented for random sampling")
    if (read_data) stop("read_data not implemented for random sampling")
    if (!is.null(use_coverage)) stop("use_coverage not implemented for random sampling")
    if (!is.null(add_noise)) stop("add_noise not implemented for random sampling")
  }
  
  if (train_type %in% c("label_rds", "lm_rds") & format != "rds") {
    warning(paste("train_type is", train_type, "but format is not 'rds'"))
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
    #gen <- generator_dummy(model, ifelse(is.null(set_learning), batch_size, new_batch_size))
    gen <- generator_dummy(model, batch_size)
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
      concat_seq = concat_seq,
      reshape_xy = reshape_xy)
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
                              file_filter = file_filter, use_coverage = use_coverage, return_int = return_int,
                              sample_by_file_size = sample_by_file_size, add_noise = add_noise,
                              reshape_xy = reshape_xy)
  }
  
  # label by folder
  if (train_type %in% c("label_folder", "masked_lm") & random_sampling) {
    
    gen <- generator_random(
      train_type = train_type,
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
      concat_seq = concat_seq,
      reshape_xy = reshape_xy)
  } 
  
  if (train_type == "label_folder" & !random_sampling) {
    
    gen_list <- generator_initialize(directories = path, format = format, batch_size = batch_size, maxlen = maxlen, vocabulary = vocabulary,
                                     verbose = verbose, shuffle_file_order = shuffle_file_order, step = step, seed = seed[1],
                                     shuffle_input = shuffle_input, file_limit = file_limit, skip_amb_nuc = skip_amb_nuc,
                                     path_file_log = path_file_log, reverse_complement = reverse_complement,
                                     reverse_complement_encoding = reverse_complement_encoding, return_int = return_int,
                                     ambiguous_nuc = ambiguous_nuc, proportion_per_seq = proportion_per_seq,
                                     read_data = read_data, use_quality_score = use_quality_score, val = val,
                                     padding = padding, max_samples = max_samples, concat_seq = concat_seq,
                                     added_label_path = added_label_path, add_input_as_seq = add_input_as_seq, use_coverage = use_coverage,
                                     set_learning = set_learning, proportion_entries = proportion_entries,
                                     sample_by_file_size = sample_by_file_size, n_gram = n_gram, n_gram_stride = n_gram_stride,
                                     add_noise = add_noise, reshape_xy = reshape_xy)
    
    gen <- generator_fasta_label_folder_wrapper(val = val, path = path, 
                                                batch_size = batch_size, voc_len = length(vocabulary),
                                                gen_list = gen_list,
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
                                        add_noise = add_noise,
                                        reshape_xy = reshape_xy) 
  }
  
  
  if ((train_type == "label_csv" | train_type == "label_header") & !random_sampling) {
    
    gen <- generator_fasta_label_header_csv(path_corpus = path, format = format, batch_size = batch_size, maxlen = maxlen,
                                            vocabulary = vocabulary, verbose = verbose, shuffle_file_order = shuffle_file_order, step = step,
                                            seed = seed[1], shuffle_input = shuffle_input, return_int = return_int,
                                            path_file_log = path_file_log, vocabulary_label = vocabulary_label, reverse_complement = reverse_complement,
                                            ambiguous_nuc = ambiguous_nuc, proportion_per_seq = proportion_per_seq,
                                            read_data = read_data, use_quality_score = use_quality_score, padding = padding,
                                            added_label_path = added_label_path, add_input_as_seq = add_input_as_seq,
                                            skip_amb_nuc = skip_amb_nuc, max_samples = max_samples, concat_seq = concat_seq,
                                            target_from_csv = target_from_csv, target_split = target_split, file_filter = file_filter,
                                            use_coverage = use_coverage, proportion_entries = proportion_entries,
                                            sample_by_file_size = sample_by_file_size, n_gram = n_gram, n_gram_stride = n_gram_stride,
                                            add_noise = add_noise, reverse_complement_encoding = reverse_complement_encoding,
                                            reshape_xy = reshape_xy)
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
      concat_seq = concat_seq,
      reshape_xy = reshape_xy)
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
                         delete_used_files = delete_used_files, reshape_xy = reshape_xy)
    
  }
  
  return(gen)
  
}
