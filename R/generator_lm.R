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
#' @param return_int Whether to return integer encoding or one-hot encoding.
#' @param reshape_xy Can be a list of functions to apply to input and/or target. List elements (containing the reshape functions)
#'  must be called x for input or y for target. 
#' @rawNamespace import(data.table, except = c(first, last, between))
#' @importFrom magrittr %>%
#' @examples
#' # create dummy fasta files
#' path_input_1 <- tempfile()
#' dir.create(path_input_1)
#' create_dummy_data(file_path = path_input_1,
#'                   num_files = 2,
#'                   seq_length = 8,
#'                   num_seq = 1,
#'                   vocabulary = c("a", "c", "g", "t"))
#' 
#' gen <- generator_fasta_lm(path_corpus = path_input_1, batch_size = 2,
#'                                    maxlen = 7)
#' z <- gen()
#' dim(z[[1]])
#' z[[2]]
#' 
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
                               return_int = FALSE,
                               reshape_xy = NULL) {
  
  
  ##TODO: add check for n-gram and option for stride
  # if (!is.null(n_gram) & !(any(n_gram_stride == c(n_gram, 1)))) {
  #   stop("When using language model with n_gram encoding, n_gram_stride must be 1 or equal to n_gram")
  # } 
  if (!is.null(n_gram)) {
    # maxlen_n_gram <- ceiling((maxlen - n_gram + 1)/n_gram_stride)
    # target_len_n_gram <- ceiling((target_len - n_gram + 1)/n_gram_stride) 
    if (!n_gram_stride == n_gram) {
      stop("When using train_type='lm' with n_gram encoding, n_gram_stride must be equal to n_gram.")
    }  
  } # else {
  #   maxlen_n_gram <- maxlen
  #   target_len_n_gram <- target_len 
  # }
  
  if (!is.null(reshape_xy)) {
    reshape_xy_bool <- TRUE
    reshape_x_bool <- ifelse(is.null(reshape_xy$x), FALSE, TRUE)
    reshape_y_bool <- ifelse(is.null(reshape_xy$y), FALSE, TRUE)
  } else {
    reshape_xy_bool <- FALSE
  }
  
  
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
                               # maxlen_n_gram = maxlen_n_gram,
                               # target_len_n_gram = target_len_n_gram, 
                               n_gram_stride = n_gram_stride,
                               total_seq_len = total_seq_len,
                               return_int = return_int)
    
    if (reshape_xy_bool) {
      if (reshape_x_bool) xy_list$x <- reshape_xy$x(xy_list$x)
      if (reshape_y_bool) xy_list$y <- reshape_xy$y(xy_list$y)
    }
    
    if (is.null(added_label_path)) {
      return(xy_list)
    } else {
      return(list(append(added_input, list(xy_list$x)), xy_list$y))
    }
    
    # add dim for batch size 1
    
  }
}
