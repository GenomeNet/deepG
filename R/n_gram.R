#' Get distribution of n-grams
#' 
#' Get distribution of next character given previous n nucleotides.
#'
#' @inheritParams generator_fasta_lm
#' @param path_input Path to folder containing fasta files or single fasta file.
#' @param n Size of n gram.
#' @param vocabulary Vector of allowed characters, samples outside vocabulary get discarded.
#' @param file_sample If integer, size of random sample of files in \code{path_input}.
#' @param nuc_dist Nucleotide distribution.
#' @return Returns a matrix with distributions of nucleotides given the previous n nucleotides.
#' @examples
#' temp_dir <- tempfile()
#' dir.create(temp_dir)
#' create_dummy_data(file_path = temp_dir,
#'                   num_files = 3,
#'                   seq_length = 80,
#'                   vocabulary = c("A", "C", "G", "T"),
#'                   num_seq = 2)
#' 
#' m <- n_gram_dist(path_input = temp_dir,
#'                  n = 3,
#'                  step = 1,
#'                  nuc_dist = FALSE)
#' head(round(m, 2))
#' @returns A data frame of n-gram predictions.
#' @export
n_gram_dist <- function(path_input,
                        n = 2,
                        vocabulary = c("A", "C", "G", "T"),
                        format = "fasta",
                        file_sample = NULL,
                        step = 1,
                        nuc_dist = FALSE) {
  
  if (endsWith(path_input, paste0(".", format))) {
    num_files <- 1
    fasta_files <- path_input
  } else {
    fasta_files <- list.files(
      path = path_input,
      pattern = paste0("\\.", format, "$"),
      full.names = TRUE)
    num_files <- length(fasta_files)
  }
  
  # take random subset of files
  if (!is.null(file_sample)){
    fasta_files <- sample(fasta_files)[1:min(file_sample, length(fasta_files))]
    num_files <- length(fasta_files)
  }
  
  l <- vector("list")
  for (i in 1:n){
    l[[i]] <-  vocabulary
  }
  label_df <- apply(expand.grid(l), 2, as.character)
  labels <- vector("character")
  for (i in 1:nrow(label_df)){
    labels[i] <- paste(label_df[i, ], collapse = "")
  }
  #labels
  
  targets <- vector("character")
  for (i in 1:length(vocabulary)){
    targets <- c(targets, rep(vocabulary[i], length(labels)))
  }
  gram <- rep(labels, length(vocabulary))
  freq <- rep(0, length(labels) * length(vocabulary))
  freq_df <- data.frame(gram, targets, freq)
  nuc_table <- vector("list")
  
  for (i in 1:num_files) {
    
    if (format == "fasta") {
      fasta_file <-  microseq::readFasta(fasta_files[i])
      
    } 
    if (format == "fastq") {
      fasta_file <-  microseq::readFastq(fasta_files[i])
    } 
    
    seq_vector <- fasta_file$Sequence
    start_ind <- get_start_ind(seq_vector = seq_vector,
                               length_vector = nchar(seq_vector),
                               maxlen = n, step = step, train_mode = "lm")
    nuc_seq <- paste(seq_vector, collapse = "")
    split_seq <- strsplit(nuc_seq, "")[[1]]
    nuc_seq_length <- nchar(nuc_seq)
    gram <- split_seq[1 : (nuc_seq_length - n)]
    if (n > 1){
      for (j in 2:n){
        gram <- paste0(gram, split_seq[j : (nuc_seq_length - n + j - 1)])
      }
    }
    targets <- split_seq[(n + 1) : nuc_seq_length]
    
    # remove sequences with overlapping fasta entries
    gram <- gram[start_ind]
    targets <- targets[start_ind]
    
    # remove sequences with ambiguous nucleotides
    amb_pos_gram <- c(1:(length(gram)))[stringr::str_detect(gram, paste0("[^", paste0(vocabulary, collapse = ""), "]"))]
    amb_pos_targets <- c(1:(length(gram)))[stringr::str_detect(targets, paste0("[^", paste0(vocabulary, collapse = ""), "]"))]
    amb_pos <- union(amb_pos_gram, amb_pos_targets)
    if (length(amb_pos) > 0){
      gram <- gram[-amb_pos]
      targets <- targets[-amb_pos]
    }
    
    gram_df <- data.frame(gram = factor(gram, levels = labels),
                          targets = factor(targets, levels = vocabulary))
    table_df <- as.data.frame(table(gram_df))
    
    stopifnot(all(freq_df$gram == table_df$gram) & all(freq_df$targets == table_df$targets))
    
    freq_df$freq <- freq_df$freq + table_df$Freq
  }
  
  dist_matrix <- df_to_distribution_matrix(freq_df, vocabulary = vocabulary)
  dist_matrix
}

df_to_distribution_matrix <- function(freq_df, vocabulary = c("A", "C", "G", "T")) {
  
  stopifnot(names(freq_df) == c("gram", "targets", "freq"))
  gram_levels <- levels(factor(freq_df$gram))
  num_levels <- length(gram_levels)
  dist_matrix <- matrix(0, nrow = num_levels, ncol = length(vocabulary))
  dist_matrix <- as.data.frame(dist_matrix)
  rownames(dist_matrix) <- as.character(freq_df$gram[1:nrow(dist_matrix)])
  colnames(dist_matrix) <- vocabulary
 
  for (nuc in vocabulary){
    nuc_column <- freq_df %>% dplyr::filter(targets == nuc) %>% dplyr::select(gram, freq)
    stopifnot(nuc_column$gram == rownames(dist_matrix))
    dist_matrix[ , nuc] <- nuc_column$freq
  }
  dist_matrix$sum <- apply(dist_matrix, 1, sum)
  non_zero <- dist_matrix$sum != 0
  for (nuc in vocabulary) {
    dist_matrix[non_zero, nuc] <- dist_matrix[non_zero, nuc]/dist_matrix$sum[non_zero]
  }
  dist_matrix[ , vocabulary]
}

#' Predict the next nucleotide using n-gram
#'
#' Predict the next nucleotide using n-gram. 
#'
#' @inheritParams generator_fasta_lm
#' @param path_input Path to folder containing fasta files or single fasta file.
#' @param distribution_matrix A data frame containing frequency of next nucleotide given the previous n nucleotides (output of \code{\link{n_gram_dist}} function).
#' @param default_pred Either character from vocabulary or `"random"`. Will be used as prediction if certain n-gram did not appear before.
#' If `"random"` assign random prediction.
#' @param vocabulary Vector of allowed characters, samples outside vocabulary get discarded.
#' @param file_sample If integer, size of random sample of files in \code{path_input}.
#' @param return_data_frames Boolean, whether to return data frame with input, predictions, target position and true target.
#'
#' @examples
#' # create dummy fasta files
#' temp_dir <- tempfile()
#' dir.create(temp_dir)
#' create_dummy_data(file_path = temp_dir,
#'                   num_files = 3,
#'                   seq_length = 8,
#'                   vocabulary = c("A", "C", "G", "T"),
#'                   num_seq = 2)
#' 
#' m <- n_gram_dist(path_input = temp_dir,
#'                  n = 3,
#'                  step = 1,
#'                  nuc_dist = FALSE)
#' 
#' # use distribution matrix to make predictions for one file
#' predictions <- predict_with_n_gram(path_input = list.files(temp_dir, full.names = TRUE)[1], 
#'                                    distribution_matrix = m)
#' 
#' # show accuracy
#' predictions[[1]]
#' 
#' @returns List of prediction evaluations.
#' @export
predict_with_n_gram <- function(path_input, distribution_matrix, default_pred = "random", vocabulary = c("A", "C", "G", "T"),
                                file_sample = NULL, format = "fasta", return_data_frames = FALSE, step = 1) {
  
  n <- nchar(rownames(distribution_matrix)[1])
  pred_int <- apply(distribution_matrix, 1, which.max)
  # predict most common nucleotide if gram did not appear before
  sum_columns <- apply(distribution_matrix, 2, sum)
  zero_rows <- which(sum_columns == 0)
  if (default_pred == "random") {
    random_pred <- sample(1:length(vocabulary), length(zero_rows), replace = TRUE)
    pred_int[zero_rows] <- random_pred
  } else {
    pred_int[zero_rows] <- which(vocabulary == default_pred)
  }
  # integer to nucleotide
  pred <- vector("character")
  for (i in 1:length(pred_int)){
    pred[i] <- vocabulary[pred_int[i]]
  }
  
  model <- data.frame(gram = rownames(distribution_matrix), pred = pred)
  
  if (endsWith(path_input, paste0(".", format))) {
    num_files <- 1
    fasta_files <- path_input
  } else {
    fasta_files <- list.files(
      path = path_input,
      pattern = paste0("\\.", format, "$"),
      full.names = TRUE)
    num_files <- length(fasta_files)
  }
  
  # take random subset of files
  if (!is.null(file_sample)){
    fasta_files <- sample(fasta_files)[1 : min(file_sample, length(fasta_files))]
    num_files <- length(fasta_files)
  }
  
  labels <- rownames(distribution_matrix)
  
  pred_df_list <- vector("list")
  
  for (i in 1:num_files) {
    
    if (format == "fasta") {
      fasta_file <-  microseq::readFasta(fasta_files[i])
      
    } 
    if (format == "fastq") {
      fasta_file <-  microseq::readFastq(fasta_files[i])
    } 
    
    seq_vector <- fasta_file$Sequence
    start_ind <- get_start_ind(seq_vector = seq_vector,
                               length_vector = nchar(seq_vector),
                               maxlen = n, step = step, train_mode = "lm")
    nuc_seq <- paste(seq_vector, collapse = "")
    split_seq <- strsplit(nuc_seq, "")[[1]]
    
    nuc_seq_length <- nchar(nuc_seq)
    gram <- split_seq[1 : (nuc_seq_length - n)]
    if (n > 1){
      for (j in 2:n){
        gram <- paste0(gram, split_seq[j : (nuc_seq_length - n + j - 1)])
      }
    }
    targets <- split_seq[(n + 1) : nuc_seq_length]
    
    # remove sequences with overlapping fasta entries
    gram <- gram[start_ind]
    targets <- targets[start_ind]
    gram_df <- data.frame(gram = factor(gram, levels = labels),
                          targets = factor(targets, levels = vocabulary),
                          target_pos = start_ind + n)
    
    # remove sequences with ambiguous nucleotides
    gram_df <- gram_df[stats::complete.cases(gram_df), ]
    
    pred_df <- dplyr::left_join(gram_df, model, by = "gram")
    names(pred_df)[2] <- "true"
    if (return_data_frames) {
      pred_df_list[[i]] <- list(pred_df, accuracy = sum(pred_df$true == pred_df$pred)/nrow(pred_df))
    } else {
      pred_df_list[[i]] <- list(accuracy = sum(pred_df$true == pred_df$pred)/nrow(pred_df))
    }
  }
  
  return(pred_df_list)
}
