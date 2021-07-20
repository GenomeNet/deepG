#' Get distribution of n-grams
#' 
#' @param fasta_path Path to folder containing fasta files or single fasta file.
#' @param n Size of n gram.
#' @param format "fasta" or "fastq".
#' @param vocabulary Vector of allowed characters, samples outside vocabulary get discarded.
#' @param file_sample If integer, size of random sample of files in \code{fasta_path}.
#' @param step Frequency of samples.
#' @return Returns a matrix with distributions of nucleotides given the previous n nucleotides.
#' @export
n_gram <- function(fasta_path,
                   n = 2,
                   vocabulary = c("A", "C", "G", "T"),
                   format = "fasta", 
                   file_sample = NULL,
                   step = 1,
                   nuc_dist = FALSE) { 
  
  if (endsWith(fasta_path, paste0(".", format))) {
    num_files <- 1 
    fasta_files <- fasta_path   
  } else {
    fasta_files <- list.files(
      path = xfun::normalize_path(fasta_path),
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
    
    fasta_file <-  Biostrings::readDNAStringSet(fasta_files[i], format = format)
    seq_vector <- as.data.frame(fasta_file)$x
    start_ind <- getStartInd(seq_vector = seq_vector,
                             length_vector = as.data.frame(fasta_file@ranges)$width,
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
    
    # TODO: check if freq correctly allignes with grams and targets
    stopifnot(all(freq_df$gram == table_df$gram) & all(freq_df$targets == table_df$targets))
    
    freq_df$freq <- freq_df$freq + table_df$Freq
  }
  
  dist_matrix <- df_to_distribution_matrix(freq_df, vocabulary = vocabulary)
  dist_matrix
}

#' Helper function for \code{n_gram}
#'
#' @export
df_to_distribution_matrix <- function(freq_df, vocabulary = c("A", "C", "G", "T")) {
  stopifnot(names(freq_df) == c("gram", "targets", "freq"))
  gram_levels <- levels(freq_df$gram)
  num_levels <- length(gram_levels)
  dist_matrix <- matrix(0, nrow = num_levels, ncol = length(vocabulary))
  dist_matrix <- as.data.frame(dist_matrix)
  #freq_df$percentage <- freq_df$freq/sum(freq_df$freq) 
  rownames(dist_matrix) <- as.character(freq_df$gram[1:nrow(dist_matrix)])
  colnames(dist_matrix) <- vocabulary
  # for (gram_level in gram_levels) {
  #   df_subset <- freq_df[freq_df$gram == gram_level, ]
  #   for (nuc in vocabulary) {
  #     dist_matrix[gram_level, nuc] <- df_subset %>% dplyr::filter(targets == nuc) %>% dplyr::select(freq)
  #   }
  # }
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

#' predict the next nucleotide using n-gram
#' 
#' @param fasta_path Path to folder containing fasta files or single fasta file.
#' @param distribution_matrix Output of \code{n_gram} function.
#' @param default_pred Either character from vocabulary or "random". Will be used as prediction if certain n-gram did not appear before.
#' If "random" assign random prediction.   
#' @param vocabulary Vector of allowed characters, samples outside vocabulary get discarded.
#' @param file_sample If integer, size of random sample of files in \code{fasta_path}.
#' @param format "fasta" or "fastq".
#' @param return_data_frames Boolean, whether to return data frame with input, predictions, target position and true target. 
#' 
#' @examples 
#' \dontrun{
#' n_gram_dist <- n_gram(fasta_path = "/net/sgi/genomenet/data/ncbi_genomes/train",
#'                       n = 2,
#'                       file_sample = NULL,
#'                       nuc_dist = TRUE)
#' 
#' # use distribution matrix to make predictions for one file
#' predictions <- predict_with_n_gram(fasta_path = "/net/sgi/genomenet/data/ncbi_genomes/test/GCF_001941465.1_ASM194146v1_genomic.fasta",
#'                                    distribution_matrix = n_gram_dist)
#' 
#' # show accuracy
#' predictions$accuracy
#' }
#' @export
predict_with_n_gram <- function(fasta_path, distribution_matrix, default_pred = "random", vocabulary = c("A", "C", "G", "T"),
                                file_sample = NULL, format = "fasta", return_data_frames = FALSE) {
  
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
  
  if (endsWith(fasta_path, paste0(".", format))) {
    num_files <- 1 
    fasta_files <- fasta_path   
  } else {
    fasta_files <- list.files(
      path = xfun::normalize_path(fasta_path),
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
    
    fasta_file <-  Biostrings::readDNAStringSet(fasta_files[i], format = format)
    seq_vector <- as.data.frame(fasta_file)$x
    start_ind <- getStartInd(seq_vector = seq_vector,
                             length_vector = as.data.frame(fasta_file@ranges)$width,
                             maxlen = n, step = 1, train_mode = "lm")
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
    gram_df <- gram_df[complete.cases(gram_df), ]
    
    pred_df <- dplyr::left_join(gram_df, model, by = "gram")
    names(pred_df)[2] <- "true"
    if (return_data_frames) {
      pred_df_list[[i]] <- list(pred_df, accurcy = sum(pred_df$true == pred_df$pred)/nrow(pred_df))
    } else {
      pred_df_list[[i]] <- list(accurcy = sum(pred_df$true == pred_df$pred)/nrow(pred_df))
    }
  }
  return(pred_df_list)
}


