#' Write random sequences to fasta file
#' 
#' Create random sequences from predefined vocabulary and write to fasta file.
#' 
#' @param file_path Output directory; can also be a file name but only possible if \code{write_to_file_path = TRUE} and 
#' \code{num_files = 1}).
#' @param num_files Number of files to create.
#' @param header Fasta header name. 
#' @param seq_length Length of one sequence.
#' @param num_seq Number of sequences per file.
#' @param fasta_name_start Beginning string of file name. Output files are named fasta_name_start + _i.fasta where i is an integer index.
#' @param write_to_file_path Whether to write output directly to \code{file_path}, i.e. file_path is not a directory.
#' @param prob Probabiltiy of each character in the \code{vocabulary} to be sampled. If `NULL` each character has same probability.
#' @param vocabulary Set of characters to sample sequences from. 
#' @examples
#' path_output <- tempfile()
#' dir.create(path_output)
#' create_dummy_data(file_path = path_output,
#'                   num_files = 3,
#'                   seq_length = 11, 
#'                   num_seq = 5,                   
#'                   vocabulary = c("a", "c", "g", "t"))
#' list.files(path_output)                   
#' @export
create_dummy_data <- function(file_path,
                              num_files,
                              header = "header", 
                              seq_length, 
                              num_seq,
                              fasta_name_start = "file",
                              write_to_file_path = FALSE,
                              prob = NULL,
                              vocabulary = c("a", "c", "g", "t")) {
  
  if (!is.null(prob)) {
    stopifnot(length(prob) == length(vocabulary))
    stopifnot(sum(prob) == 1)
  }
  
  if (write_to_file_path) {
    stopifnot(num_files == 1)
  }
  
  if (!dir.exists(file_path) & !write_to_file_path) {
    dir.create(file_path)
  }
  
  for (i in 1:num_files){
    df <- data.frame(Header = NULL, Sequence = NULL)
    for (j in 1:num_seq) {
      nuc_seq <- sample(vocabulary, seq_length, replace = TRUE, prob = prob) %>% paste(collapse = "")
      new_row <- data.frame(Header = header, Sequence = nuc_seq)
      df <- rbind(df, new_row)
    }
    df$Header <- as.character(df$Header)
    df$Sequence <- as.character(df$Sequence)
    if (write_to_file_path) {
      out.file <- file_path
    } else {
      file_name <- paste0(fasta_name_start, "_", i, ".fasta")
      out.file <- file.path(file_path, file_name)
    }
    microseq::writeFasta(fdta = dplyr::as_tibble(df), 
                         out.file =  out.file)
  }
}
