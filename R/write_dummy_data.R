#' create random sequences and write to fasta file
#' 
#' @export
create_dummy_data <- function(file_path, num_files, seq_length, num_seq, vocabulary) {
  dir.create(file_path)
  for (i in 1:num_files){
    df <- data.frame(Header = NULL, Sequence = NULL)
    for (j in 1:num_seq) {
      nuc_seq <- sample(vocabulary, seq_length, replace = TRUE) %>% paste(collapse = "")
      header <- paste0("label_", ifelse(j %% 2 == 0, "a", "b"))
      new_row <- data.frame(Header = header, Sequence = nuc_seq)
      df <- rbind(df, new_row)
    }
    df$Header <- as.character(df$Header)
    df$Sequence <- as.character(df$Sequence)
    file_name <- paste0("file_", i, ".fasta")
    microseq::writeFasta(fdta = dplyr::as_tibble(df), 
                         out.file =  file.path(file_path, file_name))
  }
}
