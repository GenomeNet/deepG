library(Biostrings)
ecoli <- Biostrings::readDNAStringSet("data-raw/ecoli.fasta")[[1]]
ecoli <- paste0("|", paste(ecoli, collapse = "-"),"|") 
use_data(ecoli, overwrite = TRUE)
