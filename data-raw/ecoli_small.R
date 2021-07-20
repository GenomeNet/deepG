library(Biostrings)
ecoli_small <- Biostrings::readDNAStringSet("data-raw/ecoli_small.fasta")[[1]]
ecoli_small <- paste0("|", paste(ecoli_small, collapse = "-"),"|") 
use_data(ecoli_small, overwrite = TRUE)
