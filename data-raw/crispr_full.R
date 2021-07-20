library(readr)

crispr_full <- read_lines("crispr_1000.txt")
use_data(crispr_full, overwrite = TRUE)