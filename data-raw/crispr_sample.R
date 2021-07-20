library(readr)

crispr_sample <- read_lines("crispr_sample.txt")
use_data(crispr_sample, overwrite = TRUE)