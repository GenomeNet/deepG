#' CRISPR data
#' 
#' Example training dataset consisting of a sequence of nucleotides of CRISPR loci
#' Filtered for unambiguous characters and contains only characters in the vocabulary {A,G,G,T
#' }
#' Can be loaded to workspace via data(crispr_sample)
#' @format Large character of 442.41 kB
#' @references \url{http://github.com/philippmuench}
"crispr_sample"

#' CRISPR data subset
#'
#' Example training dataset consisting of a sequence of nucleotides of CRISPR loci
#' Filtered for unambiguous characters and contains only characters in the vocabulary {A,G,G,T
#' }
#' contain all CRISPR loci found in NCBI representative genomes using CRT
#' Can be loaded to workspace via data(crispr_full)
#' @format Large character of 7.50 MB
#' @references \url{http://github.com/philippmuench}
"crispr_full"

#' Parenthesis data
#' 
#' Training dataset of synthetic parenthesis language
#' Can be loaded to workspace via data(parenthesis)
#' @format Large character of 1.00 MB
#' @references \url{http://github.com/philippmuench}
"parenthesis"

#' Ecoli genome
#' 
#' E. coli genome for evaluation
#' Can be loaded to workspace via data(ecoli)
#' @format Large character of 4.64 MB
#' @references \url{https://science.sciencemag.org/content/277/5331/1453.long}
"ecoli"

#' Ecoli subset
#' 
#' subset of the E. coli genome for evaluation
#' Can be loaded to workspace via data(ecoli_small)
#' @format character 326.73 kB
#' @references \url{https://science.sciencemag.org/content/277/5331/1453.long}
"ecoli_small"

#' Genomenet model maxlen 150
#' 
#' Serialized keras model for "bacteria", "virus-no-phage","virus-phage" classification.
#' Can be loaded to workspace via data(self_genomenet_model_maxlen_150).
#' Use keras::unserialize_model to transform to keras model or load with
#' load_model_self_genomenet(maxlen = 150) function call.
#' @format raw
"model_self_genomenet_maxlen_150"

#' Genomenet model maxlen 10000
#' 
#' Serialized keras model for "bacteria", "virus-no-phage","virus-phage" classification.
#' Can be loaded to workspace via data(self_genomenet_model_maxlen_10k)
#' Use keras::unserialize_model to transform to keras model or load with
#' load_model_self_genomenet(maxlen = 10000) function call.
#' @format raw
"model_self_genomenet_maxlen_10k"
