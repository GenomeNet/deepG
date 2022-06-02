#' Example training dataset consiting of a sequence of nucleotides of CRISPR loci
#' Filtered for unambigous characters and contains only characters in the vocabulary {A,G,G,T
#' }
#' Can be loaded to workspace via data(crispr_sample)
#' @format Large character of 4.9 Mb
#' @references \url{http://github.com/philippmuench}
"crispr_sample"

#' Example training dataset consiting of a sequence of nucleotides of CRISPR loci
#' Filtered for unambigous characters and contains only characters in the vocabulary {A,G,G,T
#' }
#' contain all CRISPR loci found in NCBI representative genomes using CRT
#' Can be loaded to workspace via data(crispr_full)
#' @format Large character of 4.9 Mb
#' @references \url{http://github.com/philippmuench}
"crispr_full"

#' Training dataset of synthetic parenthesis language
#' Can be loaded to workspace via data(parenthesis)
#' @format Large character of 4.9 Mb
#' @references \url{http://github.com/philippmuench}
"parenthesis"

#' E. coli genome for evaluation
#' Can be loaded to workspace via data(ecoli)
#' @format Large character of 4.4 Mb
#' @references \url{https://science.sciencemag.org/content/277/5331/1453.long}
"ecoli"

#' subset of the E. coli genome for evaluation
#' Can be loaded to workspace via data(ecoli_small)
#' @format character
#' @references \url{https://science.sciencemag.org/content/277/5331/1453.long}
"ecoli_small"

#' Serialized keras model for "bacteria", "virus-no-phage","virus-phage" classification.
#' Can be loaded to workspace via data(self_genomenet_model_maxlen_150).
#' Use keras::unserialize_model to transform to keras model or load with
#' load_model_self_genomenet(maxlen = 150) function call.
#' @format raw
"model_self_genomenet_maxlen_150"

#' Serialized keras model for "bacteria", "virus-no-phage","virus-phage" classification.
#' Can be loaded to workspace via data(self_genomenet_model_maxlen_10k)
#' Use keras::unserialize_model to transform to keras model or load with
#' load_model_self_genomenet(maxlen = 10000) function call.
#' @format raw
"model_self_genomenet_maxlen_10k"
