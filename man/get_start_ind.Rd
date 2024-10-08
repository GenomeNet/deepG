% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{get_start_ind}
\alias{get_start_ind}
\title{Computes start position of samples}
\usage{
get_start_ind(
  seq_vector,
  length_vector,
  maxlen,
  step,
  train_mode = "label",
  discard_amb_nuc = FALSE,
  vocabulary = c("A", "C", "G", "T")
)
}
\arguments{
\item{seq_vector}{Vector of character sequences.}

\item{length_vector}{Length of sequences in \code{seq_vector}.}

\item{maxlen}{Length of one predictor sequence.}

\item{step}{Distance between samples from one entry in \code{seq_vector}.}

\item{train_mode}{Either \code{"lm"} for language model or \code{"label"} for label classification.}

\item{discard_amb_nuc}{Whether to discard all samples that contain characters outside vocabulary.}

\item{vocabulary}{Vector of allowed characters. Characters outside vocabulary get encoded as specified in \code{ambiguous_nuc}.}
}
\value{
A numeric vector.
}
\description{
Helper function for data generators.
Computes start positions in sequence where samples can be extracted, given maxlen, step size and ambiguous nucleotide constraints.
}
\examples{
seq_vector <- c("AAACCCNNNGGGTTT")
get_start_ind(
  seq_vector = seq_vector,
  length_vector = nchar(seq_vector),
  maxlen = 4,
  step = 2,
  train_mode = "label",
  discard_amb_nuc = TRUE,
  vocabulary = c("A", "C", "G", "T"))
  
}
