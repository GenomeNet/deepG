% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{start_ind_ignore_amb}
\alias{start_ind_ignore_amb}
\title{Helper function for get_start_ind, extracts the start positions of all potential samples (considering step size and vocabulary)}
\usage{
start_ind_ignore_amb(
  seq_vector,
  length_vector,
  maxlen,
  step,
  vocabulary,
  train_mode = "lm"
)
}
\arguments{
\item{seq_vector}{Vector of character sequences.}

\item{length_vector}{Length of sequences in \code{seq_vector}.}

\item{maxlen}{Length of one sample.}

\item{step}{How often to take a sample.}

\item{vocabulary}{Vector of allowed characters in samples.}

\item{train_mode}{"lm" or "label".}
}
\description{
Helper function for get_start_ind, extracts the start positions of all potential samples (considering step size and vocabulary)
}
\keyword{internal}
