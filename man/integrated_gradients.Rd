% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{integrated_gradients}
\alias{integrated_gradients}
\title{Compute integrated gradients}
\usage{
integrated_gradients(
  m_steps = 50,
  baseline_type = "zero",
  input_seq,
  target_class_idx,
  model,
  pred_stepwise = FALSE,
  num_baseline_repeats = 1
)
}
\arguments{
\item{m_steps}{Number of steps between baseline and original input.}

\item{baseline_type}{Baseline sequence, either \code{"zero"} for all zeros or \code{"shuffle"} for random permutation of \code{input_seq}.}

\item{input_seq}{Input tensor.}

\item{target_class_idx}{Index of class to compute gradient for}

\item{model}{Model to compute gradient for.}

\item{pred_stepwise}{Whether to do predictions with batch size 1 rather than all at once. Can be used if
input is too big to handle at once. Only supported for single input layer.}

\item{num_baseline_repeats}{Number of different baseline estimations if baseline_type is \code{"shuffle"} (estimate integrated
gradient repeatedly for different shuffles). Final result is average of \code{num_baseline} single calculations.}
}
\value{
A tensorflow tensor.
}
\description{
Computes integrated gradients scores for model and an input sequence.
This can be used to visualize what part of the input is import for the models decision.
Code is R implementation of python code from \href{https://www.tensorflow.org/tutorials/interpretability/integrated_gradients}{here}.
Tensorflow implementation is based on this \href{https://arxiv.org/abs/1703.01365}{paper}.
}
\examples{
\dontshow{if (reticulate::py_module_available("tensorflow")) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
library(reticulate)
model <- create_model_lstm_cnn(layer_lstm = 8, layer_dense = 3, maxlen = 20, verbose = FALSE)
random_seq <- sample(0:3, 20, replace = TRUE)
input_seq <- array(keras::to_categorical(random_seq), dim = c(1, 20, 4))
integrated_gradients(
  input_seq = input_seq,
  target_class_idx = 3,
  model = model)
  
\dontshow{\}) # examplesIf}
}
