% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_model.R
\name{create_model_lstm_cnn_time_dist}
\alias{create_model_lstm_cnn_time_dist}
\title{Create LSTM/CNN network for combining multiple sequences}
\usage{
create_model_lstm_cnn_time_dist(
  maxlen = 50,
  dropout_lstm = 0,
  recurrent_dropout_lstm = 0,
  layer_lstm = NULL,
  layer_dense = c(4),
  solver = "adam",
  learning_rate = 0.001,
  vocabulary_size = 4,
  bidirectional = FALSE,
  stateful = FALSE,
  batch_size = NULL,
  compile = TRUE,
  kernel_size = NULL,
  filters = NULL,
  strides = NULL,
  pool_size = NULL,
  padding = "same",
  dilation_rate = NULL,
  gap_time_dist = NULL,
  use_bias = TRUE,
  zero_mask = FALSE,
  label_smoothing = 0,
  label_noise_matrix = NULL,
  last_layer_activation = "softmax",
  loss_fn = "categorical_crossentropy",
  auc_metric = FALSE,
  f1_metric = FALSE,
  samples_per_target,
  batch_norm_momentum = 0.99,
  verbose = TRUE,
  model_seed = NULL,
  aggregation_sum = TRUE,
  transformer_args = NULL,
  lstm_time_dist = NULL,
  mixed_precision = FALSE,
  mirrored_strategy = NULL
)
}
\arguments{
\item{maxlen}{Length of predictor sequence.}

\item{dropout_lstm}{Fraction of the units to drop for inputs.}

\item{recurrent_dropout_lstm}{Fraction of the units to drop for recurrent state.}

\item{layer_lstm}{Number of cells per network layer. Can be a scalar or vector.}

\item{layer_dense}{Vector specifying number of neurons per dense layer after last LSTM or CNN layer (if no LSTM used).}

\item{solver}{Optimization method, options are \verb{"adam", "adagrad", "rmsprop"} or \code{"sgd"}.}

\item{learning_rate}{Learning rate for optimizer.}

\item{vocabulary_size}{Number of unique character in vocabulary.}

\item{bidirectional}{Use bidirectional wrapper for lstm layers.}

\item{stateful}{Boolean. Whether to use stateful LSTM layer.}

\item{batch_size}{Number of samples that are used for one network update. Only used if \code{stateful = TRUE}.}

\item{compile}{Whether to compile the model.}

\item{kernel_size}{Size of 1d convolutional layers. For multiple layers, assign a vector. (e.g, \code{rep(3,2)} for two layers and kernel size 3)}

\item{filters}{Number of filters. For multiple layers, assign a vector.}

\item{strides}{Stride values. For multiple layers, assign a vector.}

\item{pool_size}{Integer, size of the max pooling windows. For multiple layers, assign a vector.}

\item{padding}{Padding of CNN layers, e.g. \verb{"same", "valid"} or \code{"causal"}.}

\item{dilation_rate}{Integer, the dilation rate to use for dilated convolution.}

\item{gap_time_dist}{Pooling or flatten method after last time distribution wrapper. Same options as for \code{flatten_method} argument
in \link{create_model_transformer} function.}

\item{use_bias}{Boolean. Usage of bias for CNN layers.}

\item{zero_mask}{Boolean, whether to apply zero masking before LSTM layer. Only used if model does not use any CNN layers.}

\item{label_smoothing}{Float in [0, 1]. If 0, no smoothing is applied. If > 0, loss between the predicted
labels and a smoothed version of the true labels, where the smoothing squeezes the labels towards 0.5.
The closer the argument is to 1 the more the labels get smoothed.}

\item{label_noise_matrix}{Matrix of label noises. Every row stands for one class and columns for percentage of labels in that class.
If first label contains 5 percent wrong labels and second label no noise, then

\code{label_noise_matrix <- matrix(c(0.95, 0.05, 0, 1), nrow = 2, byrow = TRUE )}}

\item{last_layer_activation}{Either \code{"sigmoid"} or \code{"softmax"}.}

\item{loss_fn}{Either \code{"categorical_crossentropy"} or \code{"binary_crossentropy"}. If \code{label_noise_matrix} given, will use custom \code{"noisy_loss"}.}

\item{auc_metric}{Whether to add AUC metric.}

\item{f1_metric}{Whether to add F1 metric.}

\item{samples_per_target}{Number of samples to combine for one target.}

\item{batch_norm_momentum}{Momentum for the moving mean and the moving variance.}

\item{verbose}{Boolean.}

\item{model_seed}{Set seed for model parameters in tensorflow if not \code{NULL}.}

\item{aggregation_sum}{Whether to add representations.}

\item{transformer_args}{List of arguments for transformer blocks; see \link{layer_transformer_block_wrapper}.
Additionally, list can contain \code{pool_flatten} argument to apply global pooling or flattening after last transformer block (same options
as \code{flatten_method} argument in \link{create_model_transformer} function).}

\item{lstm_time_dist}{Vector containing number of units per LSTM cell. Applied after time distribution part.}

\item{mixed_precision}{Whether to use mixed precision (https://www.tensorflow.org/guide/mixed_precision).}

\item{mirrored_strategy}{Whether to use distributed mirrored strategy. If NULL, will use distributed mirrored strategy only if >1 GPU available.}
}
\description{
Creates a network consisting of an arbitrary number of CNN, LSTM and dense layers.
Input is a 4D tensor, where axis correspond to:
\enumerate{
\item batch size
\item number of samples in one batch
\item length of one sample
\item size of vocabulary
}
After LSTM/CNN part all representations get aggregated by summation.
Can be used to make single prediction for combination of multiple input sequences. Architecture
is equivalent to \code{\link[=create_model_lstm_cnn_multi_input]{create_model_lstm_cnn_multi_input()}} but instead of multiple input layers with 3D input,
input here in one 4D tensor.
}
\examples{
create_model_lstm_cnn_time_dist(
  maxlen = 50,
  vocabulary_size = 4,
  samples_per_target = 7,
  kernel_size = c(10, 10),
  filters = c(64, 128),
  pool_size = c(2, 2),
  layer_lstm = c(32),
  layer_dense = c(64, 2),
  learning_rate = 0.001)
}
