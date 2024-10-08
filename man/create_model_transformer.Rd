% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_model_transformer.R
\name{create_model_transformer}
\alias{create_model_transformer}
\title{Create transformer model}
\usage{
create_model_transformer(
  maxlen,
  vocabulary_size = 4,
  embed_dim = 64,
  pos_encoding = "embedding",
  head_size = 4L,
  num_heads = 5L,
  ff_dim = 8,
  dropout = 0,
  n = 10000,
  layer_dense = 2,
  dropout_dense = NULL,
  flatten_method = "flatten",
  last_layer_activation = "softmax",
  loss_fn = "categorical_crossentropy",
  solver = "adam",
  learning_rate = 0.01,
  label_noise_matrix = NULL,
  bal_acc = FALSE,
  f1_metric = FALSE,
  auc_metric = FALSE,
  label_smoothing = 0,
  verbose = TRUE,
  model_seed = NULL,
  mixed_precision = FALSE,
  mirrored_strategy = NULL
)
}
\arguments{
\item{maxlen}{Length of predictor sequence.}

\item{vocabulary_size}{Number of unique character in vocabulary.}

\item{embed_dim}{Dimension for token embedding. No embedding if set to 0. Should be used when input is not one-hot encoded
(integer sequence).}

\item{pos_encoding}{Either \code{"sinusoid"} or \code{"embedding"}. How to add positional information.
If \code{"sinusoid"}, will add sine waves of different frequencies to input.
If \code{"embedding"}, model learns positional embedding.}

\item{head_size}{Dimensions of attention key.}

\item{num_heads}{Number of attention heads.}

\item{ff_dim}{Units of first dense layer after attention blocks.}

\item{dropout}{Vector of dropout rates after attention block(s).}

\item{n}{Frequency of sine waves for positional encoding. Only applied if \code{pos_encoding = "sinusoid"}.}

\item{layer_dense}{Vector specifying number of neurons per dense layer after last LSTM or CNN layer (if no LSTM used).}

\item{dropout_dense}{Dropout for dense layers.}

\item{flatten_method}{How to process output of last attention block. Can be \code{"max_ch_first"}, \code{"max_ch_last"}, \code{"average_ch_first"},
\code{"average_ch_last"}, \code{"both_ch_first"}, \code{"both_ch_last"}, \code{"all"}, \code{"none"} or \code{"flatten"}.
If \code{"average_ch_last"} /  \code{"max_ch_last"}  or \code{"average_ch_first"} / \code{"max_ch_first"}, will apply global average/max pooling.
\verb{_ch_first} / \verb{_ch_last} to decide along which axis. \code{"both_ch_first"} / \code{"both_ch_last"} to use max and average together. \code{"all"} to use all 4
global pooling options together. If \code{"flatten"}, will flatten output after last attention block. If \code{"none"} no flattening applied.}

\item{last_layer_activation}{Activation function of output layer(s). For example \code{"sigmoid"} or \code{"softmax"}.}

\item{loss_fn}{Either \code{"categorical_crossentropy"} or \code{"binary_crossentropy"}. If \code{label_noise_matrix} given, will use custom \code{"noisy_loss"}.}

\item{solver}{Optimization method, options are \verb{"adam", "adagrad", "rmsprop"} or \code{"sgd"}.}

\item{learning_rate}{Learning rate for optimizer.}

\item{label_noise_matrix}{Matrix of label noises. Every row stands for one class and columns for percentage of labels in that class.
If first label contains 5 percent wrong labels and second label no noise, then

\code{label_noise_matrix <- matrix(c(0.95, 0.05, 0, 1), nrow = 2, byrow = TRUE )}}

\item{bal_acc}{Whether to add balanced accuracy.}

\item{f1_metric}{Whether to add F1 metric.}

\item{auc_metric}{Whether to add AUC metric.}

\item{label_smoothing}{Float in [0, 1]. If 0, no smoothing is applied. If > 0, loss between the predicted
labels and a smoothed version of the true labels, where the smoothing squeezes the labels towards 0.5.
The closer the argument is to 1 the more the labels get smoothed.}

\item{verbose}{Boolean.}

\item{model_seed}{Set seed for model parameters in tensorflow if not \code{NULL}.}

\item{mixed_precision}{Whether to use mixed precision (https://www.tensorflow.org/guide/mixed_precision).}

\item{mirrored_strategy}{Whether to use distributed mirrored strategy. If NULL, will use distributed mirrored strategy only if >1 GPU available.}
}
\value{
A keras model implementing transformer architecture.
}
\description{
Creates transformer network for classification. Model can consist of several stacked attention blocks.
}
\examples{
\dontshow{if (reticulate::py_module_available("tensorflow")) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
model <- create_model_transformer(maxlen = 50,
                                  head_size=c(10,12),
                                  num_heads=c(7,8),
                                  ff_dim=c(5,9), 
                                  dropout=c(0.3, 0.5))
\dontshow{\}) # examplesIf}
}
