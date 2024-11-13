#' Create transformer model
#'
#' Creates transformer network for classification. Model can consist of several stacked attention blocks.
#' 
#' @inheritParams keras::layer_multi_head_attention
#' @inheritParams create_model_lstm_cnn
#' @param pos_encoding Either `"sinusoid"` or `"embedding"`. How to add positional information.
#' If `"sinusoid"`, will add sine waves of different frequencies to input.
#' If `"embedding"`, model learns positional embedding.
#' @param embed_dim Dimension for token embedding. No embedding if set to 0. Should be used when input is not one-hot encoded
#' (integer sequence).
#' @param head_size Dimensions of attention key.
#' @param n Frequency of sine waves for positional encoding. Only applied if `pos_encoding = "sinusoid"`.
#' @param ff_dim Units of first dense layer after attention blocks.
#' @param dropout Vector of dropout rates after attention block(s). 
#' @param dropout_dense Dropout for dense layers.
#' @param flatten_method How to process output of last attention block. Can be `"max_ch_first"`, `"max_ch_last"`, `"average_ch_first"`,
#' `"average_ch_last"`, `"both_ch_first"`, `"both_ch_last"`, `"all"`, `"none"` or `"flatten"`.
#' If `"average_ch_last"` /  `"max_ch_last"`  or `"average_ch_first"` / `"max_ch_first"`, will apply global average/max pooling.
#' `_ch_first` / `_ch_last` to decide along which axis. `"both_ch_first"` / `"both_ch_last"` to use max and average together. `"all"` to use all 4 
#' global pooling options together. If `"flatten"`, will flatten output after last attention block. If `"none"` no flattening applied.
#' @examples
#' 
#' maxlen <- 50
#' \donttest{
#' library(keras)
#' model <- create_model_transformer(maxlen = maxlen,
#'                                   head_size=c(10,12),
#'                                   num_heads=c(7,8),
#'                                   ff_dim=c(5,9),
#'                                   dropout=c(0.3, 0.5))
#' }
#' @returns A keras model implementing transformer architecture.
#' @export
create_model_transformer <- function(maxlen,
                                     vocabulary_size = 4,
                                     embed_dim = 64,
                                     pos_encoding = "embedding",
                                     head_size = 4L,
                                     num_heads = 5L,
                                     ff_dim = 8,
                                     dropout=0,
                                     n = 10000, # pos emb frequency
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
                                     mirrored_strategy = NULL) {
  
  if (mixed_precision) tensorflow::tf$keras$mixed_precision$set_global_policy("mixed_float16")
  
  if (is.null(mirrored_strategy)) mirrored_strategy <- ifelse(count_gpu() > 1, TRUE, FALSE)
  if (mirrored_strategy) {
    mirrored_strategy <- tensorflow::tf$distribute$MirroredStrategy()
    with(mirrored_strategy$scope(), { 
      argg <- as.list(environment())
      argg$mirrored_strategy <- FALSE
      model <- do.call(create_model_transformer, argg)
    })
    return(model)
  }
  
  stopifnot(length(head_size) == length(num_heads))
  stopifnot(length(head_size) == length(dropout))
  stopifnot(length(head_size) == length(ff_dim))
  stopifnot(flatten_method %in% c("max_ch_first", "max_ch_last", "average_ch_first",
                                  "average_ch_last", "both_ch_first", "both_ch_last", "all", "none", "flatten"))
  stopifnot(pos_encoding %in% c("sinusoid", "embedding"))
  num_dense_layers <- length(layer_dense)
  head_size <- as.integer(head_size)
  num_heads <- as.integer(num_heads)
  maxlen <- as.integer(maxlen)
  num_attention_blocks <- length(num_heads)
  vocabulary_size <-  as.integer(vocabulary_size)
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  
  if (embed_dim == 0) {
    input_tensor <- keras::layer_input(shape = c(maxlen, vocabulary_size))
  } else {
    input_tensor <- keras::layer_input(shape = c(maxlen))
  }
  
  # positional encoding 
  if (pos_encoding == "sinusoid") { 
    pos_enc_layer <- layer_pos_sinusoid_wrapper(maxlen = maxlen, vocabulary_size = vocabulary_size,
                                                n = n, embed_dim = embed_dim)
  } 
  if (pos_encoding == "embedding") { 
    pos_enc_layer <- layer_pos_embedding_wrapper(maxlen = maxlen, vocabulary_size = vocabulary_size,
                                                 embed_dim = embed_dim)
  }
  output_tensor <- input_tensor %>% pos_enc_layer
  
  # attention blocks
  for (i in 1:num_attention_blocks) {
    attn_block <- layer_transformer_block_wrapper(
      num_heads = num_heads[i],
      head_size = head_size[i],
      dropout_rate = dropout[i],
      ff_dim = ff_dim[i],
      embed_dim = embed_dim,
      vocabulary_size = vocabulary_size,
      load_r6 = FALSE)
    output_tensor <- output_tensor %>% attn_block
  }
  
  if (flatten_method != "none") {
    output_tensor <- pooling_flatten(global_pooling = flatten_method, output_tensor = output_tensor)
  }
  
  # dense layers
  if (num_dense_layers > 1) {
    for (i in 1:(num_dense_layers - 1)) {
      output_tensor <- output_tensor %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
      if (!is.null(dropout_dense)) {
        output_tensor <- output_tensor %>% keras::layer_dropout(rate = dropout_dense[i])
      }
    }
  }  
  
  output_tensor <- output_tensor %>% keras::layer_dense(units = layer_dense[length(layer_dense)], activation = last_layer_activation, dtype = "float32")
  
  # create model
  model <- keras::keras_model(inputs = input_tensor, outputs = output_tensor)
  
  model <- compile_model(model = model, label_smoothing = label_smoothing, layer_dense = layer_dense,
                         solver = solver, learning_rate = learning_rate, loss_fn = loss_fn, 
                         num_output_layers = 1, label_noise_matrix = label_noise_matrix,
                         bal_acc = bal_acc, f1_metric = f1_metric, auc_metric = auc_metric)
  
  if (verbose) print(model)
  
  argg <- c(as.list(environment()))
  model <- add_hparam_list(model, argg)
  reticulate::py_set_attr(x = model, name = "hparam", value = model$hparam)
  
  model
  
}
