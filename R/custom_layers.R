#' Aggregation layer 
#' 
#' Aggregate output of time distribution representations using sum, max and/or mean function.
#' 
#' @param load_r6 Whether to load the R6 layer class.
#' @param method At least one of the options, `"sum", "max"` or `"mean"`.
#' @param multi_in Whether to aggregate for a model with multiple inputs (and shared weights).
#' @examples
#' l <- layer_aggregate_time_dist_wrapper() 
#' 
#' @returns A keras layer applying pooling operation(s).
#' @export
layer_aggregate_time_dist_wrapper <- function(load_r6 = FALSE, method = "sum", multi_in = FALSE) {
  
  layer_aggregate_time_dist <- keras::new_layer_class(
    "layer_aggregate_time_dist",
    
    initialize = function(method, ...) {
      super$initialize(...)
      self$method <- method
      self$axis <- ifelse(multi_in, 0L, 1L)
    },
    
    call = function(inputs, mask = NULL) {
      out <- list()
      if ("sum" %in% self$method) {
        out <- c(out, tensorflow::tf$math$reduce_sum(inputs, axis = self$axis))
      }
      if ("mean" %in% self$method) {
        out <- c(out, tensorflow::tf$math$reduce_mean(inputs, axis = self$axis))
      }
      if ("max" %in% self$method) {
        out <- c(out, tensorflow::tf$math$reduce_max(inputs, axis = self$axis))
      }
      
      if (length(out) > 1) {
        out <- tensorflow::tf$concat(out, axis = -1L)
      } else {
        out <- out[[1]]
      }
      
      out
    },
    
    get_config = function() {
      config <- super$get_config()
      config$method <- self$method
      config$multi_in <- self$multi_in
      config
    }
  )
  
  if (load_r6) {
    return(layer_aggregate_time_dist)
  } else {
    return(layer_aggregate_time_dist(method = method))
  }
  
}


#' Layer for positional embedding
#' 
#' Positional encoding layer with learned embedding.
#' 
#' @inheritParams create_model_transformer
#' @param load_r6 Whether to load the R6 layer class.
#' @examples
#' l <- layer_pos_embedding_wrapper()
#' 
#' @returns A keras layer implementing positional embedding.
#' @export
layer_pos_embedding_wrapper <- function(maxlen = 100, vocabulary_size = 4, load_r6 = FALSE, embed_dim = 64) {
  
  layer_pos_embedding <- keras::new_layer_class(
    "layer_pos_embedding",
    
    initialize = function(maxlen=100, vocabulary_size=4, embed_dim=64, ...) {
      super$initialize(...)
      if (embed_dim != 0) {
        self$token_emb <- tensorflow::tf$keras$layers$Embedding(input_dim = as.integer(vocabulary_size),
                                                                output_dim = as.integer(embed_dim))
        self$position_embeddings <- tensorflow::tf$keras$layers$Embedding(input_dim = as.integer(maxlen),
                                                                          output_dim = as.integer(embed_dim))
      } else {
        self$position_embeddings <- tensorflow::tf$keras$layers$Embedding(input_dim = as.integer(maxlen),
                                                                          output_dim = as.integer(vocabulary_size))
      }
      self$embed_dim <- as.integer(embed_dim)
      self$maxlen <- as.integer(maxlen)
      self$vocabulary_size <- as.integer(vocabulary_size)
    },
    
    call = function(inputs) {
      positions <- tensorflow::tf$range(self$maxlen, dtype = "int32") 
      embedded_positions <- self$position_embeddings(positions)
      if (self$embed_dim != 0) inputs <- self$token_emb(inputs)
      inputs + embedded_positions
    },
    
    get_config = function() {
      config <- super$get_config()
      config$maxlen <- self$maxlen
      config$vocabulary_size <- self$vocabulary_size
      config$embed_dim <- self$embed_dim
      config
    }
  )
  
  if (load_r6) {
    return(layer_pos_embedding)
  } else {
    return(layer_pos_embedding(maxlen=maxlen, vocabulary_size=vocabulary_size, embed_dim=embed_dim))
  }
  
}

#' Layer for positional encoding
#' 
#' Positional encoding layer with sine/cosine matrix of different frequencies.
#' 
#' @inheritParams create_model_transformer
#' @param load_r6 Whether to load the R6 layer class.
#' @examples
#' l <- layer_pos_sinusoid_wrapper() 
#' 
#' @returns A keras layer implementing positional encoding using sine/cosine waves.
#' @export
layer_pos_sinusoid_wrapper <- function(maxlen = 100, vocabulary_size = 4, n = 10000, load_r6 = FALSE, embed_dim = 64) {
  
  layer_pos_sinusoid <- keras::new_layer_class(
    "layer_pos_sinusoid",
    initialize = function(maxlen, vocabulary_size, n, embed_dim, ...) {
      super$initialize(...)
      self$maxlen <- as.integer(maxlen)
      self$vocabulary_size <- vocabulary_size
      self$n <- as.integer(n)
      self$pe_matrix <- positional_encoding(seq_len = maxlen,
                                            d_model = ifelse(embed_dim == 0,
                                                             as.integer(vocabulary_size),
                                                             as.integer(embed_dim)),  
                                            n = n)
      
      if (embed_dim != 0) {
        self$token_emb <- tensorflow::tf$keras$layers$Embedding(input_dim = vocabulary_size, output_dim = as.integer(embed_dim))
      }
      self$embed_dim <- as.integer(embed_dim)
      
    },
    
    call = function(inputs) {
      if (self$embed_dim != 0) {
        inputs <- self$token_emb(inputs)
      } 
      inputs + self$pe_matrix
    },
    
    get_config = function() {
      config <- super$get_config()
      config$maxlen <- self$maxlen
      config$vocabulary_size <- self$vocabulary_size
      config$n <- self$n
      config$embed_dim <- self$embed_dim
      config$pe_matrix <- self$pe_matrix
      config
    }
  )
  
  if (load_r6) {
    return(layer_pos_sinusoid)
  } else {
    return(layer_pos_sinusoid(maxlen=maxlen, vocabulary_size=vocabulary_size, n=n,
                              embed_dim = embed_dim))
  }
  
}



#' Transformer block
#' 
#' Create transformer block. Consists of self attention, dense layers, layer normalization, recurrent connection and dropout.
#' 
#' @inheritParams create_model_transformer
#' @param dropout_rate Rate to randomly drop out connections.
#' @param load_r6 Whether to return the layer class.
#' @examples
#' l <- layer_transformer_block_wrapper()
#' 
#' @returns A keras layer implementing a transformer block.
#' @export
layer_transformer_block_wrapper <- function(num_heads = 2, head_size = 4, dropout_rate = 0, ff_dim = 64,  
                                            vocabulary_size = 4, load_r6 = FALSE, embed_dim = 64) {
  
  layer_transformer_block <- keras::new_layer_class(
    "layer_transformer_block",
    initialize = function(num_heads=2, head_size=4, dropout_rate=0, ff_dim=64L, vocabulary_size=4, embed_dim=64, ...) {
      super$initialize(...)
      self$num_heads <- num_heads
      self$head_size <- head_size
      self$dropout_rate <- dropout_rate
      self$ff_dim <- ff_dim
      self$embed_dim <- as.integer(embed_dim)
      self$vocabulary_size <- vocabulary_size
      self$att <- tensorflow::tf$keras$layers$MultiHeadAttention(num_heads=as.integer(num_heads),
                                                                 key_dim=as.integer(head_size))
      
      self$ffn <- keras::keras_model_sequential() %>% keras::layer_dense(units=as.integer(ff_dim), activation="relu") %>%
        keras::layer_dense(units=ifelse(embed_dim == 0, as.integer(vocabulary_size), as.integer(embed_dim)))
      
      self$layernorm1 <- keras::layer_layer_normalization(epsilon=1e-6)
      self$layernorm2 <- keras::layer_layer_normalization(epsilon=1e-6)
      self$dropout1 <- keras::layer_dropout(rate=dropout_rate)
      self$dropout2 <- keras::layer_dropout(rate=dropout_rate)
    },
    
    call = function(inputs) {
      attn_output <- self$att(inputs, inputs, inputs)
      attn_output <- self$dropout1(attn_output)
      out1 <- self$layernorm1(inputs + attn_output)
      ffn_output <- self$ffn(out1)
      ffn_output <- self$dropout2(ffn_output)
      seq_output <- self$layernorm2(out1 + ffn_output)
      return(seq_output)
    },
    
    get_config = function() {
      config <- super$get_config()
      config$num_heads <- self$num_heads
      config$head_size <- self$head_size
      config$dropout_rate <- self$dropout_rate
      config$ff_dim <- self$ff_dim
      config$vocabulary_size <- self$vocabulary_size
      config$embed_dim <- self$embed_dim
      config
    }
  )
  
  if (load_r6) {
    return(layer_transformer_block)
  } else {
    return(layer_transformer_block(num_heads=num_heads,
                                   head_size=head_size,
                                   dropout_rate=dropout_rate,
                                   vocabulary_size=vocabulary_size,
                                   embed_dim=embed_dim,
                                   ff_dim=ff_dim))
  }
  
}
