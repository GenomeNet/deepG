get_layer_names <- function(model) {
  
  n <- length(model$layers)
  layer_names <- vector("character", n)
  for (i in 1:n) {
    layer_names[i] <- model$layers[[i]]$name
  }
  layer_names
} 

#' Compile model
#' 
#' @inheritParams create_model_lstm_cnn
#' @param model A keras model.
#' @examplesIf reticulate::py_module_available("tensorflow")
#' model <- create_model_lstm_cnn(layer_lstm = 8, compile = FALSE)
#' model <- compile_model(model = model,
#'                        solver = 'adam',
#'                        learning_rate = 0.01,
#'                        loss_fn = 'categorical_crossentropy')
#' 
#' @returns A compiled keras model.
#' @export
compile_model <- function(model, solver, learning_rate, loss_fn, label_smoothing = 0,
                          num_output_layers = 1, label_noise_matrix = NULL,
                          bal_acc = FALSE, f1_metric = FALSE, auc_metric = FALSE,
                          layer_dense = NULL) {
  
  optimizer <- set_optimizer(solver, learning_rate) 
  if (num_output_layers == 1) num_targets <- model$output_shape[[2]]
  
  #add metrics
  if (loss_fn == "binary_crossentropy") {
    model_metrics <- c(tf$keras$metrics$BinaryAccuracy(name = "acc"))
  } else if (loss_fn == "sparse_categorical_crossentropy") {
    model_metrics <- tensorflow::tf$keras$metrics$SparseCategoricalAccuracy(name = "acc")
  } else {
    model_metrics <- c("acc")
  } 
  
  cm_dir <- NULL
  if (num_output_layers == 1) {
    cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
    while (dir.exists(cm_dir)) {
      cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
    }
    dir.create(cm_dir)
    model$cm_dir <- cm_dir
    
    if (loss_fn == "sparse_categorical_crossentropy") {
      
      if (length(model$outputs) == 1 & length(model$output_shape) == 3) {
        loss_fn <- tensorflow::tf$keras$losses$SparseCategoricalCrossentropy(
          reduction=tensorflow::tf$keras$losses$Reduction$NONE
        )
      }
    }
    
    if (loss_fn == "categorical_crossentropy") {
      
      if (length(model$outputs) == 1 & length(model$output_shape) == 3) {
        loss_fn <- tensorflow::tf$keras$losses$CategoricalCrossentropy(
          label_smoothing=label_smoothing,
          reduction=tensorflow::tf$keras$losses$Reduction$NONE,
          name='categorical_crossentropy'
        )
      }
    }  
    
    if (loss_fn == "categorical_crossentropy" | loss_fn == "sparse_categorical_crossentropy") {
      
      if (bal_acc) {
        macro_average_cb <- balanced_acc_wrapper(num_targets, cm_dir)
        model_metrics <- c(macro_average_cb, "acc")
      }
      
      if (f1_metric) {
        f1 <- f1_wrapper(num_targets)
        model_metrics <- c(model_metrics, f1)
      }
    }
    
    if (auc_metric) {
      auc <- auc_wrapper(model_output_size = layer_dense[length(layer_dense)],
                         loss = loss_fn)
      model_metrics <- c(model_metrics, auc)
    }
    
  }
  
  if (label_smoothing > 0 & !is.null(label_noise_matrix)) {
    stop("Can not use label smoothing and label noise at the same time. Either set label_smoothing = 0 or label_noise_matrix = NULL")
  }
  
  if (label_smoothing > 0) {
    if (loss_fn == "categorical_crossentropy") {
      smooth_loss <- tensorflow::tf$losses$CategoricalCrossentropy(label_smoothing = label_smoothing, name = "smooth_loss")
    }
    if (loss_fn == "binary_crossentropy") {
      smooth_loss <- tensorflow::tf$losses$BinaryCrossentropy(label_smoothing = label_smoothing, name = "smooth_loss")
    }
    model %>% keras::compile(loss = smooth_loss,
                             optimizer = optimizer, metrics = model_metrics)
  } else if (!is.null(label_noise_matrix)) {
    row_sums <- rowSums(label_noise_matrix)
    if (!all(row_sums == 1)) {
      warning("Sum of noise matrix rows don't add up to 1")
    }
    noisy_loss <- noisy_loss_wrapper(solve(label_noise_matrix))
    model %>% keras::compile(loss =  noisy_loss,
                             optimizer = optimizer, metrics = model_metrics)
  } else {
    model %>% keras::compile(loss = loss_fn,
                             optimizer = optimizer, metrics = model_metrics)
  }
  
  model
  
}


#' Load checkpoint 
#' 
#' Load checkpoint from directory. Chooses best checkpoint based on some condition. Condition
#' can be best accuracy, best loss, last epoch number or a specified epoch number.
#' 
#' @inheritParams create_model_lstm_cnn
#' @param cp_path A directory containing checkpoints or a single checkpoint file. 
#' If a directory, choose checkpoint based on `cp_filter` or `ep_index`.
#' @param cp_filter Condition to choose checkpoint if `cp_path` is a directory.
#' Either "acc" for bast validation accuracy, "loss" for best validation loss or "last_ep" for last epoch.
#' @param ep_index Load checkpoint from specific epoch number. If not `NULL`, has priority over `cp_filter`.
#' @param compile Whether to load compiled model.
#' @param re_compile Whether to compile model with parameters from `learning_rate`,
#' `solver` and `loss`.  
#' @param add_custom_object Named list of custom objects.
#' @param verbose Whether to print chosen checkpoint path.
#' @param loss Loss function. Only used if model gets compiled.
#' @param margin Margin for contrastive loss, see \link{loss_cl}.
#' @examplesIf reticulate::py_module_available("tensorflow")
#' model <- create_model_lstm_cnn(layer_lstm = 8)
#' checkpoint_folder <- tempfile()
#' dir.create(checkpoint_folder)
#' keras::save_model_hdf5(model, file.path(checkpoint_folder, 'Ep.007-val_loss11.07-val_acc0.6.hdf5'))
#' keras::save_model_hdf5(model, file.path(checkpoint_folder, 'Ep.019-val_loss8.74-val_acc0.7.hdf5'))
#' keras::save_model_hdf5(model, file.path(checkpoint_folder, 'Ep.025-val_loss0.03-val_acc0.8.hdf5'))
#' model <- load_cp(cp_path = checkpoint_folder, cp_filter = "last_ep")
#' 
#' @returns A keras model loaded from a checkpoint.
#' @export
load_cp <- function(cp_path, cp_filter = "last_ep", ep_index = NULL, compile = FALSE,
                    learning_rate = 0.01, solver = "adam", re_compile = FALSE,
                    loss = "categorical_crossentropy",
                    add_custom_object = NULL, margin = 1,
                    verbose = TRUE, mirrored_strategy = FALSE) {
  
  if (is.null(mirrored_strategy)) mirrored_strategy <- ifelse(count_gpu() > 1, TRUE, FALSE)
  if (mirrored_strategy) {
    mirrored_strategy <- tensorflow::tf$distribute$MirroredStrategy()
    with(mirrored_strategy$scope(), { 
      argg <- as.list(environment())
      argg$mirrored_strategy <- FALSE
      model <- do.call(load_cp, argg)
    })
    return(model)
  }
  
  # custom objects to load transformer architecture
  custom_objects <- list(
    "layer_pos_embedding" = layer_pos_embedding_wrapper(),
    "layer_pos_sinusoid" = layer_pos_sinusoid_wrapper(),
    "layer_transformer_block" = layer_transformer_block_wrapper(),
    "layer_euc_dist" = layer_euc_dist_wrapper(),
    "layer_cosine_sim" = layer_cosine_sim_wrapper(),
    "layer_aggregate_time_dist" = layer_aggregate_time_dist_wrapper(),
    "loss_cl_margin___margin_" = loss_cl(margin=margin)
  )
  
  if (!is.null(add_custom_object)) {
    for (i in 1:length(add_custom_object)) {
      custom_objects[[names(add_custom_object)[i]]] <- add_custom_object[[i]]
    }
  }
  
  cp <- get_cp(cp_path = cp_path, cp_filter = cp_filter,
               ep_index = ep_index, verbose = verbose) 
  
  model <- keras::load_model_hdf5(cp, compile = compile, custom_objects = custom_objects)
  
  if (re_compile) {
    optimizer <- set_optimizer(solver, learning_rate)
    model %>% keras::compile(loss = loss,
                             optimizer = optimizer,
                             metrics = model$metrics)
  }
  
  return(model)
  
}

get_cp  <- function(cp_path, cp_filter = "last_ep", ep_index = NULL, verbose = TRUE) {
  
  if (!is.null(cp_filter)) {
    stopifnot(cp_filter %in% c("acc", "loss", "last_ep"))
    if (!is.null(ep_index)) {
      cp_filter <- NULL
    }
  } 
  
  is_directory <- dir.exists(cp_path)
  if (is_directory) {
    cps <- list.files(cp_path, full.names = TRUE)
    files_basename <- basename(cps)
    stopifnot(xor(is.null(cp_filter), is.null(ep_index)))
  } else {
    stopifnot(file.exists(cp_path))
    cp <- cp_path
  } 
  
  if (is_directory & !is.null(cp_filter)) {
    
    if (cp_filter == "acc") {
      if (!all(stringr::str_detect(files_basename, "acc"))) {
        stop("No accuracy information in checkpoint names ('acc' string), use other metric.")
      }
      acc_scores <- files_basename %>% stringr::str_extract("acc\\d++\\.\\d++") %>% 
        stringr::str_remove("acc") %>% as.numeric()
      # use later epoch for ties
      index <- which.max(rank(acc_scores, ties.method = "last"))
    }
    
    if (cp_filter == "loss") {
      if (!all(stringr::str_detect(files_basename, "loss"))) {
        stop("No loss information in checkpoint names ('loss' string), use other metric.")
      }
      loss_scores <- files_basename %>% stringr::str_extract("loss\\d++\\.\\d++") %>% 
        stringr::str_remove("loss") %>% as.numeric()
      index <- which.min(rank(loss_scores, ties.method = "last"))
    }
    
    if (cp_filter == "last_ep") {
      ep_scores <- files_basename %>% stringr::str_extract("Ep\\.\\d++") %>% 
        stringr::str_remove("Ep\\.") %>% as.numeric()
      index <- which.max(ep_scores)
    }
    
  }
  
  if (is_directory & !is.null(ep_index)) {
    ep_scores <- files_basename %>% stringr::str_extract("Ep\\.\\d++") %>% 
      stringr::str_remove("Ep\\.") %>% as.numeric()
    index <- which(ep_scores == ep_index)
  }
  
  if (is_directory) {
    cp <- cps[index]
  }
  
  if (verbose) {
    cat("Using checkpoint", cp, "\n")
    cat("(Date created:", as.character(file.info(cp)$mtime), ")\n")
  }
  
  return(cp)
  
}


# temporary fix for metric bugs
manage_metrics <- function(model, compile = FALSE) {
  
  dummy_gen <- generator_dummy(model, batch_size = 1)
  z <- dummy_gen()
  x <- z[[1]]
  y <- z[[2]]
  
  if (length(model$metrics) == 0) {
    suppressMessages(
      eval <- model$evaluate(x, y, verbose = 0L)
    )
  }
  
  if (compile) {
    metric_names <- vector("character", length(model$metrics))
    for (i in 1:length(model$metrics)) {
      metric_names[i] <-  model$metrics[[i]]$name
    }
    
    duplicated_index <- duplicated(metric_names)
    loss_index <- stringr::str_detect(metric_names, "loss")
    index <- duplicated_index | loss_index
    
    # remove double metrics
    model <- model %>% keras::compile(loss = model$loss,
                                      optimizer = model$optimizer,
                                      metrics = model$metrics[!index])
    suppressMessages(
      eval <- model$evaluate(x, y, verbose = 0L)
    )
  }
  
  return(model)
  
}


#' Get activation functions of output layers
#' 
#' Get activation functions of output layers.
#' 
#' @param model A keras model.
#' @examplesIf reticulate::py_module_available("tensorflow")
#' model <-  create_model_lstm_cnn(
#'   maxlen = 50,
#'   layer_lstm = 8,
#'   layer_dense = c(64, 2),
#'   verbose = FALSE)
#' get_output_activations(model)
#' 
#' @returns Character vector with names of activation functions of model outputs.
#' @export
get_output_activations <- function(model) {
  
  out_names <- model$output_names
  act_vec <- vector("character", length(out_names))
  count <- 1
  
  for (layer_name in out_names) {
    act_name <- model$get_layer(layer_name)$get_config()$activation
    if (is.null(act_name)) act_name <- "linear"
    act_vec[count] <- act_name
    count <- count + 1
  }
  return(act_vec)
}


set_optimizer <- function(solver = "adam", learning_rate = 0.01) {
  
  stopifnot(solver %in% c("adam", "adagrad", "rmsprop", "sgd"))
  
  named_lr <- "lr" %in% names(formals(keras::optimizer_adam))
  if (named_lr) {
    arg_list <- list(lr = learning_rate)
  } else {
    arg_list <- list(learning_rate = learning_rate)
  }
  
  if (solver == "adam")
    keras_optimizer <- do.call(keras::optimizer_adam, arg_list)
  if (solver == "adagrad")
    keras_optimizer <- do.call(keras::optimizer_adagrad, arg_list)
  if (solver == "rmsprop")
    keras_optimizer <- do.call(keras::optimizer_rmsprop, arg_list)
  if (solver == "sgd")
    keras_optimizer <- do.call(keras::optimizer_sgd, arg_list)
  
  return(keras_optimizer)
  
}


#' Get solver and learning_rate from model.
#'
#' @returns Keras optimizer.
#' @noRd
get_optimizer <- function(model) {
  solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
  learning_rate <- keras::k_eval(model$optimizer$lr)
  if (solver == "adam") {
    optimizer <-  keras::optimizer_adam(learning_rate)
  }
  if (solver == "adagrad") {
    optimizer <- keras::optimizer_adagrad(learning_rate)
  }
  if (solver == "rmsprop") {
    optimizer <- keras::optimizer_rmsprop(learning_rate)
  }
  if (solver == "sgd") {
    optimizer <- keras::optimizer_sgd(learning_rate)
  }
  return(optimizer)
}

#' Replace input layer
#'
#' Replace first layer of model with new input layer of different shape. Only works for sequential models that
#' use CNN and LSTM layers.
#'
#' @param model A keras model.
#' @param input_shape The new input shape vector (without batch size).
#' @examplesIf reticulate::py_module_available("tensorflow")
#' model_1 <-  create_model_lstm_cnn(
#'   maxlen = 50,
#'   kernel_size = c(10, 10),
#'   filters = c(64, 128),
#'   pool_size = c(2, 2),
#'   layer_lstm = c(32),
#'   verbose = FALSE,
#'   layer_dense = c(64, 2))
#' model <- reshape_input(model_1, input_shape = c(120, 4))
#' model
#' 
#' @returns A keras model with changed input shape of input model.
#' @export
reshape_input <- function(model, input_shape) {
  
  in_layer <- layer_input(shape = input_shape)
  for (i in 2:length(model$layers)) {
    layer_name <- model$layers[[i]]$name
    if (i == 2) {
      out_layer <- in_layer %>% model$get_layer(layer_name)()
    } else {
      out_layer <- out_layer %>% model$get_layer(layer_name)()
    }
  }
  new_model <- tensorflow::tf$keras$Model(in_layer, out_layer)
  return(new_model)
}


#' Check if layer is in model
#' 
#' @returns Error message if model does not contain layer of certain name.
#' @noRd
check_layer_name <- function(model, layer_name) {
  num_layers <- length(model$layers)
  layer_names <- vector("character")
  for (i in 1:num_layers) {
    layer_names[i] <- model$layers[[i]]$name
  }
  if (!(layer_name %in% layer_names)) {
    message <- paste0("Model has no layer named ", "'", layer_name, "'")
    stop(message)
  }
}


#' Remove layers from model and add dense layers
#' 
#' Function takes a model as input and removes all layers after a certain layer, specified in \code{layer_name} argument.
#' Optional to add dense layers on top of pruned model. Model can have multiple output layers with separate loss/activation functions.
#' You can freeze all the weights of the pruned model by setting \code{freeze_base_model = TRUE}.
#'
#' @inheritParams create_model_lstm_cnn
#' @param layer_name Name of last layer to use from old model.
#' @param model A keras model. 
#' @param dense_layers List of vectors specifying number of units for each dense layer. If this is a list of length > 1, model
#' has multiple output layers.
#' @param shared_dense_layers Vector with number of units for dense layer. These layers will be connected on top of layer in 
#' argument `layer_name`. Can be used to have shared dense layers, before model has multiple output layers. Don't use if model has just one output layer 
#' (use only `dense_layers`).   
#' @param last_activation List of activations for last entry for each list entry from \code{dense_layers}. Either `"softmax"`, `"sigmoid"` or `"linear"`.
#' @param output_names List of names for each output layer.
#' @param losses List of loss function for each output.
#' @param verbose Boolean.
#' @param dropout List of vectors with dropout rates for each new dense layer.
#' @param dropout_shared Vectors of dropout rates for dense layer from `shared_dense_layers`.
#' @param freeze_base_model Whether to freeze all weights before new dense layers.
#' @param compile Boolean, whether to compile the new model.
#' @param flatten Whether to add flatten layer before new dense layers.
#' @param learning_rate Learning rate if `compile = TRUE`, default learning rate of the old model.
#' @param global_pooling "max_ch_first" for global max pooling with channel first
#' ([keras docs](https://keras.io/api/layers/pooling_layers/global_average_pooling1d/)),
#' "max_ch_last" for global max pooling with channel last, "average_ch_first" for global average pooling with channel first, 
#' "average_ch_last" for global average pooling with channel last or `NULL` for no global pooling. 
#' "both_ch_first" or "both_ch_last" to combine average and max pooling. "all" for all 4 options at once.
#' @examplesIf reticulate::py_module_available("tensorflow")
#' model_1 <- create_model_lstm_cnn(layer_lstm = c(64, 64),
#'                                  maxlen = 50,
#'                                  layer_dense = c(32, 4), 
#'                                  verbose = FALSE)
#' # get name of second to last layer 
#' num_layers <- length(model_1$get_config()$layers)
#' layer_name <- model_1$get_config()$layers[[num_layers-1]]$name
#' # add dense layer with multi outputs and separate loss/activation functions
#' model_2 <- remove_add_layers(model = model_1,
#'                              layer_name = layer_name,
#'                              dense_layers = list(c(32, 16, 1), c(8, 1), c(12, 5)),
#'                              losses = list("binary_crossentropy", "mae",
#'                                            "categorical_crossentropy"),
#'                              last_activation = list("sigmoid", "linear", "softmax"),
#'                              freeze_base_model = TRUE,
#'                              output_names = list("out_1_binary_classsification", 
#'                                                  "out_2_regression", 
#'                                                  "out_3_classification")
#' ) 
#' 
#' @returns A keras model; added and/or removed layers from some base model. 
#' @export
remove_add_layers <- function(model = NULL,
                              layer_name = NULL,
                              dense_layers = NULL,
                              shared_dense_layers = NULL,
                              last_activation = list("softmax"),
                              output_names = NULL,
                              losses = NULL,
                              verbose = TRUE,
                              dropout = NULL,
                              dropout_shared = NULL,
                              freeze_base_model = FALSE,
                              compile = FALSE,
                              learning_rate = 0.001,
                              solver = "adam",
                              flatten = FALSE,
                              global_pooling = NULL,
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
      model <- do.call(remove_add_layers, argg)
    })
    return(model)
  }
  
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  if (!is.null(layer_name)) check_layer_name(model, layer_name)
  if (!is.null(shared_dense_layers) & is.null(dense_layers)) {
    stop("You need to specify output layers in dense_layers argument")
  }
  if (!is.null(shared_dense_layers) & length(dense_layers) == 1) {
    stop("If your model has just one output layer, use only dense_layers argument (and set shared_dense_layers = NULL).")
  }
  if (!is.null(global_pooling)) {
    stopifnot(global_pooling %in% c("max_ch_first", "max_ch_last", "average_ch_first", "average_ch_last", "both_ch_first", "both_ch_last", "all"))
  }
  
  if (!is.list(dense_layers)) {
    dense_layers <- list(dense_layers)
  }
  
  if (!is.null(dropout) && !is.list(dropout)) {
    dropout <- list(dropout)
  }
  
  if (is.null(losses)) {
    losses <- list()
    for (i in 1:length(last_activation)) {
      if (last_activation[[i]] == "softmax") loss <- "categorical_crossentropy"
      if (last_activation[[i]] == "sigmoid") loss <- "binary_crossentropy"
      if (last_activation[[i]] == "linear") loss <- "mse"
      losses[[i]] <- loss
    }
  }
  
  if (is.null(output_names)) {
    output_names <- vector("list", length = length(dense_layers))
  }
  
  if (length(dense_layers) != length(last_activation)) {
    stop("Length of dense_layers and last_activation must be the same")
  }
  
  if (length(dense_layers) != length(output_names)) {
    stop("Length of dense_layers and output_names must be the same")
  }
  
  if (!is.null(dropout)) {
    for (i in 1:length(dense_layers)) {
      stopifnot(length(dropout[[i]]) == length(dense_layers[[i]]))
    }
  }
  
  if (verbose) {
    print("Original model: ")
    print(model$summary())
  }
  
  is_sequential <- any(stringr::str_detect(class(model), "sequential"))
  
  if (!is.null(layer_name)) {
    model_new <- tensorflow::tf$keras$Model(model$input, model$get_layer(layer_name)$output)
    
    if (freeze_base_model) {
      keras::freeze_weights(model_new)
    }
    
    if (!is.null(global_pooling)) {
      if (global_pooling == "max_ch_first") {
        out <- model_new$output %>% keras::layer_global_max_pooling_1d(data_format="channels_first")
      } else if (global_pooling == "max_ch_last") {
        out <- model_new$output %>% keras::layer_global_max_pooling_1d(data_format="channels_last")
      } else if (global_pooling ==  "average_ch_first") {
        out <- model_new$output %>% keras::layer_global_average_pooling_1d(data_format="channels_first")
      } else if (global_pooling ==  "average_ch_last") { 
        out <- model_new$output %>% keras::layer_global_average_pooling_1d(data_format="channels_last")
      } else if (global_pooling ==  "both_ch_last") { 
        out1 <- model_new$output %>% keras::layer_global_average_pooling_1d(data_format="channels_last")
        out2 <- model_new$output %>% keras::layer_global_max_pooling_1d(data_format="channels_last")
        out <- keras::layer_concatenate(list(out1, out2))
      } else if (global_pooling ==  "both_ch_first") {
        out1 <- model_new$output %>% keras::layer_global_average_pooling_1d(data_format="channels_first")
        out2 <- model_new$output %>% keras::layer_global_max_pooling_1d(data_format="channels_first")
        out <- keras::layer_concatenate(list(out1, out2))
      } else {
        out1 <- model_new$output %>% keras::layer_global_average_pooling_1d(data_format="channels_first")
        out2 <- model_new$output %>% keras::layer_global_max_pooling_1d(data_format="channels_first")
        out3 <- model_new$output %>% keras::layer_global_average_pooling_1d(data_format="channels_last")
        out4 <- model_new$output %>% keras::layer_global_max_pooling_1d(data_format="channels_last")
        out <- keras::layer_concatenate(list(out1, out2, out3, out4))
      }       
      model_new <- tensorflow::tf$keras$Model(model_new$input, out)
    }
    
    if (flatten) {
      out <- model_new$output %>% keras::layer_flatten()
      model_new <- tensorflow::tf$keras$Model(model_new$input, out)
    }
    
    if (!is.null(shared_dense_layers)) {
      out <- model_new$output 
      for (i in 1:length(shared_dense_layers)) {
        if (!is.null(dropout_shared)) {
          out <- out %>% keras::layer_dropout(dropout_shared[i])
        }
        out <- out %>% keras::layer_dense(shared_dense_layers[i], activation = "relu")
      }
      model_new <- tensorflow::tf$keras$Model(model_new$input, out)
    }
    
    output_list <- list()
    name_dense_index <- 1
    
    if (!is.null(dense_layers[[1]])) {
      for (output_num in 1:length(dense_layers)) {
        for (i in 1:length(dense_layers[[output_num]])) {
          if (i == length(dense_layers[[output_num]])) {
            activation <- last_activation[[output_num]]
            dtype <- "float32"
          } else {
            activation <- "relu"
            dtype <- NULL
          }
          
          if (i == length(dense_layers[[output_num]])) {
            layer_name <- output_names[[output_num]]
          } else {
            layer_name <- paste0("dense_new_", name_dense_index)
            name_dense_index <- name_dense_index + 1
          }
          
          if (is.null(dropout)) {
            if (i == 1) {
              output_list[[output_num]] <- model_new$output %>%
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name, dtype = dtype)
            } else {
              output_list[[output_num]] <- output_list[[output_num]] %>%
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name, dtype = dtype)
            }
          } else {
            if (i == 1) {
              output_list[[output_num]] <- model_new$output %>%
                keras::layer_dropout(rate = dropout[[output_num]][i]) %>%
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name, dtype = dtype)
            } else {
              output_list[[output_num]] <- output_list[[output_num]] %>%
                keras::layer_dropout(rate = dropout[[output_num]][i]) %>%
                keras::layer_dense(units = dense_layers[[output_num]][i], activation = activation,
                                   name = layer_name, dtype = dtype)
            }
          }
        }
      }
      model_new <- tensorflow::tf$keras$Model(model_new$input, output_list)
    }
  }
  
  if (verbose) {
    print("New model: ")
    print(model_new$summary())
    for (i in 1:length(model_new$layers)) {
      cat(model_new$layers[[i]]$name , "trainable:" , model_new$layers[[i]]$trainable, "\n")
    }
  }
  
  if (compile) {
    if (is.null(learning_rate)) {
      learning_rate <- keras::k_eval(model$optimizer$lr)
    } else {
      learning_rate <- learning_rate
    }
    
    if (is.null(solver)) {
      solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
    } 
    
    optimizer <- set_optimizer(solver, learning_rate) 
    
    metric_list <- list()
    for (i in 1:length(losses)) {
      metric_list[[i]] <- ifelse(losses[[i]] == "binary_crossentropy", "binary_accuracy", "acc")
    }
    
    model_new %>% keras::compile(loss = losses,
                                 optimizer = optimizer,
                                 metrics = metric_list)
  }
  
  model_new
}


#' Merge two models
#' 
#' Combine two models at certain layers and add dense layer(s) afterwards.
#'
#' @param models List of two models.
#' @param layer_names Vector of length 2 with names of layers to merge.
#' @param freeze_base_model Boolean vector of length 2. Whether to freeze weights of individual models.
#' @inheritParams create_model_lstm_cnn
#' @examplesIf reticulate::py_module_available("tensorflow")
#' model_1 <- create_model_lstm_cnn(layer_lstm = c(64, 64), maxlen = 50, layer_dense = c(32, 4),
#'                                  verbose = FALSE)
#' model_2 <- create_model_lstm_cnn(layer_lstm = c(32), maxlen = 40, 
#'                                  layer_dense = c(8, 2), verbose = FALSE)
#' # get names of second to last layers
#' num_layers_1 <- length(model_1$get_config()$layers)
#' layer_name_1 <- model_1$get_config()$layers[[num_layers_1 - 1]]$name
#' num_layers_2 <- length(model_2$get_config()$layers)
#' layer_name_2 <- model_2$get_config()$layers[[num_layers_2 - 1]]$name
#' # merge models
#' model <- merge_models(models = list(model_1, model_2),
#'                       layer_names = c(layer_name_1, layer_name_2),
#'                       layer_dense = c(6, 2), 
#'                       freeze_base_model = c(FALSE, FALSE)) 
#' 
#' @returns A keras model merging two input models.                        
#' @export
merge_models <- function(models, layer_names, layer_dense, solver = "adam",
                         learning_rate = 0.0001,
                         freeze_base_model = c(FALSE, FALSE),
                         model_seed = NULL) {
  
  if (!is.null(model_seed)) tensorflow::tf$random$set_seed(model_seed)
  
  model_1 <- remove_add_layers(model = models[[1]],
                               layer_name = layer_names[1],
                               dense_layers = NULL,
                               verbose = FALSE,
                               dropout = NULL,
                               freeze_base_model = freeze_base_model[1],
                               compile = FALSE,
                               learning_rate = NULL)
  
  model_2 <- remove_add_layers(model = models[[2]],
                               layer_name = layer_names[2],
                               dense_layers = NULL,
                               verbose = FALSE,
                               dropout = NULL,
                               freeze_base_model = freeze_base_model[2],
                               compile = FALSE,
                               learning_rate = NULL)
  
  # choose optimization method
  optimizer <- set_optimizer(solver, learning_rate) 
  
  output_tensor <- keras::layer_concatenate(c(model_1$output, model_2$output))
  num_targets <- layer_dense[length(layer_dense)]
  
  if (length(layer_dense) > 1) {
    for (i in 1:(length(layer_dense) - 1)) {
      output_tensor <- output_tensor %>% keras::layer_dense(units = layer_dense[i], activation = "relu")
    }
  }
  
  output_tensor <- output_tensor %>%
    keras::layer_dense(units = num_targets, activation = "softmax")
  
  model <- keras::keras_model(inputs = list(model_1$input, model_2$input), outputs = output_tensor)
  model %>% keras::compile(loss = "categorical_crossentropy",
                           optimizer = optimizer, metrics = c("acc"))
  model
}


#' Extract hyperparameters from model
#'
#' @param model A keras model.
#' @returns List of hyperparameters.
#' @noRd
get_hyper_param <- function(model) {
  layers.lstm <- 0
  use.cudnn <- FALSE
  bidirectional <- FALSE
  use.codon.cnn <- FALSE
  learning_rate <- keras::k_eval(model$optimizer$lr)
  solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
  
  layerList <- keras::get_config(model)["layers"]
  for (i in 1:length(layerList)) {
    layer_class_name <- layerList[[i]]$class_name
    
    if (layer_class_name == "Conv1D") {
      use.codon.cnn <- TRUE
    }
    
    if (layer_class_name == "MaxPooling1D") {
    }
    
    if (layer_class_name == "BatchNormalization") {
    }
    
    if (layer_class_name == "CuDNNLSTM") {
      layers.lstm <- layers.lstm + 1
      use.cudnn <- TRUE
      lstm_layer_size <- layerList[[i]]$config$units
      recurrent_dropout_lstm <- 0
      dropout_lstm <- 0
    }
    
    if (layer_class_name == "LSTM") {
      layers.lstm <- layers.lstm + 1
      lstm_layer_size <- layerList[[i]]$config$units
      recurrent_dropout_lstm <- layerList[[i]]$config$recurrent_dropout_lstm
      dropout_lstm <- layerList[[i]]$config$dropout
    }
    # TODO: wrong output since bidirectional is layer wrapper (?)
    if (layer_class_name == "Bidirectional") {
      bidirectional <- TRUE
      if (layerList[[i]]$config$layer$class_name == "LSTM") {
        use.cudnn <- FALSE
        layers.lstm <- layers.lstm + 1
        lstm_layer_size <- layerList[[i]]$config$layer$config$units
        recurrent_dropout_lstm <- layerList[[i]]$config$layer$config$recurrent_dropout
        dropout_lstm <- layerList[[i]]$config$layer$config$dropout
      } else {
        use.cudnn <- FALSE
        layers.lstm <- layers.lstm + 1
        lstm_layer_size <- layerList[[i]]$config$layer$config$units
      }
    }
    
    if (layer_class_name == "Dense") {
    }
    
    if (layer_class_name == "Activation") {
    }
  }
  
  list(dropout = dropout_lstm,
       recurrent_dropout = recurrent_dropout_lstm,
       lstm_layer_size =  lstm_layer_size,
       solver = solver,
       use.cudnn = use.cudnn,
       layers.lstm = layers.lstm,
       learning_rate = learning_rate,
       use.codon.cnn = use.codon.cnn,
       bidirectional = bidirectional
  )
}
