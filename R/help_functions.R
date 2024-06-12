#' Stride length calculation
#' 
#' Compute the optimal length for Stride.
#'
#'@param maxlen Length of the input sequence.
#'@param plen Length of a patch.
#' @keywords internal
stridecalc <- function(maxlen, plen) {
  vec <- c()
  for (i in ceiling(plen / 3):(floor(plen / 2) - 1)) {
    if ((maxlen - plen) %% i == 0) {
      vec <- c(vec, i)
    }
  }
  return(vec)
}


#' Number of Patches calculation
#' 
#' Compute the Number of Patches.
#'
#'@param plen Length of a patch.
#'@param maxlen Length of the input sequence.
#'@param stride Stride.
#' @keywords internal
nopatchescalc <- function(plen, maxlen, stride) {
  ((maxlen - plen)/stride) + 1
}

maxlencalc <- function(plen, nopatches, stride) {
  (nopatches - 1) * stride + plen
}


#' Checkpoints saving function
#'
#'@param cp Type of the checkpoint.
#'@param runname Name of the run. Name will be used to identify output from callbacks.
#'@param model A keras model.
#'@param optimizer A keras optimizer.
#'@param history A keras history object.
#' @keywords internal
savechecks <- function(cp, runname, model, optimizer, history) {
  np = import("numpy", convert = F)
  ## define path for saved objects
  modpath <-
    paste("model_results/models", runname, cp , sep = "/")
  ## save model object
  model %>% keras::save_model_hdf5(paste0(modpath, "mod_temp.h5"))
  file.rename(paste0(modpath, "mod_temp.h5"),
              paste0(modpath, "mod.h5"))
  ## save optimizer object
  np$save(
    paste0(modpath, "opt.npy"),
    np$array(keras::backend(FALSE)$batch_get_value(optimizer$weights),
             dtype = "object"),
    allow_pickle = TRUE
  )
  ## save history object
  saveRDS(history, paste0(modpath, "history_temp.rds"))
  file.rename(paste0(modpath, "history_temp.rds"),
              paste0(modpath, "history.rds"))
  ## print when finished
  cat(paste0("---------- New ", cp, " model saved\n"))
}

#' Tensorboard Writer
#' 
#' Writes the loss and the accuracy for a given epoch to the tensorboard.
#'
#' @param writer Name of the tensorboard writer function.
#' @param loss Computed loss for a given epoch.
#' @param acc Computed accracy for a given epoch.
#' @param epoch Epoch, for which the values shall be written to the tensorboard.
#' @keywords internal
TB_loss_acc <- function(writer, loss, acc, epoch) {
  with(writer$as_default(), {
    tensorflow::tf$summary$scalar('epoch_loss',
                                  loss$result(),
                                  step = tensorflow::tf$cast(epoch, "int64"))
    tensorflow::tf$summary$scalar('epoch_accuracy',
                                  acc$result(),
                                  step = tensorflow::tf$cast(epoch, "int64"))
  })
}


#' Step function 
#'
#'@param trainvaldat A data generator.
#'@param model A keras model.
#'@param train_type Either `"cpc"`, `"Self-GenomeNet"`.
#'@param training Boolean. Whether this step is a training step.
#' @keywords internal
modelstep <-
  function(trainvaldat,
           model,
           train_type = "cpc",
           training = F) {
    ## get batch
    a <- trainvaldat$x %>% tensorflow::tf$convert_to_tensor()
    if (train_type == "Self-GenomeNet") {
      ## get complement 
      a_complement <-
        tensorflow::tf$convert_to_tensor(array(as.array(a)[, (dim(a)[2]):1, 4:1], dim = c(dim(a)[1], dim(a)[2], dim(a)[3])))
      a <- tensorflow::tf$concat(list(a, a_complement), axis = 0L)
    }
    ## insert data in model
    model(a, training = training)
  }


#' Reading Pretrained Model function
#'
#'@param pretrained_model The path to a saved keras model.
#' @keywords internal
ReadOpt <- function(pretrained_model) {
  ## Read configuration
  optconf <-
    readRDS(paste(sub("/[^/]+$", "", pretrained_model),
                  "optconfig.rds",
                  sep = "/"))
  ## Read optimizer
  optimizer <- tensorflow::tf$optimizers$Adam$from_config(optconf)
  # Initialize optimizer
  with(
    backend()$name_scope(optimizer$`_name`),
    with(tensorflow::tf$python$framework$ops$init_scope(), {
      optimizer$iterations
      optimizer$`_create_hypers`()
      optimizer$`_create_slots`(model$trainable_weights)
    })
  )
  # Read optimizer weights
  wts2 <-
    np$load(paste(
      sub("/[^/]+$", "", pretrained_model),
      "/",
      utils::tail(stringr::str_remove(
        strsplit(pretrained_model, "/")[[1]], "mod.h5"
      ), 1),
      "opt.npy",
      sep = ""
    ), allow_pickle = TRUE)
  
  # Set optimizer weights
  optimizer$set_weights(wts2)
  return(optimizer)
}

#' Learning Rate Schedule - Parameter Check
#' 
#' Checks, whether all necessary parameters for a defined learning rate schedule are given.
#'
#' @param lr_schedule The name of a learning rate schedule.
#' @keywords internal
LRstop <- function(lr_schedule) {
  # cosine annealing
  if ("cosine_annealing" %in% lr_schedule) {
    if (!isTRUE(all.equal(sort(names(lr_schedule)), sort(
      c("schedule", "lrmin", "lrmax", "restart", "mult")
    )))) {
      stop(
        "Please define lrmin, lrmax, restart, and mult within the list to use cosine annealing"
      )
    }
    # step decay
  } else if ("step_decay" %in% lr_schedule) {
    if (!isTRUE(all.equal(sort(names(lr_schedule)), sort(
      c("schedule", "lrmax", "newstep", "mult")
    )))) {
      stop("Please define lrmax, newstep, and mult within the list to use step decay")
    }
    # exponential decay
  } else if ("exp_decay" %in% lr_schedule) {
    if (!isTRUE(all.equal(sort(names(lr_schedule)), sort(c(
      "schedule", "lrmax", "mult"
    ))))) {
      stop("Please define lrmax, and mult within the list to use exponential decay")
    }
  }
}

#' Learning Rate Calculator
#' 
#' Computes the learning rate for a given epoch.
#'
#'@param lr_schedule The name of a learning rate schedule.
#'@param epoch Epoch, for which the learning rate shall be calculated.
#' @keywords internal
getEpochLR <- function(lr_schedule, epoch) {
  if (lr_schedule$schedule == "cosine_annealing") {
    # cosine annealing
    sgdr(
      lrmin = lr_schedule$lrmin,
      restart = lr_schedule$restart,
      lrmax = lr_schedule$lrmax,
      mult = lr_schedule$mult,
      epoch = epoch
    )
  } else if (lr_schedule$schedule == "step_decay") {
    # step decay
    stepdecay(
      newstep = lr_schedule$newstep,
      lrmax = lr_schedule$lrmax,
      mult = lr_schedule$mult,
      epoch = epoch
    )
    
  } else if (lr_schedule$schedule == "exp_decay") {
    # exponential decay
    exp_decay(
      lrmax = lr_schedule$lrmax,
      mult = lr_schedule$mult,
      epoch = epoch
    )
  }
}


########################################################################################################
########################################### Parameter Lists ############################################
########################################################################################################
#
GenParams <- function(maxlen,
                      batch_size,
                      step,
                      proportion_per_seq,
                      max_samples) {
  checkmate::assertInt(maxlen, lower = 1)
  checkmate::assertInt(batch_size, lower = 1)
  checkmate::assertInt(step, lower = 1)
  checkmate::assertInt(max_samples, lower = 1, null.ok = T)
  checkmate::assertNumber(
    proportion_per_seq,
    lower = 0,
    upper = 1,
    null.ok = T
  )
  
  structure(
    list(
      maxlen = maxlen,
      batch_size = batch_size,
      step = step,
      proportion_per_seq = proportion_per_seq,
      max_samples = max_samples
    ),
    class = "Params"
  )
}


GenTParams <- function(path,
                       shuffle_file_orderTrain,
                       path_file_log,
                       seed) {
  checkmate::assertLogical(shuffle_file_orderTrain)
  checkmate::assertInt(seed)
  
  structure(
    list(
      path_corpus = path,
      shuffle_file_order = shuffle_file_orderTrain,
      path_file_log = path_file_log,
      seed = seed
    ),
    class = "Params"
  )
}

GenVParams <- function(path_val,
                       shuffle_file_orderVal) {
  checkmate::assertLogical(shuffle_file_orderVal)
  
  structure(list(path_corpus = path_val[[1]],
                 shuffle_file_order = shuffle_file_orderVal),
            class = "Params")
}

# add list of hyperparameters to model
add_hparam_list <- function(model, argg) {
  
  argg["model_metrics"] <- NULL
  argg["model"] <- NULL
  argg["i"] <- NULL
  argg["optimizer"] <- NULL
  argg["layer_lstm"] <- paste(as.character(argg$layer_lstm), collapse = " ")
  argg["filters"] <- paste(as.character(argg$filters), collapse = " ")
  argg["kernel_size"] <- paste(as.character(argg$kernel_size), collapse = " ")
  argg["pool_size"] <- paste(as.character(argg$pool_size), collapse = " ")
  argg["strides"] <- paste(as.character(argg$strides), collapse = " ")
  argg["residual_block"] <- paste(as.character(argg$residual_block), collapse = " ")
  argg["residual_block_length"] <- paste(as.character(argg$residual_block_length), collapse = " ")
  argg["size_reduction_1Dconv"] <- paste(as.character(argg$size_reduction_1Dconv), collapse = " ")
  argg["layer_dense"] <- paste(as.character(argg$layer_dense), collapse = " ")
  argg["padding"] <- paste(as.character(argg$padding), collapse = " ")
  argg["use_bias"] <- paste(as.character(argg$use_bias), collapse = " ")
  argg["input_label_list"] <- paste(as.character(argg$layer_dense), collapse = " ")
  argg["num_heads"] <- paste(as.character(argg$num_heads), collapse = " ")
  argg["head_size"] <- paste(as.character(argg$head_size), collapse = " ")
  argg["dropout"] <- paste(as.character(argg$dropout), collapse = " ")
  argg["input_tensor"] <- NULL
  argg["label_inputs"] <- NULL
  argg["f1"] <- NULL
  argg["multi_acc"] <- NULL
  argg[["number_model_params"]] <- model$count_params()
  for (i in 1:length(argg$label_input)) {
    argg[paste0("input_tensor_", i)] <- NULL
    argg[paste0("label_input_layer_", i)] <- NULL
  }
  argg["output_tensor"] <- NULL
  argg["output_list"] <- NULL
  argg["residual_layer"] <- NULL
  argg["label_noise_matrix"] <- NULL
  argg["smooth_loss"] <- NULL
  argg["noisy_loss"] <- NULL
  argg["col_sums"] <- NULL
  argg["auc"] <- NULL
  argg["multi_label"] <- NULL
  argg["macro_average_cb"] <- NULL
  argg["verbose"] <- NULL
  argg["embedded_indices"] <- NULL
  argg["position_indices"] <- NULL
  
  argg["optimizer"] <- NULL
  argg["residual_blocks"] <- paste(as.character(argg$residual_blocks), collapse = " ")
  
  argg["model_metrics"] <- NULL
  argg["i"] <- NULL
  argg["optimizer"] <- NULL
  argg["model"] <- NULL
  argg["input_tensor_1"] <- NULL
  argg["input_tensor_2"] <- NULL
  argg["input_label_list"] <- NULL
  for (i in 1:length(argg$label_input)) {
    argg[paste0("input_tensor_", i)] <- NULL
    argg[paste0("label_input_layer_", i)] <- NULL
  }
  argg[["number_model_params"]] <- model$count_params()
  argg["label_input"] <- NULL
  argg["label_inputs"] <- NULL
  argg["maxlen_1"] <- NULL
  argg["maxlen_2"] <- NULL
  argg["f1"] <- NULL
  argg["output_tensor"] <- NULL
  argg["output_tensor_1"] <- NULL
  argg["output_tensor_2"] <- NULL
  argg[["attn_block"]] <- NULL
  argg["feature_ext_model"] <- NULL
  argg["ff_dim"] <- NULL
  argg["pos_enc_layer"] <- NULL
  argg["number_of_cnn_layers"] <- paste(as.character(argg$number_of_cnn_layers), collapse = " ")
  argg["feature_ext_model"] <- NULL
  argg["pe_matrix"] <- NULL
  argg["position_embedding_layer"] <- NULL 
  argg["layer_add_td"] <- NULL 
  
  model$hparam <- argg
  model
}


get_maxlen <- function(model, set_learning, target_middle, read_data, return_int = FALSE,
                       n_gram = NULL) {
  if (is.null(set_learning)) {
    num_in_layers <- length(model$inputs)
    if (num_in_layers == 1) {
      maxlen <- model$input$shape[[2]]
    } else {
      if (!target_middle & !read_data) {
        maxlen <- model$input[[num_in_layers]]$shape[[2]]
      } else {
        maxlen <- model$inputs[[num_in_layers - 1]]$shape[[2]] + model$inputs[[num_in_layers]]$shape[[2]]
      }
    }
    
    if (!is.null(n_gram)) {
      maxlen <- maxlen + n_gram - 1
    }
    
  } else {
    maxlen <- set_learning$maxlen
  }
  return(maxlen)
}

# combine lists containing x, y and sample weight subsets
reorder_masked_lm_lists <- function(array_lists, include_sw = NULL) {
  
  if (is.null(include_sw)) include_sw <- FALSE
  x <- list()
  y <- list()
  sw <- list()
  for (i in 1:length(array_lists)) {
    x[[i]] <- array_lists[[i]]$x
    y[[i]] <- array_lists[[i]]$y
    if (include_sw) sw[[i]] <- array_lists[[i]]$sample_weight
  }
  x <- abind::abind(x, along = 1)
  y <- abind::abind(y, along = 1)
  if (include_sw) sw <- abind::abind(sw, along = 1)
  if (include_sw) {
    return(list(x=x, y=y, sw=sw))
  } else {
    return(list(x=x, y=y))
  }
  
}

# stack 
create_x_y_tensors_lm <- function(sequence_list, nuc_dist_list, target_middle,
                                  maxlen, vocabulary, ambiguous_nuc,
                                  start_index_list, quality_list, target_len,
                                  coverage_list, use_coverage, max_cov,n_gram,
                                  n_gram_stride, output_format, wavenet_format) {
  
  if (!wavenet_format) {
    
    array_list <- purrr::map(1:length(sequence_list),
                             ~seq_encoding_lm(sequence_list[[.x]], nuc_dist = nuc_dist_list[[.x]], adjust_start_ind = TRUE,
                                              maxlen = maxlen, vocabulary = vocabulary, ambiguous_nuc = ambiguous_nuc,
                                              start_ind =  start_index_list[[.x]], 
                                              quality_vector = quality_list[[.x]], target_len = target_len,
                                              cov_vector = coverage_list[[.x]], use_coverage = use_coverage, max_cov = max_cov,
                                              n_gram = n_gram, n_gram_stride = n_gram_stride, output_format = output_format)
    )
    
    if (!is.list(array_list[[1]][[2]])) {
      if (!target_middle) {
        x <- array_list[[1]][[1]]
        y <- array_list[[1]][[2]]
        if (length(array_list) > 1) {
          for (i in 2:length(array_list)) {
            x <- abind::abind(x, array_list[[i]][[1]], along = 1)
            y <- rbind(y, array_list[[i]][[2]])
          }
        }
        
        # coerce y type to matrix
        if (dim(x)[1] == 1) {
          if (is.null(n_gram)) {
            dim(y) <-  c(1, length(vocabulary))
          } else {
            dim(y) <-  c(1, length(vocabulary)^n_gram)
          }
        }
      } else {
        x_1 <- array_list[[1]][[1]][[1]]
        x_2 <- array_list[[1]][[1]][[2]]
        y <- array_list[[1]][[2]]
        if (length(array_list) > 1) {
          for (i in 2:length(array_list)) {
            x_1 <- abind::abind(x_1, array_list[[i]][[1]][[1]], along = 1)
            x_2 <- abind::abind(x_2, array_list[[i]][[1]][[2]], along = 1)
            y <- rbind(y, array_list[[i]][[2]])
          }
        }
        x <- list(x_1, x_2)
        
        # coerce y type to matrix
        if (dim(x_1)[1] == 1) {
          if (is.null(n_gram)) {
            dim(y) <-  c(1, length(vocabulary))
          } else {
            dim(y) <-  c(1, length(vocabulary)^n_gram)
          }
        }
      }
    } else {
      if (!target_middle) {
        x <- array_list[[1]][[1]]
        y <- array_list[[1]][[2]]
        if (length(array_list) > 1) {
          for (i in 2:length(array_list)) {
            x <- abind::abind(x, array_list[[i]][[1]], along = 1)
            for (j in 1:length(y)) {
              y[[j]] <- rbind(y[[j]], array_list[[i]][[2]][[j]] )
            }
          }
        }
        
        # coerce y type to matrix
        if (dim(x)[1] == 1) {
          for (i in 1:length(y)) {
            if (is.null(n_gram)) {
              dim(y[[i]]) <-  c(1, length(vocabulary))
            } else {
              dim(y[[i]]) <-  c(1, length(vocabulary)^n_gram)
            }
          }
        }
      } else {
        x_1 <- array_list[[1]][[1]][[1]]
        x_2 <- array_list[[1]][[1]][[2]]
        y <- array_list[[1]][[2]]
        if (length(array_list) > 1) {
          for (i in 2:length(array_list)) {
            x_1 <- abind::abind(x_1, array_list[[i]][[1]][[1]], along = 1)
            x_2 <- abind::abind(x_2, array_list[[i]][[1]][[2]], along = 1)
            for (j in 1:length(y)) {
              y[[j]] <- rbind(y[[j]], array_list[[i]][[2]][[j]] )
            }
          }
        }
        x <- list(x_1, x_2)
        
        # coerce y type to matrix
        if (dim(x_1)[1] == 1) {
          for (i in 1:length(y)) {
            if (is.null(n_gram)) {
              dim(y[[i]]) <-  c(1, length(vocabulary))
            } else {
              dim(y[[i]]) <-  c(1, length(vocabulary)^n_gram)
            }
          }
        }
      }
    }
    
    # wavenet format
  } else {
    
    if (target_len > 1) {
      stop("target_len must be 1 when using wavenet_format")
    }
    
    # one hot encode strings collected in sequence_list and connect arrays
    array_list <- purrr::map(1:length(sequence_list),
                             ~seq_encoding_lm(sequence_list[[.x]], ambiguous_nuc = ambiguous_nuc, adjust_start_ind = TRUE,
                                              maxlen = maxlen, vocabulary = vocabulary, nuc_dist = nuc_dist_list[[.x]],
                                              start_ind =  start_index_list[[.x]], 
                                              quality_vector = quality_list[[.x]], n_gram = n_gram,
                                              cov_vector = coverage_list[[.x]], use_coverage = use_coverage, max_cov = max_cov,
                                              output_format = output_format)
    )
    
    x <- array_list[[1]][[1]]
    y <- array_list[[1]][[2]]
    if (length(array_list) > 1) {
      for (i in 2:length(array_list)) {
        x <- abind::abind(x, array_list[[i]][[1]], along = 1)
        y <- abind::abind(y, array_list[[i]][[2]], along = 1)
      }
    }
  }
  return(list(x, y))
  
}

# 
slice_tensor_lm <- function(xy, output_format, target_len, n_gram,
                            n_gram_stride, 
                            # maxlen_n_gram = NULL,
                            # target_len_n_gram = NULL, 
                            total_seq_len, return_int) {
  
  xy_dim <- dim(xy)
  
  if (!is.null(n_gram)) {
    target_len <- floor(target_len/n_gram)
  }
  
  if (output_format == "target_right") {
    x_index <- get_x_index(xy_dim, output_format, target_len)
    if (return_int) {
      x <- xy[ , x_index, drop=FALSE]
      y <- xy[ , -x_index]
    } else {
      x <- xy[ , x_index, , drop=FALSE]
      y <- xy[ , -x_index, ]
    }
  }
  
  if (output_format == "wavenet") {
    
    if (target_len != 1) {
      stop("Target length must be 1 for wavenet model")
    }
    x_index <- 1:(xy_dim[2] - target_len)
    y_index <- 2:dim(xy)[2]
    if (return_int) {
      x <- xy[ , x_index, drop=FALSE]
      y <- xy[ , y_index, drop=FALSE]
    } else {
      x <- xy[ , x_index, , drop=FALSE]
      y <- xy[ , y_index, , drop=FALSE]
    }
    
  }
  
  if (output_format == "target_middle_cnn") {
    
    seq_middle <- ceiling(xy_dim[2]/2)
    y_index <- (1:target_len) + (seq_middle - ceiling(target_len/2))
    if (return_int) {
      x <- xy[ , -y_index, drop=FALSE]
      y <- xy[ , y_index]
    } else {
      x <- xy[ , -y_index, , drop=FALSE]
      y <- xy[ , y_index, ]
    }
    
  }
  
  if (output_format == "target_middle_lstm") {
    
    seq_middle <- ceiling(xy_dim[2]/2)
    y_index <- (1:target_len) + (seq_middle - ceiling(target_len/2))
    
    if (return_int) {
      x1 <- xy[ , 1:(min(y_index) - 1), drop=FALSE]
      # reverse order of x2
      x2 <- xy[ , xy_dim[2] : (max(y_index) + 1), drop=FALSE]
      y <- xy[ , y_index]
    } else {
      x1 <- xy[ , 1:(min(y_index) - 1), , drop=FALSE]
      # reverse order of x2
      x2 <- xy[ ,  xy_dim[2] : (max(y_index) + 1), , drop=FALSE]
      y <- xy[ , y_index, ]
    }
    
    x <- list(x1, x2)
    
  }
  
  if (target_len == 1 & xy_dim[1] == 1 & output_format != "wavenet") {
    y <- matrix(y, nrow = 1)
  }
  
  if (target_len > 1 & xy_dim[1] == 1) {
    y <- array(y, dim = c(1, dim(y)))
  }
  
  return(list(x=x, y=y))
  
}

get_x_index <- function(xy_dim, output_format, target_len) {
  
  if (output_format == "target_right") {
    x_index <- 1:(xy_dim[2] - target_len)
  }
  
  #TODO: subset for other formats with n_gram/stride
  
  return(x_index)
}


add_dim <- function(x) {
  
  if (is.null(dim(x))) {
    return(matrix(x, nrow = 1))
  } else {
    return(array(x, dim = c(1, dim(x))))
  }
  
}

shuffle_batches <- function(x, shuffle_index) {
  
  if (!is.list(x)) {
    dim_len <- length(dim(x))
    x <- shuffle_sample(x, dim_len, shuffle_index)
  } else {
    dim_len <- length(dim(x[[1]]))
    for (i in 1:length(x)) {
      x[[i]] <- shuffle_sample(x[[i]], dim_len, shuffle_index)
    }
  }
  
}  

shuffle_sample <- function(x, dim_len, shuffle_index) {
  
  if (is.null(dim_len) | dim_len == 1) {
    x <- x[shuffle_index]
  }
  
  if (dim_len == 2) {
    x <- x[shuffle_index, ]
  }
  
  if (dim_len == 3) {
    x <- x[shuffle_index, , ]
  }
  
  return(x)
}

count_gpu <- function() {
  
  pd <- tensorflow::tf$config$list_physical_devices()
  count <- 0
  for (i in 1:length(pd)) {
    if (pd[[i]]$device_type == "GPU") count <- count + 1
  }
  return(count)
}

to_time_dist <- function(x, samples_per_target) {
  x_dim <- dim(x)
  x_dim_td <- c(x_dim[1], samples_per_target, x_dim[2]/samples_per_target, x_dim[3])
  x_td <- keras::k_reshape(x, shape = x_dim_td)
  keras::k_eval(x_td)
}


#' Plot confusion matrix
#' 
#' Plot confusion matrix, either with absolute numbers or percentages per column (true labels).
#' 
#' @param cm A confusion matrix
#' @param perc Whether to use absolute numbers or percentages.
#' @param cm_labels Labels correspoding to confusion matrix entries.
#' @param round_dig How to round numbers.
#' @param text_size Size of text annotations.
#' @param highlight_diag Whether to highlight entries in diagonal.
#' @examples
#' cm <- matrix(c(90, 1, 0, 2, 7, 1, 8, 3, 1), nrow = 3, byrow = TRUE)
#' plot_cm(cm, perc = TRUE, cm_labels = paste0('label_', 1:3), text_size = 8)
#' 
#' @export
plot_cm <- function(cm, perc = FALSE, cm_labels, round_dig = 2, text_size = 1, highlight_diag = TRUE) {
  
  if (perc) cm <- cm_perc(cm, round_dig)
  cm <- create_conf_mat_obj(cm, cm_labels)
  
  cm_plot <- ggplot2::autoplot(cm, type = "heatmap") +
    ggplot2::scale_fill_gradient(low="#D6EAF8", high = "#2E86C1")  +
    ggplot2::theme(axis.text.x =
                     ggplot2::element_text(angle=90, hjust=1, size = text_size)) +
    ggplot2::theme(axis.text.y =
                     ggplot2::element_text(size = text_size))
  
  if (highlight_diag) {
    #diagonal_data <- data.frame(x = levels(y_true), y = levels(y_pred))
    diagonal_data <- data.frame(x = cm_labels, y = cm_labels)
    cm_plot <- cm_plot + ggplot2::geom_tile(data = diagonal_data, ggplot2::aes(x = x, y = y),
                                            fill = "red", colour = "white", size = 1,
                                            alpha = 0.000001) 
  }
  
  
  # TODO: add conf mat with ComplexHeatmap for bigger sizes
  cm_plot
  
}



#' Remove checkpoints
#' 
#' Remove all but n 'best' checkpoints, based on some condition. Condition can be 
#' accuracy, loss or epoch number.
#' 
#' @param cp_dir Directory containing checkpoints.
#' @param metric Either `"acc"`, `"loss"` or `"last_ep"`. Condition which checkpoints to keep.
#' @param best_n Number of checkpoints to keep.
#' @param ask_before_remove Whether to show files to keep before deleting rest.
#' @examples
#' model <- create_model_lstm_cnn(layer_lstm = 8)
#' checkpoint_folder <- tempfile()
#' dir.create(checkpoint_folder)
#' keras::save_model_hdf5(model, file.path(checkpoint_folder, 'Ep.007-val_loss11.07-val_acc0.6.hdf5'))
#' keras::save_model_hdf5(model, file.path(checkpoint_folder, 'Ep.019-val_loss8.74-val_acc0.7.hdf5'))
#' keras::save_model_hdf5(model, file.path(checkpoint_folder, 'Ep.025-val_loss0.03-val_acc0.8.hdf5'))
#' remove_checkpoints(cp_dir = checkpoint_folder, metric = "acc", best_n = 2,
#'                    ask_before_remove = FALSE)
#' list.files(checkpoint_folder)
#'  
#' @export 
remove_checkpoints <- function(cp_dir, metric = "acc", best_n = 1, ask_before_remove = TRUE) {
  
  stopifnot(metric %in% c("acc", "loss", "last_ep"))
  stopifnot(dir.exists(cp_dir))
  stopifnot(best_n >= 1)
  
  files <- list.files(cp_dir, full.names = TRUE)
  if (length(files) == 0) {
    stop("Directory is empty")
  }
  files_basename <- basename(files)
  num_cp <- length(files)
  
  if (metric == "acc") {
    if (!all(stringr::str_detect(files_basename, "acc"))) {
      stop("No accuracy information in checkpoint names ('acc' string), use other metric.")
    }
    acc_scores <- files_basename %>% stringr::str_extract("acc\\d++\\.\\d++") %>% 
      stringr::str_remove("acc") %>% as.numeric()
    rank_order <- rank(acc_scores, ties.method = "last")
    index <- rank_order > (num_cp - best_n)
  }
  
  if (metric == "loss") {
    if (!all(stringr::str_detect(files_basename, "loss"))) {
      stop("No loss information in checkpoint names ('loss' string), use other metric.")
    }
    loss_scores <- files_basename %>% stringr::str_extract("loss\\d++\\.\\d++") %>% 
      stringr::str_remove("loss") %>% as.numeric()
    rank_order <- rank(loss_scores, ties.method = "last")
    index <- rank_order <= best_n
  }
  
  if (metric == "last_ep") {
    ep_scores <- files_basename %>% stringr::str_extract("Ep\\.\\d++") %>% 
      stringr::str_remove("Ep\\.") %>% as.numeric()
    rank_order <- rank(ep_scores)
    index <- rank_order > (num_cp - best_n)
  }
  
  if (ask_before_remove) {
    cat("Deleting", sum(!index), paste0(ifelse(sum(!index) == 1, "file", "files") , "."),
        "Only keep \n", paste0(basename(files[index]), collapse = ",\n "), "\n")
    remove_cps <- utils::askYesNo("")
  } else {
    remove_cps <- TRUE
  }
  
  if (is.na(remove_cps)) return(NULL)
  
  if (remove_cps) {
    invisible(file.remove(files[!index]))
  }
  
}

pooling_flatten <- function(global_pooling = NULL, output_tensor) {
  
  if (!is.null(global_pooling)) {
    stopifnot(global_pooling %in% c("max_ch_first", "max_ch_last", "average_ch_first",
                                    "average_ch_last", "both_ch_first", "both_ch_last", "all", "none", "flatten"))
  } else {
    out <- output_tensor %>% keras::layer_flatten()
    return(out)
  }
  
  #if (!is.null(global_pooling) & global_pooling != "flatten") {
  if (stringr::str_detect(global_pooling, "_ch_")) {
    
    if (global_pooling == "max_ch_first") {
      out <- output_tensor %>% keras::layer_global_max_pooling_1d(data_format="channels_first")
    } 
    if (global_pooling == "max_ch_last") {
      out <- output_tensor %>% keras::layer_global_max_pooling_1d(data_format="channels_last")
    } 
    if (global_pooling ==  "average_ch_first") {
      out <- output_tensor %>% keras::layer_global_average_pooling_1d(data_format="channels_first")
    } 
    if (global_pooling ==  "average_ch_last") { 
      out <- output_tensor %>% keras::layer_global_average_pooling_1d(data_format="channels_last")
    } 
    if (global_pooling ==  "both_ch_last") { 
      out1 <- output_tensor %>% keras::layer_global_average_pooling_1d(data_format="channels_last")
      out2 <- output_tensor %>% keras::layer_global_max_pooling_1d(data_format="channels_last")
      out <- keras::layer_concatenate(list(out1, out2))
    } 
    if (global_pooling ==  "both_ch_first") {
      out1 <- output_tensor %>% keras::layer_global_average_pooling_1d(data_format="channels_first")
      out2 <- output_tensor %>% keras::layer_global_max_pooling_1d(data_format="channels_first")
      out <- keras::layer_concatenate(list(out1, out2))
    } 
    if (global_pooling == "all") {
      out1 <- output_tensor %>% keras::layer_global_average_pooling_1d(data_format="channels_first")
      out2 <- output_tensor %>% keras::layer_global_max_pooling_1d(data_format="channels_first")
      out3 <- output_tensor %>% keras::layer_global_average_pooling_1d(data_format="channels_last")
      out4 <- output_tensor %>% keras::layer_global_max_pooling_1d(data_format="channels_last")
      out <- keras::layer_concatenate(list(out1, out2, out3, out4))
    }       
  }
  
  if (global_pooling == "flatten") {
    out <- output_tensor %>% keras::layer_flatten()
  }
  
  if (global_pooling == "none") {
    return(output_tensor)
  }
  
  return(out)
}


pooling_flatten_time_dist <- function(global_pooling = NULL, output_tensor) {
  
  if (!is.null(global_pooling)) {
    stopifnot(global_pooling %in% c("max_ch_first", "max_ch_last", "average_ch_first",
                                    "average_ch_last", "both_ch_first", "both_ch_last", "all", "none", "flatten"))
  } else {
    out <- output_tensor %>% keras::layer_flatten()
    return(out)
  }
  
  #if (!is.null(global_pooling) & global_pooling != "flatten") {
  if (stringr::str_detect(global_pooling, "_ch_")) {
    
    if (global_pooling == "max_ch_first") {
      out <- output_tensor %>% keras::time_distributed(keras::layer_global_max_pooling_1d(data_format="channels_first"))
    } 
    if (global_pooling == "max_ch_last") {
      out <- output_tensor %>% keras::time_distributed(keras::layer_global_max_pooling_1d(data_format="channels_last"))
    } 
    if (global_pooling ==  "average_ch_first") {
      out <- output_tensor %>% keras::time_distributed(keras::layer_global_average_pooling_1d(data_format="channels_first"))
    } 
    if (global_pooling ==  "average_ch_last") { 
      out <- output_tensor %>% keras::time_distributed(keras::layer_global_average_pooling_1d(data_format="channels_last"))
    } 
    if (global_pooling ==  "both_ch_last") { 
      out1 <- output_tensor %>% keras::time_distributed(keras::layer_global_average_pooling_1d(data_format="channels_last"))
      out2 <- output_tensor %>% keras::time_distributed(keras::layer_global_max_pooling_1d(data_format="channels_last"))
      out <- keras::layer_concatenate(list(out1, out2))
    } 
    if (global_pooling ==  "both_ch_first") {
      out1 <- output_tensor %>% keras::time_distributed(keras::layer_global_average_pooling_1d(data_format="channels_first"))
      out2 <- output_tensor %>% keras::time_distributed(keras::layer_global_max_pooling_1d(data_format="channels_first"))
      out <- keras::layer_concatenate(list(out1, out2))
    } 
    if (global_pooling == "all") {
      out1 <- output_tensor %>% keras::time_distributed(keras::layer_global_average_pooling_1d(data_format="channels_first"))
      out2 <- output_tensor %>% keras::time_distributed(keras::layer_global_max_pooling_1d(data_format="channels_first"))
      out3 <- output_tensor %>% keras::time_distributed(keras::layer_global_average_pooling_1d(data_format="channels_last"))
      out4 <- output_tensor %>% keras::time_distributed(keras::layer_global_max_pooling_1d(data_format="channels_last"))
      out <- keras::layer_concatenate(list(out1, out2, out3, out4))
    }       
  }
  
  if (global_pooling == "flatten") {
    out <- output_tensor %>% keras::time_distributed(keras::layer_flatten())
  }
  
  if (global_pooling == "none") {
    return(output_tensor)
  }
  
  return(out)
}



get_pooling_flatten_layer <- function(global_pooling = NULL) {
  
  if (!is.null(global_pooling)) {
    stopifnot(global_pooling %in% c("max_ch_first", "max_ch_last", "average_ch_first",
                                    "average_ch_last", "both_ch_first", "both_ch_last", "all", "none", "flatten"))
  }
  
  if (stringr::str_detect(global_pooling, "_ch_")) {
    
    if (global_pooling == "max_ch_first") {
      out <- keras::layer_global_max_pooling_1d(data_format="channels_first")
    } 
    if (global_pooling == "max_ch_last") {
      out <- keras::layer_global_max_pooling_1d(data_format="channels_last")
    } 
    if (global_pooling ==  "average_ch_first") {
      out <- keras::layer_global_average_pooling_1d(data_format="channels_first")
    } 
    if (global_pooling ==  "average_ch_last") { 
      out <- keras::layer_global_average_pooling_1d(data_format="channels_last")
    } 
    if (global_pooling ==  "both_ch_last") { 
      out1 <- keras::layer_global_average_pooling_1d(data_format="channels_last")
      out2 <- keras::layer_global_max_pooling_1d(data_format="channels_last")
      out <- keras::layer_concatenate(list(out1, out2))
    } 
    if (global_pooling == "both_ch_first") {
      out1 <- keras::layer_global_average_pooling_1d(data_format="channels_first")
      out2 <- keras::layer_global_max_pooling_1d(data_format="channels_first")
      out <- keras::layer_concatenate(list(out1, out2))
    } 
    if (global_pooling == "all") {
      out1 <- keras::layer_global_average_pooling_1d(data_format="channels_first")
      out2 <- keras::layer_global_max_pooling_1d(data_format="channels_first")
      out3 <- keras::layer_global_average_pooling_1d(data_format="channels_last")
      out4 <- keras::layer_global_max_pooling_1d(data_format="channels_last")
      out <- keras::layer_concatenate(list(out1, out2, out3, out4))
    }     
    return(out)
  }
  
  if (is.null(global_pooling) | global_pooling == "flatten") {
    out <- keras::layer_flatten()
    return(out)
  }
  
  return(keras::layer_flatten())
}
