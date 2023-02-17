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
    if ((len - plen) %% i == 0) {
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
    a <- trainvaldat$X %>% tensorflow::tf$convert_to_tensor()
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
      tail(str_remove(
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
  argg[["trainable_params"]] <- model$count_params()
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
  argg[["trainable_params"]] <- model$count_params()
  argg["label_input"] <- NULL
  argg["label_inputs"] <- NULL
  argg["maxlen_1"] <- NULL
  argg["maxlen_2"] <- NULL
  argg["f1"] <- NULL
  argg["output_tensor"] <- NULL
  argg["output_tensor_1"] <- NULL
  argg["output_tensor_2"] <- NULL

  argg["feature_ext_model"] <- NULL

  argg["number_of_cnn_layers"] <- paste(as.character(argg$number_of_cnn_layers), collapse = " ")
  argg["feature_ext_model"] <- NULL
  argg["pe_matrix"] <- NULL
  argg["position_embedding_layer"] <- NULL 
  
  model$hparam <- argg
  model
}
