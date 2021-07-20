#' Stop training after specified time
#' 
#' @export
early_stopping_time_cb <- function(stop_time = NULL) {
  
  early_stopping_time_cb_py_class <- reticulate::PyClass("early_stopping_time_cb",
                                                         inherit = tensorflow::tf$keras$callbacks$Callback,                              
                                                         list(
                                                           
                                                           `__init__` = function(self, stop_time) {
                                                             self$start_time <- Sys.time()
                                                             self$stop_time <- stop_time
                                                             NULL
                                                           },
                                                           
                                                           on_batch_end = function(self, epoch, logs) {
                                                             time_passed <- as.double(difftime(Sys.time(), self$start_time, units = "secs"))
                                                             if (time_passed > self$stop_time) {
                                                               self$model$stop_training <- TRUE
                                                             }
                                                           }
                                                           
                                                         ))
  
  early_stopping_time_cb_py_class(stop_time = stop_time)
  
}

#' early stopping callback
#' 
#' @param early_stopping_time Time in minutes after which to stop training.
#' @param early_stopping_patience Stop training if val_loss does not improve for \code{early_stopping_patience}.
#' @param by_time Whether to use time or patience as metric.
#' 
#' @export
early_stopping_cb <- function(early_stopping_patience, early_stopping_time, by_time = TRUE) {
  
  if (by_time) {
    early_stopping_time_cb(stop_time = early_stopping_time)
  } else {   
    keras::callback_early_stopping(patience = early_stopping_patience)
  }
}  

#' log callback
#' 
#' @export
log_cb <- function(run.name) {
  keras::callback_csv_logger(
    paste0(run.name, "_log.csv"),
    separator = ";",
    append = TRUE)
}  

#' learning rate callback
#' 
#' @inheritParams trainNetwork
#' @export
reduce_lr_cb <- function(patience,
                         cooldown,
                         lr.plateau.factor,
                         monitor = "val_acc") {    
  keras::callback_reduce_lr_on_plateau(
    monitor = "val_acc",
    factor = lr.plateau.factor,
    patience = patience,
    cooldown = cooldown)
}

#' checkpoint callback
#' 
#' @inheritParams trainNetwork
#' @export
checkpoint_cb <- function(filepath_checkpoints,
                          save_weights_only,
                          save_best_only) {
  keras::callback_model_checkpoint(filepath = filepath_checkpoints,
                                   save_weights_only = save_weights_only,
                                   save_best_only = save_best_only,
                                   verbose = 1,
                                   monitor = "val_acc")
  
}                

#' hyperparameter callback
#' 
#' @inheritParams trainNetwork
#' @export
hyper_param_model_outside_cb <- function(tensorboard.log, run.name, wavenet_format, cnn_format, model, vocabulary, path, reverseComplements,
                                         labelVocabulary, maxlen, epochs, max.queue.size, lr.plateau.factor, batch.size,
                                         patience, cooldown, steps.per.epoch, step, randomFiles) {
  
  train_hparams <- list(
    run.name = run.name,
    vocabulary = paste(vocabulary, collapse = ","),
    path = paste(unlist(path), collapse = ", "),
    reverseComplements = paste(reverseComplements), 
    labelVocabulary = paste(labelVocabulary, collapse = ", "),
    epochs = epochs, 
    max.queue.size = max.queue.size,
    lr.plateau.factor = lr.plateau.factor,
    batch.size = batch.size,
    patience = patience, 
    cooldown = cooldown,
    steps.per.epoch = steps.per.epoch,
    step = step,
    randomFiles = randomFiles
  )
  #hparams$update(model$hparam)
  model_hparams <- vector("list")
  for (i in names(model$hparam)) {
    model_hparams[[i]] <- model$hparam[i]
  }
  
  hparams_R <- c(train_hparams, model_hparams)
  
  for (i in 1:length(hparams_R)) {
    if (length(hparams_R[[i]]) > 1) { # length(hparams_R[[i]]) == 0 || 
      hparams_R[[i]] <- paste(hparams_R[[i]], collapse = " ")
    }
  }
  
  hparams <- reticulate::dict(hparams_R)
  hp <- reticulate::import("tensorboard.plugins.hparams.api")
  hp$KerasCallback(file.path(tensorboard.log, run.name), hparams, trial_id = run.name) 
}

#' hyperparameter callback
#' 
#' @inheritParams trainNetwork
#' @export
hyper_param_with_model_cb <- function(default_arguments, model, tensorboard.log, run.name, train_type, model_path, path, validation.split, batch.size,
                                      epochs, max.queue.size, lr.plateau.factor,
                                      patience, cooldown, steps.per.epoch, step, randomFiles, initial_epoch, vocabulary, learning.rate,
                                      shuffleFastaEntries, labelVocabulary, solver, numberOfFiles, reverseComplements, wavenet_format, cnn_format) {
  
  model_hparam <- vector("list")
  model_hparam_names <- vector("list")
  for (i in 1:length(default_arguments)) {
    if (is.null(default_arguments[[i]])) {
      model_hparam[i] <- "NULL"
    } else {
      model_hparam[i] <- default_arguments[i]
    }
  }
  names(model_hparam) <- names(default_arguments)
  # hparam from trainNetwork 
  learning.rate <- keras::k_eval(model$optimizer$lr)
  solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
  
  train_hparam_names <- c("train_type", "model_path", "path", "validation.split", "run.name", "batch.size", "epochs", "max.queue.size", "lr.plateau.factor",
                          "patience", "cooldown", "steps.per.epoch", "step", "randomFiles", "initial_epoch", "vocabulary", "learning.rate",
                          "shuffleFastaEntries", "labelVocabulary", "solver", "numberOfFiles", "reverseComplements", "wavenet_format", "cnn_format")
  train_hparam <- vector("list")
  for (i in 1:length(train_hparam_names)) {
    if (is.null(eval(parse(text=train_hparam_names[i])))) {
      train_hparam[[i]] <- "NULL"
    } else if (length(eval(parse(text=train_hparam_names[i])) > 1)) {
      train_hparam[[i]] <- toString(eval(parse(text=train_hparam_names[i])))
      if (length(train_hparam[[i]]) > 1) {
        train_hparam[[i]] <- paste(train_hparam[[i]], collapse = " ")
      }
    } else {
      train_hparam[[i]] <- eval(parse(text=train_hparam_names[i]))
      if (length(train_hparam[[i]]) > 1) {
        train_hparam[[i]] <- paste(train_hparam[[i]], collapse = " ")
      }
    }
  }
  names(train_hparam) <- train_hparam_names
  hparams_R <- c(train_hparam, model_hparam)
  hparams <- reticulate::dict(hparams_R)
  hp <- reticulate::import("tensorboard.plugins.hparams.api")
  return(hp$KerasCallback(file.path(tensorboard.log, run.name), hparams, trial_id = run.name))
}

#' tensorboard callback
#' 
#' @inheritParams trainNetwork
#' @export
tensorboard_cb <- function(tensorboard.log, run.name) {
  keras::callback_tensorboard(file.path(tensorboard.log, run.name),
                              write_graph = TRUE, 
                              histogram_freq = 1,
                              write_images = TRUE,
                              write_grads = TRUE)
}

#' function arguments callback
#' 
#' @inheritParams trainNetwork
#' @export
function_args_cb <- function(argumentList, tensorboard.log, run.name) {
  
  argAsChar <- as.character(argumentList)
  argText <- vector("character")
  if (length(argumentList$path) > 1) {
    
    argsInQuotes <- c("model_path", "checkpoint_path", "run.name", "solver", "format", "output_format",
                      "tensorboard.log", "fileLog", "train_type", "ambiguous_nuc", "added_label_path", "added_label_names")
  } else {
    argsInQuotes <- c("model_path", "path", "path.val", "checkpoint_path", "run.name", "solver", "output_format",
                      "tensorboard.log", "fileLog", "train_type", "ambiguous_nuc", "format", "added_label_path", "added_label_names")
  }
  argText[1] <- "trainNetwork("
  for (i in 2:(length(argumentList) - 1)) {
    arg <- argAsChar[[i]]
    if (names(argumentList)[i] %in% argsInQuotes) {  
      if (arg == "NULL") {
        argText[i] <- paste0(names(argumentList)[i], " = ", arg, ",")
      } else {
        argText[i] <- paste0(names(argumentList)[i], " = ", '\"', arg, '\"', ",")
      }
    } else {
      argText[i] <- paste0(names(argumentList)[i], " = ", arg, ",")
    }
  }
  i <- length(argumentList)
  if (names(argumentList)[i] %in% argsInQuotes) {
    if (arg == "NULL") {
      argText[i] <- paste0(names(argumentList)[i], " = ", argAsChar[[i]], ")")
    } else {
      argText[i] <- paste0(names(argumentList)[i], " = ", '\"', argAsChar[[i]], '\"', ")")
    }
  } else {
    argText[i] <- paste0(names(argumentList)[i], " = ", argAsChar[[i]], ")")
  }
  
  # write function arguments as text in tensorboard
  trainNetworkArguments <- keras::callback_lambda(
    on_train_begin = function(logs) {
      file.writer <- tensorflow::tf$summary$create_file_writer(file.path(tensorboard.log, run.name))
      file.writer$set_as_default()
      tensorflow::tf$summary$text(name="Arguments",  data = argText, step = 0L)
      file.writer$flush()
    }
  )
  trainNetworkArguments
}

#' tensorboard callback
#' 
#' @inheritParams trainNetwork
#' @export
tensorboard_complete_cb <- function(default_arguments, model, tensorboard.log, run.name, train_type, model_path, path, validation.split, batch.size,
                                    epochs, max.queue.size, lr.plateau.factor, patience, cooldown, steps.per.epoch, step, randomFiles, initial_epoch, vocabulary, learning.rate,
                                    shuffleFastaEntries, labelVocabulary, solver, numberOfFiles, reverseComplements, wavenet_format, cnn_format, create_model_function, vocabulary.size, gen_cb,
                                    argumentList, maxlen, labelGen, labelByFolder, label.vocabulary.size, tb_images = FALSE, stateful, target_middle, num_train_files, fileLog,
                                    proportion_per_file, skip_amb_nuc, max_samples) {
  l <- vector("list")
  
  l[[1]] <- hyper_param_model_outside_cb(tensorboard.log = tensorboard.log, run.name = run.name, wavenet_format = wavenet_format, cnn_format = cnn_format, model = model,
                                         vocabulary = vocabulary, path = path, reverseComplements = reverseComplements, labelVocabulary = labelVocabulary,
                                         maxlen = maxlen, epochs = epochs, max.queue.size = max.queue.size, lr.plateau.factor = lr.plateau.factor,
                                         batch.size = batch.size, patience = patience, cooldown = cooldown, steps.per.epoch = steps.per.epoch,
                                         step = step, randomFiles = randomFiles)
  
  l[[2]] <- tensorboard_cb(tensorboard.log = tensorboard.log, run.name = run.name)
  l[[3]] <- function_args_cb(argumentList = argumentList, tensorboard.log = tensorboard.log, run.name = run.name)
  
  proportion_training_files_cb <- reticulate::PyClass("proportion_training_files_cb",
                                                      inherit = tensorflow::tf$keras$callbacks$Callback,                              
                                                      list(
                                                        
                                                        `__init__` = function(self, num_train_files, fileLog, tensorboard.log, run.name, labelVocabulary, path, train_type, 
                                                                              start_index, proportion_per_file, max_samples, step) {
                                                          self$num_train_files <- num_train_files
                                                          self$fileLog <- fileLog
                                                          self$tensorboard.log <- tensorboard.log
                                                          self$run.name <- run.name
                                                          self$labelVocabulary <- labelVocabulary
                                                          self$path <- path
                                                          self$train_type <- train_type
                                                          self$proportion_per_file <- proportion_per_file
                                                          self$max_samples <-  max_samples
                                                          self$step <- step
                                                          self$start_index <- 1
                                                          self$first_epoch <- TRUE
                                                          self$description <- ""
                                                          NULL
                                                        },
                                                        
                                                        on_epoch_end = function(self, epoch, logs) {
                                                          file.writer <- tensorflow::tf$summary$create_file_writer(file.path(self$tensorboard.log, self$run.name))
                                                          file.writer$set_as_default()
                                                          files_used <- read.csv(self$fileLog, stringsAsFactors = FALSE, header = FALSE)
                                                          if (self$train_type == "label_folder") {
                                                            if (self$first_epoch) {
                                                              if (length(self$step == 1)) self$step <- rep(self$step, length(labelVocabulary))
                                                              if (length(self$proportion_per_file) == 1) {
                                                                self$proportion_per_file <- rep(self$proportion_per_file, length(self$labelVocabulary))
                                                              }
                                                              
                                                              for (i in 1:length(self$labelVocabulary)) {
                                                                if (is.null(self$max_samples)) {
                                                                  self$description[i] <- paste0("Using step size ", self$step[i], " and ",
                                                                                                ifelse(is.null(self$proportion_per_file[i]), 1,
                                                                                                       self$proportion_per_file[i]) * 100, "% per file")
                                                                } else {
                                                                  self$description[i] <- paste0("Using step size ", self$step[i], ", ",
                                                                                                ifelse(is.null(self$proportion_per_file[i]), 1,
                                                                                                       self$proportion_per_file[i]) * 100, "% per file and maximum of ",
                                                                                                self$max_samples, " samples per file")
                                                                }
                                                              }
                                                              self$first_epoch <- FALSE
                                                            }
                                                            
                                                            for (i in 1:length(self$labelVocabulary)) {
                                                              files_of_class <-  sum(stringr::str_detect(
                                                                files_used[ , 1], paste(unlist(self$path[[i]]), collapse = "|")
                                                              ))
                                                              files_percentage <- 100 * files_of_class/self$num_train_files[i]
                                                              tensorflow::tf$summary$scalar(name = paste0("training files seen (%): '",
                                                                                                          self$labelVocabulary[i], "'"), data = files_percentage, step = epoch,
                                                                                            description = self$description[i])
                                                            }
                                                          } else {
                                                            files_percentage <- 100 * nrow(files_used)/self$num_train_files
                                                            if (is.null(self$max_samples)) {
                                                              description <- paste0("Using step size ", step, " and ", 
                                                                                    ifelse(is.null(self$proportion_per_file), 1,
                                                                                           self$proportion_per_file) * 100, "% per file")
                                                            } else {
                                                              description <- paste0("Using step size ", step, ", ",
                                                                                    ifelse(is.null(self$proportion_per_file), 1, 
                                                                                           self$proportion_per_file) * 100, "% per file and maximum of ",
                                                                                    self$max_samples, " samples per file")
                                                              
                                                            } 
                                                            tensorflow::tf$summary$scalar(name = paste("training files seen (%)"), data = files_percentage, step = epoch,
                                                                                          description = description)
                                                          }  
                                                          
                                                          file.writer$flush()
                                                        }
                                                        
                                                      ))
  
  
  l[[4]] <- proportion_training_files_cb(num_train_files = num_train_files, fileLog = fileLog, tensorboard.log = tensorboard.log, run.name = run.name,
                                         labelVocabulary = labelVocabulary, path = path, train_type = train_type, proportion_per_file = proportion_per_file,
                                         max_samples = max_samples, step = step)
  
  return(l)
}

#' Reset states at start/end of validation and whenever file changes.
#'   
#' @export
reset_states_cb <- function(fileLog, fileLogVal, num_files_old = 0, num_files_new = 0) {
  
  reset_states_cb_py_class <- reticulate::PyClass("reset_states_cb",
                                                  inherit = tensorflow::tf$keras$callbacks$Callback,                              
                                                  list(
                                                    
                                                    `__init__` = function(self, fileLog, fileLogVal, num_files_old, num_files_new) {
                                                      self$fileLog <- fileLog
                                                      self$fileLogVal <- fileLogVal
                                                      self$num_files_old <- 0
                                                      self$num_files_new <- 0
                                                      self$num_files_old_val <- 0
                                                      self$num_files_new_val <- 0
                                                      NULL
                                                    },
                                                    
                                                    on_test_begin = function(self, epoch, logs) {
                                                      self$model$reset_states()
                                                    },
                                                    
                                                    on_test_end = function(self, epoch, logs) {
                                                      self$model$reset_states()
                                                    },
                                                    
                                                    on_train_batch_begin = function(self, batch, logs) {
                                                      files_used <- readLines(self$fileLog)
                                                      self$num_files_new <- length(files_used)
                                                      if (self$num_files_new > self$num_files_old) {
                                                        self$model$reset_states()
                                                        self$num_files_old <- self$num_files_new
                                                      }
                                                    },
                                                    
                                                    on_test_batch_begin = function(self, batch, logs) {
                                                      files_used <- readLines(self$fileLogVal)
                                                      self$num_files_new_val <- length(files_used)
                                                      if (self$num_files_new_val > self$num_files_old_val) {
                                                        self$model$reset_states()
                                                        self$num_files_old_val <- self$num_files_new_val
                                                      }
                                                    }
                                                    
                                                  ))
  
  reset_states_cb_py_class(fileLog = fileLog, fileLogVal = fileLogVal, num_files_old = 0, num_files_new = 0)
}

#' Do validation after training finished
#' 
#' @export
validation_after_training_cb <- function(gen.val, validation_steps) {
  
  validation_after_training_cb_py_class <- reticulate::PyClass("validation_after_training_cb",
                                                               inherit = tensorflow::tf$keras$callbacks$Callback,                              
                                                               list(
                                                                 
                                                                 `__init__` = function(self, gen.val, validation_steps) {
                                                                   self$gen.val <- gen.val
                                                                   self$validation_steps <- validation_steps
                                                                   NULL
                                                                 },
                                                                 
                                                                 
                                                                 on_train_end = function(self, logs = list()) {
                                                                   validation_eval <- keras::evaluate_generator(
                                                                     object = self$model,
                                                                     generator = gen.val,
                                                                     steps = self$validation_steps,
                                                                     max_queue_size = 10,
                                                                     workers = 1,
                                                                     callbacks = NULL
                                                                   )
                                                                   self$model$val_loss <- validation_eval[["loss"]]
                                                                   self$model$val_acc <- validation_eval[["acc"]]
                                                                 }
                                                                 
                                                               ))
  
  validation_after_training_cb_py_class(gen.val = gen.val, validation_steps = validation_steps)
  
}

#' confusion matrix callback
#' 
#' @export
conf_matrix_cb <- function(path, tensorboard.log, run.name, confMatLabels) {
  
  conf_matrix_cb_py_class <- reticulate::PyClass("conf_matrix_cb",
                                                 inherit = tensorflow::tf$keras$callbacks$Callback,                              
                                                 list(
                                                   
                                                   `__init__` = function(self, path, tensorboard.log, run.name, confMatLabels, graphics = "png") {
                                                     self$path <- path
                                                     self$tensorboard.log <- tensorboard.log
                                                     self$run.name <- run.name
                                                     self$plot_path <- tempfile(pattern = "", fileext = paste0(".", graphics))
                                                     self$confMatLabels <- confMatLabels
                                                     self$epoch <- 0
                                                     self$train_images <- NULL 
                                                     self$val_images <- NULL
                                                     self$graphics <- graphics
                                                     self$text_size <- NULL
                                                     if (length(confMatLabels) < 8) {
                                                       self$text_size <- (10 - (max(nchar(confMatLabels)) * 0.15)) * (0.95^length(confMatLabels))   
                                                     }
                                                     NULL
                                                   },
                                                   
                                                   on_test_begin = function(self, logs) {
                                                     
                                                     # confusion matrix
                                                     df_true_pred <- read.table(self$path)
                                                     names(df_true_pred) <- c("y_true", "y_pred")
                                                     df_true_pred$y_true <- factor(df_true_pred$y_true, levels = 0:(length(confMatLabels) - 1), labels = confMatLabels)
                                                     df_true_pred$y_pred <- factor(df_true_pred$y_pred, levels = 0:(length(confMatLabels) - 1), labels = confMatLabels)
                                                     cm <- yardstick::conf_mat(df_true_pred, y_true, y_pred, dnn = c("Truth", "Prediction"))
                                                     col_sums <- colSums(cm[[1]])
                                                     for (i in 1:nrow(cm[[1]])) {
                                                       if (col_sums[i] == 0) {
                                                         cm[[1]][ , i] <- 0
                                                       } else {
                                                         cm[[1]][ , i] <- cm[[1]][ , i]/col_sums[i]  
                                                       }
                                                     }
                                                     cm[[1]] <- round(cm[[1]], 2)
                                                     suppressMessages(
                                                       cm_plot <- ggplot2::autoplot(cm, type = "heatmap") +
                                                         ggplot2::scale_fill_gradient(low="#D6EAF8", high = "#2E86C1")  + 
                                                         ggplot2::theme(axis.text.x = 
                                                                          ggplot2::element_text(angle=90,hjust=1, size = self$text_size)) +
                                                         ggplot2::theme(axis.text.y = 
                                                                          ggplot2::element_text(size = self$text_size))
                                                     )
                                                     
                                                     if (length(confMatLabels) > 4) {
                                                       plot_size <- (length(confMatLabels) * 1.3) + 1
                                                     } else {
                                                       plot_size <- length(confMatLabels) * 3
                                                     } 
                                                     
                                                     if (self$graphics == "png") {
                                                       
                                                       suppressMessages(ggplot2::ggsave(filename = self$plot_path, plot = cm_plot, device = "png",
                                                                                        width = plot_size,
                                                                                        height = plot_size,
                                                                                        units = "cm"))
                                                       p_cm <- png::readPNG(self$plot_path)
                                                     } else {
                                                       suppressMessages(ggplot2::ggsave(filename = self$plot_path, plot = cm_plot, device = "jpg",
                                                                                        width = plot_size,
                                                                                        height = plot_size,
                                                                                        units = "cm"))
                                                       p_cm <- jpeg::readJPEG(self$plot_path) 
                                                     }
                                                     
                                                     p_cm <- as.array(p_cm)
                                                     p_cm <- array(p_cm, dim = c(1, dim(p_cm)))
                                                     
                                                     num_images <- 1
                                                     train_images <- array(0, dim = c(num_images, dim(p_cm)[-1]))
                                                     train_images[1, , , ] <- p_cm
                                                     self$train_images <- train_images
                                                     
                                                     write.table(x = NULL, file =  path, col.names = FALSE, row.names = FALSE)
                                                     
                                                   },
                                                   
                                                   on_test_end = function(self, logs) {
                                                     
                                                     # confusion matrix
                                                     df_true_pred <- read.table(self$path)
                                                     names(df_true_pred) <- c("y_true", "y_pred")
                                                     df_true_pred$y_true <- factor(df_true_pred$y_true, levels = 0:(length(confMatLabels) - 1), labels = confMatLabels)
                                                     df_true_pred$y_pred <- factor(df_true_pred$y_pred, levels = 0:(length(confMatLabels) - 1), labels = confMatLabels)
                                                     cm <- yardstick::conf_mat(df_true_pred, y_true, y_pred, dnn = c("Truth", "Prediction"))
                                                     col_sums <- colSums(cm[[1]])
                                                     for (i in 1:nrow(cm[[1]])) {
                                                       if (col_sums[i] == 0) {
                                                         cm[[1]][ , i] <- 0
                                                       } else {
                                                         cm[[1]][ , i] <- cm[[1]][ , i]/col_sums[i]  
                                                       }
                                                     }
                                                     cm[[1]] <- round(cm[[1]], 2)
                                                     suppressMessages(
                                                       cm_plot <- ggplot2::autoplot(cm, type = "heatmap") +
                                                         ggplot2::scale_fill_gradient(low="#D6EAF8", high = "#2E86C1")  + 
                                                         ggplot2::theme(axis.text.x = 
                                                                          ggplot2::element_text(angle=90,hjust=1, size = self$text_size)) +
                                                         ggplot2::theme(axis.text.y = 
                                                                          ggplot2::element_text(size = self$text_size))
                                                     )
                                                     
                                                     if (length(confMatLabels) > 4) {
                                                       plot_size <- (length(confMatLabels) * 1.3) + 1
                                                     } else {
                                                       plot_size <- length(confMatLabels) * 3
                                                     } 
                                                     
                                                     if (self$graphics == "png") {
                                                       
                                                       suppressMessages(ggplot2::ggsave(filename = self$plot_path, plot = cm_plot, device = "png",
                                                                                        width = plot_size,
                                                                                        height = plot_size,
                                                                                        units = "cm"))
                                                       p_cm <- png::readPNG(self$plot_path)
                                                     } else {
                                                       suppressMessages(ggplot2::ggsave(filename = self$plot_path, plot = cm_plot, device = "jpg",
                                                                                        width = plot_size,
                                                                                        height = plot_size,
                                                                                        units = "cm"))
                                                       p_cm <- jpeg::readJPEG(self$plot_path) 
                                                     }
                                                     
                                                     p_cm <- as.array(p_cm)
                                                     p_cm <- array(p_cm, dim = c(1, dim(p_cm)))
                                                     num_images <- 1
                                                     val_images <- array(0, dim = c(num_images, dim(p_cm)[-1]))
                                                     val_images[1, , , ] <- p_cm
                                                     self$val_images <- val_images
                                                     
                                                     write.table(x = NULL, file =  path, col.names = FALSE, row.names = FALSE)
                                                   },
                                                   
                                                   on_epoch_end = function(self, epoch, logs) {
                                                     file.writer <- tensorflow::tf$summary$create_file_writer(file.path(self$tensorboard.log, self$run.name))
                                                     file.writer$set_as_default() 
                                                     tensorflow::tf$summary$image(name = "confusion matrix train", data = self$train_images, step = as.integer(epoch))
                                                     tensorflow::tf$summary$image(name = "confusion matrix validation", data = self$val_images, step = as.integer(epoch))
                                                     file.writer$flush()
                                                   }
                                                 ))
  conf_matrix_cb_py_class(path, tensorboard.log, run.name, confMatLabels)
  
}

#' @param inverted_noise_matrix Inverted matrix of probabilities
#' link: https://github.com/giorgiop/loss-correction/blob/15a79de3c67c31907733392085c333547c2f2b16/loss.py#L16-L21
#' If first label contains 5% wrong labels and second label no noise, then   
#' m <- matrix(c(0.95, 0.05, 0, 1), nrow = 2, byrow = TRUE )  
#' inverted_noise_matrix <- solve(m)
#' noisy_loss <- noisy_loss_wrapper(inverted_noise_matrix)
#' To use as loss, add to compile call keras::compile(loss = noisy_loss, ...)
#' @export
noisy_loss_wrapper <- function(noise_matrix) {
  inverted_noise_matrix <- solve(noise_matrix)
  inverted_noise_matrix <- tensorflow::tf$cast(inverted_noise_matrix, dtype = "float32")
  noisy_loss <- function(y_true, y_pred) {
    y_pred <- y_pred / keras::k_sum(y_pred, axis = -1, keepdims = TRUE)
    y_pred <- keras::k_clip(y_pred, tensorflow::tf$keras$backend$epsilon(), 1.0 - tensorflow::tf$keras$backend$epsilon())
    loss <- -1 * keras::k_sum(keras::k_dot(y_true, inverted_noise_matrix) * keras::k_log(y_pred), axis=-1)
    return(loss)
  } 
  noisy_loss
}

# categorical_crossentropy_plus_conf_mat <- function(csv_path = NULL) {
#   categorical_crossentropy <- function(y_true, y_pred) {
#     if (!is.null(csv_path)) {
#       true <- keras::k_argmax(y_true)
#       pred <- keras::k_argmax(y_pred)
#       if (!is.list(dim(y_true))) {
#         df <- data.frame(as.array(true), as.array(pred))
#         write.table(x = df, file = csv_path, append = TRUE, col.names = FALSE, row.names = FALSE)
#       }
#     }
#     cce <- keras::k_categorical_crossentropy(target = y_true, output = y_pred)
#     cce
#   }
#   categorical_crossentropy
#}
