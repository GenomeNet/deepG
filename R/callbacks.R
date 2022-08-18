#' Stop training callback
#' 
#' Stop training after specified time.
#'
#' @param stop_time Time in seconds after which to stop training.
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

#' Early stopping callback
#'
#' @param early_stopping_time Time in seconds after which to stop training.
#' @param early_stopping_patience Stop training if val_loss does not improve for \code{early_stopping_patience}.
#' @param by_time Whether to use time or patience as metric.
#' @keywords internal
early_stopping_cb <- function(early_stopping_patience, early_stopping_time, by_time = TRUE) {
  
  if (by_time) {
    early_stopping_time_cb(stop_time = early_stopping_time)
  } else {
    keras::callback_early_stopping(patience = early_stopping_patience)
  }
}

#' Log callback
#'
#' @param path_log Path to output directory.
#' @param run_name Name of output file is run_name + ".csv".
#' @keywords internal
log_cb <- function(path_log, run_name) {
  keras::callback_csv_logger(
    paste0(path_log, "/", run_name, ".csv"),
    separator = ";",
    append = TRUE)
}

#' Learning_rate callback
#'
#' @inheritParams train_model
#' @keywords internal
reduce_lr_cb <- function(patience,
                         cooldown,
                         lr_plateau_factor,
                         monitor = "val_acc") {
  keras::callback_reduce_lr_on_plateau(
    monitor = "val_acc",
    factor = lr_plateau_factor,
    patience = patience,
    cooldown = cooldown)
}

#' Checkpoint callback
#'
#' @inheritParams train_model
#' @keywords internal
checkpoint_cb <- function(filepath_checkpoints,
                          save_weights_only,
                          save_best_only) {
  keras::callback_model_checkpoint(filepath = filepath_checkpoints,
                                   save_weights_only = save_weights_only,
                                   save_best_only = save_best_only,
                                   verbose = 1,
                                   monitor = "val_acc")
  
}

#' Non model hyperparameter callback
#' 
#' Get hyperparameters excluding model parameters.
#'
#' @inheritParams train_model
#' @keywords internal
hyper_param_model_outside_cb <- function(path_tensorboard, run_name, wavenet_format, cnn_format, model, vocabulary, path, reverse_complement,
                                         vocabulary_label, maxlen, epochs, max_queue_size, lr_plateau_factor, batch_size,
                                         patience, cooldown, steps_per_epoch, step, shuffle_file_order) {
  
  train_hparams <- list(
    run_name = run_name,
    vocabulary = paste(vocabulary, collapse = ","),
    path = paste(unlist(path), collapse = ", "),
    reverse_complement = paste(reverse_complement),
    vocabulary_label = paste(vocabulary_label, collapse = ", "),
    epochs = epochs,
    max_queue_size = max_queue_size,
    lr_plateau_factor = lr_plateau_factor,
    batch_size = batch_size,
    patience = patience,
    cooldown = cooldown,
    steps_per_epoch = steps_per_epoch,
    step = step,
    shuffle_file_order = shuffle_file_order
  )
  #hparams$update(model$hparam)
  model_hparams <- vector("list")
  for (i in names(model$hparam)) {
    model_hparams[[i]] <- model$hparam[i]
  }
  
  hparams_R <- c(train_hparams, model_hparams)
  
  keep_entry_index <- rep(TRUE, length(hparams_R))
  for (i in 1:length(hparams_R)) {
    
    if (length(hparams_R[[i]]) == 0) { 
      keep_entry_index[i] <- FALSE
    }
    
    if (length(hparams_R[[i]]) > 1) { 
      hparams_R[[i]] <- paste(hparams_R[[i]], collapse = " ")
    }
  }
  hparams_R <- hparams_R[keep_entry_index]
  
  hparams <- reticulate::dict(hparams_R)
  hp <- reticulate::import("tensorboard.plugins.hparams.api")
  hp$KerasCallback(file.path(path_tensorboard, run_name), hparams, trial_id = run_name)
}

#' Model hyperparameter callback
#' 
#' Get model hyperparameters.
#'
#' @inheritParams train_model
#' @keywords internal
hyper_param_with_model_cb <- function(default_arguments, model, path_tensorboard, run_name, train_type, path_model, path, train_val_ratio, batch_size,
                                      epochs, max_queue_size, lr_plateau_factor,
                                      patience, cooldown, steps_per_epoch, step, shuffle_file_order, initial_epoch, vocabulary, learning_rate,
                                      shuffle_input, vocabulary_label, solver, file_limit, reverse_complement, wavenet_format, cnn_format) {
  
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
  # hparam from train_model
  learning_rate <- keras::k_eval(model$optimizer$lr)
  solver <- stringr::str_to_lower(model$optimizer$get_config()["name"])
  
  train_hparam_names <- c("train_type", "path_model", "path", "train_val_ratio", "run_name", "batch_size", "epochs", "max_queue_size", "lr_plateau_factor",
                          "patience", "cooldown", "steps_per_epoch", "step", "shuffle_file_order", "initial_epoch", "vocabulary", "learning_rate",
                          "shuffle_input", "vocabulary_label", "solver", "file_limit", "reverse_complement", "wavenet_format", "cnn_format")
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
  return(hp$KerasCallback(file.path(path_tensorboard, run_name), hparams, trial_id = run_name))
}

#' Tensorboard callback
#'
#' @inheritParams train_model
#' @keywords internal
tensorboard_cb <- function(path_tensorboard, run_name) {
  keras::callback_tensorboard(file.path(path_tensorboard, run_name),
                              write_graph = TRUE,
                              histogram_freq = 1,
                              write_images = TRUE,
                              write_grads = TRUE)
}

#' Function arguments callback
#' 
#' Print train_model call in text field of tensorboard.
#' 
#' @inheritParams train_model
#' @export
function_args_cb <- function(argumentList, path_tensorboard, run_name) {
  
  argAsChar <- as.character(argumentList)
  argText <- vector("character")
  if (length(argumentList$path) > 1) {
    
    argsInQuotes <- c("path_model", "path_checkpoint", "run_name", "solver", "format", "output_format",
                      "path_tensorboard", "path_file_log", "train_type", "ambiguous_nuc", "added_label_path", "added_label_names",
                      "train_val_split_csv", "target_from_csv")
  } else {
    argsInQuotes <- c("path_model", "path", "path_val", "path_checkpoint", "run_name", "solver", "output_format",
                      "path_tensorboard", "path_file_log", "train_type", "ambiguous_nuc", "format", "added_label_path", "added_label_names",
                      "train_val_split_csv", "target_from_csv")
  }
  argText[1] <- "train_model("
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
  trainArguments <- keras::callback_lambda(
    on_train_begin = function(logs) {
      file.writer <- tensorflow::tf$summary$create_file_writer(file.path(path_tensorboard, run_name))
      file.writer$set_as_default()
      tensorflow::tf$summary$text(name="Arguments",  data = argText, step = 0L)
      file.writer$flush()
    }
  )
  trainArguments
}

#' Tensorboard callback wrapper
#'
#' @inheritParams train_model
#' @keywords internal
tensorboard_complete_cb <- function(default_arguments, model, path_tensorboard, run_name, train_type, path_model, path, train_val_ratio, batch_size,
                                    epochs, max_queue_size, lr_plateau_factor, patience, cooldown, steps_per_epoch, step, shuffle_file_order, initial_epoch, vocabulary, learning_rate,
                                    shuffle_input, vocabulary_label, solver, file_limit, reverse_complement, wavenet_format, cnn_format, create_model_function, vocabulary_size, gen_cb,
                                    argumentList, maxlen, labelGen, labelByFolder, vocabulary_label_size, tb_images = FALSE, stateful, target_middle, num_train_files, path_file_log,
                                    proportion_per_seq, skip_amb_nuc, max_samples, proportion_entries, train_with_gen, count_files = TRUE) {
  l <- vector("list")
  
  l[[1]] <- hyper_param_model_outside_cb(path_tensorboard = path_tensorboard, run_name = run_name, wavenet_format = wavenet_format, cnn_format = cnn_format, model = model,
                                         vocabulary = vocabulary, path = path, reverse_complement = reverse_complement, vocabulary_label = vocabulary_label,
                                         maxlen = maxlen, epochs = epochs, max_queue_size = max_queue_size, lr_plateau_factor = lr_plateau_factor,
                                         batch_size = batch_size, patience = patience, cooldown = cooldown, steps_per_epoch = steps_per_epoch,
                                         step = step, shuffle_file_order = shuffle_file_order)
  
  l[[2]] <- tensorboard_cb(path_tensorboard = path_tensorboard, run_name = run_name)
  l[[3]] <- function_args_cb(argumentList = argumentList, path_tensorboard = path_tensorboard, run_name = run_name)
  
  if (train_with_gen & count_files & train_type != "dummy_gen") {
    
    proportion_training_files_cb <- reticulate::PyClass("proportion_training_files_cb",
                                                        inherit = tensorflow::tf$keras$callbacks$Callback,
                                                        list(
                                                          
                                                          `__init__` = function(self, num_train_files, path_file_log, path_tensorboard, run_name, vocabulary_label,
                                                                                path, train_type, start_index, proportion_per_seq, max_samples, step,
                                                                                proportion_entries) {
                                                            self$num_train_files <- num_train_files
                                                            self$path_file_log <- path_file_log
                                                            self$path_tensorboard <- path_tensorboard
                                                            self$run_name <- run_name
                                                            self$vocabulary_label <- vocabulary_label
                                                            self$path <- path
                                                            self$train_type <- train_type
                                                            self$proportion_per_seq <- proportion_per_seq
                                                            self$max_samples <- max_samples
                                                            self$step <- step
                                                            self$start_index <- 1
                                                            self$first_epoch <- TRUE
                                                            self$description <- ""
                                                            self$proportion_entries <- proportion_entries
                                                            NULL
                                                          },
                                                          
                                                          on_epoch_end = function(self, epoch, logs) {
                                                            if (is.null(self$proportion_entries)) self$proportion_entries <- 1
                                                            file.writer <- tensorflow::tf$summary$create_file_writer(file.path(self$path_tensorboard, self$run_name))
                                                            file.writer$set_as_default()
                                                            files_used <- read.csv(self$path_file_log, stringsAsFactors = FALSE, header = FALSE)
                                                            if (self$train_type == "label_folder") {
                                                              if (self$first_epoch) {
                                                                if (length(self$step == 1)) self$step <- rep(self$step, length(vocabulary_label))
                                                                if (length(self$proportion_per_seq) == 1) {
                                                                  self$proportion_per_seq <- rep(self$proportion_per_seq, length(self$vocabulary_label))
                                                                }
                                                                if (length(max_samples) == 1) self$max_samples <- rep(max_samples, length(vocabulary_label))
                                                                
                                                                for (i in 1:length(self$vocabulary_label)) {
                                                                  if (is.null(self$max_samples)) {
                                                                    self$description[i] <- paste0("Using step size ", self$step[i], ", proportion_entries ",
                                                                                                  self$proportion_entries * 100, "% and ",
                                                                                                  ifelse(is.null(self$proportion_per_seq[i]), 1,
                                                                                                         self$proportion_per_seq[i]) * 100, "% per sequence")
                                                                  } else {
                                                                    self$description[i] <- paste0("Using step size ", self$step[i], ", ",
                                                                                                  ifelse(is.null(self$proportion_per_seq[i]), 1,
                                                                                                         self$proportion_per_seq[i]) * 100, "% per sequence, maximum of ",
                                                                                                  self$max_samples[i], " samples per file and proportion_entries ",
                                                                                                  self$proportion_entries * 100, "%")
                                                                  }
                                                                }
                                                                self$first_epoch <- FALSE
                                                              }
                                                              
                                                              for (i in 1:length(self$vocabulary_label)) {
                                                                files_of_class <-  sum(stringr::str_detect(
                                                                  files_used[ , 1], paste(unlist(self$path[[i]]), collapse = "|")
                                                                ))
                                                                files_percentage <- 100 * files_of_class/self$num_train_files[i]
                                                                tensorflow::tf$summary$scalar(name = paste0("training files seen (%): '",
                                                                                                            self$vocabulary_label[i], "'"), data = files_percentage, step = epoch,
                                                                                              description = self$description[i])
                                                              }
                                                            } else {
                                                              files_percentage <- 100 * nrow(files_used)/self$num_train_files
                                                              if (is.null(self$max_samples)) {
                                                                description <- paste0("Using step size ", step,
                                                                                      ", proportion_entries ", self$proportion_entries * 100, "% and ",
                                                                                      ifelse(is.null(self$proportion_per_seq), 1,
                                                                                             self$proportion_per_seq) * 100, "% per sequence")
                                                              } else {
                                                                description <- paste0("Using step size ", step, ", ",
                                                                                      ifelse(is.null(self$proportion_per_seq), 1,
                                                                                             self$proportion_per_seq) * 100, "% per sequence, maximum of ",
                                                                                      self$max_samples, " samples per file and proportion_entries ",
                                                                                      self$proportion_entries * 100, "%")
                                                                
                                                              }
                                                              if (self$train_type == "label_rds") {
                                                                description <- paste0("Using step size ",
                                                                                      ifelse(is.null(self$proportion_per_seq), 1,
                                                                                             self$proportion_per_seq) * 100, "% per sequence and maximum of ",
                                                                                      self$max_samples, " samples per file.")
                                                              }
                                                              tensorflow::tf$summary$scalar(name = paste("training files seen (%)"), data = files_percentage, step = epoch,
                                                                                            description = description)
                                                            }
                                                            
                                                            file.writer$flush()
                                                          }
                                                          
                                                        ))
    
    
    l[[4]] <- proportion_training_files_cb(num_train_files = num_train_files, path_file_log = path_file_log, path_tensorboard = path_tensorboard, run_name = run_name,
                                           vocabulary_label = vocabulary_label, path = path, train_type = train_type, proportion_per_seq = proportion_per_seq,
                                           max_samples = max_samples, step = step, proportion_entries = proportion_entries)
    #names(l) <- c("hyper_param_model_outside", "tensorboard", "function_args","proportion_training_files")
  } else {
    #names(l) <- c("hyper_param_model_outside", "tensorboard", "function_args")
  }
  return(l)
}

#' Reset states callback
#' 
#' Reset states at start/end of validation and whenever file changes. Can be used for stateful LSTM.
#' 
#' @param path_file_log Path to log of training files.
#' @param path_file_logVal  Path to log of validation files.
#' @export
reset_states_cb <- function(path_file_log, path_file_logVal) {
  
  reset_states_cb_py_class <- reticulate::PyClass("reset_states_cb",
                                                  inherit = tensorflow::tf$keras$callbacks$Callback,
                                                  list(
                                                    
                                                    `__init__` = function(self, path_file_log, path_file_logVal) {
                                                      self$path_file_log <- path_file_log
                                                      self$path_file_logVal <- path_file_logVal
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
                                                      files_used <- readLines(self$path_file_log)
                                                      self$num_files_new <- length(files_used)
                                                      if (self$num_files_new > self$num_files_old) {
                                                        self$model$reset_states()
                                                        self$num_files_old <- self$num_files_new
                                                      }
                                                    },
                                                    
                                                    on_test_batch_begin = function(self, batch, logs) {
                                                      files_used <- readLines(self$path_file_logVal)
                                                      self$num_files_new_val <- length(files_used)
                                                      if (self$num_files_new_val > self$num_files_old_val) {
                                                        self$model$reset_states()
                                                        self$num_files_old_val <- self$num_files_new_val
                                                      }
                                                    }
                                                    
                                                  ))
  
  reset_states_cb_py_class(path_file_log = path_file_log, path_file_logVal = path_file_logVal)
}

#' Validation after training callback
#' 
#' Do validation only once at end of training.
#' 
#' @param gen.val Validation generator
#' @param validation_steps Number of validation steps.
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

#' Confusion matrix callback.
#' 
#' Create a confusion matrix to display under tensorboard images. 
#'
#' @inheritParams train_model
#' @param confMatLabels Names of classes.
#' @param cm_dir Directory that contains confusion matrix files.
#' @export
conf_matrix_cb <- function(path_tensorboard, run_name, confMatLabels, cm_dir, total_epochs) {
  
  conf_matrix_cb_py_class <- reticulate::PyClass("conf_matrix_cb",
                                                 inherit = tensorflow::tf$keras$callbacks$Callback,
                                                 list(
                                                   
                                                   `__init__` = function(self, cm_dir, path_tensorboard, run_name, confMatLabels, graphics = "png", total_epochs) {
                                                     self$cm_dir <- cm_dir
                                                     self$path_tensorboard <- path_tensorboard
                                                     self$run_name <- run_name
                                                     self$plot_path_train <- tempfile(pattern = "", fileext = paste0(".", graphics))
                                                     self$plot_path_val <- tempfile(pattern = "", fileext = paste0(".", graphics))
                                                     self$confMatLabels <- confMatLabels
                                                     self$epoch <- 0
                                                     self$total_epochs <- total_epochs
                                                     self$train_images <- NULL
                                                     self$val_images <- NULL
                                                     self$graphics <- graphics
                                                     self$text_size <- NULL
                                                     self$round_dig <- 3
                                                     if (length(confMatLabels) < 8) {
                                                       self$text_size <- (10 - (max(nchar(confMatLabels)) * 0.15)) * (0.95^length(confMatLabels))
                                                     }
                                                     self$cm_display_percentage <- TRUE
                                                     NULL
                                                   },
                                                   
                                                   on_epoch_begin = function(self, epoch, logs) {
                                                     suppressMessages(library(yardstick))
                                                     if (epoch > 0 & epoch < self$total_epochs) {
                                                       
                                                       cm_train <- readRDS(file.path(self$cm_dir, paste0("cm_train_", epoch-1, ".rds")))
                                                       cm_val <- readRDS(file.path(self$cm_dir, paste0("cm_val_", epoch, ".rds")))
                                                       if (self$cm_display_percentage) {
                                                         cm_train <- cm_perc(cm_train, self$round_dig)
                                                         cm_val <- cm_perc(cm_val, self$round_dig)
                                                       }
                                                       cm_train <- create_conf_mat_obj(cm_train, self$confMatLabels)
                                                       cm_val <- create_conf_mat_obj(cm_val, self$confMatLabels)
                                                       
                                                       
                                                       suppressMessages(
                                                         cm_plot_train <- ggplot2::autoplot(cm_train, type = "heatmap") +
                                                           ggplot2::scale_fill_gradient(low="#D6EAF8", high = "#2E86C1")  +
                                                           ggplot2::theme(axis.text.x =
                                                                            ggplot2::element_text(angle=90,hjust=1, size = self$text_size)) +
                                                           ggplot2::theme(axis.text.y =
                                                                            ggplot2::element_text(size = self$text_size))
                                                       )
                                                       
                                                       suppressMessages(
                                                         cm_plot_val <- ggplot2::autoplot(cm_val, type = "heatmap") +
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
                                                         
                                                         suppressMessages(ggplot2::ggsave(filename = self$plot_path_train, plot = cm_plot_train, device = "png",
                                                                                          width = plot_size,
                                                                                          height = plot_size,
                                                                                          units = "cm"))
                                                         p_cm_train <- png::readPNG(self$plot_path_train)
                                                         
                                                         suppressMessages(ggplot2::ggsave(filename = self$plot_path_val, plot = cm_plot_val, device = "png",
                                                                                          width = plot_size,
                                                                                          height = plot_size,
                                                                                          units = "cm"))
                                                         p_cm_val <- png::readPNG(self$plot_path_val)
                                                         
                                                       } else {
                                                         
                                                         suppressMessages(ggplot2::ggsave(filename = self$plot_path_train, plot = cm_plot_train, device = "jpg",
                                                                                          width = plot_size,
                                                                                          height = plot_size,
                                                                                          units = "cm"))
                                                         p_cm_train <- jpeg::readJPEG(self$plot_path_train)
                                                         
                                                         suppressMessages(ggplot2::ggsave(filename = self$plot_path_val, plot = cm_plot_val, device = "jpg",
                                                                                          width = plot_size,
                                                                                          height = plot_size,
                                                                                          units = "cm"))
                                                         p_cm_train <- jpeg::readJPEG(self$plot_path_val)
                                                       }
                                                       
                                                       p_cm_train <- as.array(p_cm_train)
                                                       p_cm_train <- array(p_cm_train, dim = c(1, dim(p_cm_train)))
                                                       p_cm_val <- as.array(p_cm_val)
                                                       p_cm_val <- array(p_cm_val, dim = c(1, dim(p_cm_val)))
                                                       
                                                       num_images <- 1
                                                       train_images <- array(0, dim = c(num_images, dim(p_cm_train)[-1]))
                                                       train_images[1, , , ] <- p_cm_train
                                                       self$train_images <- train_images
                                                       
                                                       val_images <- array(0, dim = c(num_images, dim(p_cm_val)[-1]))
                                                       val_images[1, , , ] <- p_cm_val
                                                       self$val_images <- val_images
                                                       file.writer <- tensorflow::tf$summary$create_file_writer(file.path(self$path_tensorboard, self$run_name))
                                                       file.writer$set_as_default()
                                                       tensorflow::tf$summary$image(name = "confusion matrix train", data = self$train_images, step = as.integer(epoch-1))
                                                       tensorflow::tf$summary$image(name = "confusion matrix validation", data = self$val_images, step = as.integer(epoch-1))
                                                       file.writer$flush()
                                                       self$epoch <- epoch
                                                     }
                                                   },
                                                   
                                                   on_train_end = function(self, logs) {
                                                     
                                                     epoch <- self$epoch + 1
                                                     cm_train <- readRDS(file.path(self$cm_dir, paste0("cm_train_", epoch-1, ".rds")))
                                                     #cm_val <- readRDS(file.path(self$cm_dir, paste0("cm_val_", epoch, ".rds")))
                                                     # extract last val confusion metric from custom metric
                                                     for (i in 1:length(model$metrics)) {
                                                       if (model$metrics[[i]]$name == "balanced_acc") {
                                                         bal_acc_index <- i
                                                         break
                                                       }
                                                     }
                                                     cm_val <- as.array(model$metrics[[bal_acc_index]]$cm)
                                                     
                                                     if (self$cm_display_percentage) {
                                                       cm_train <- cm_perc(cm_train, self$round_dig)
                                                       cm_val <- cm_perc(cm_val, self$round_dig)
                                                     }
                                                     cm_train <- create_conf_mat_obj(cm_train, self$confMatLabels)
                                                     cm_val <- create_conf_mat_obj(cm_val, self$confMatLabels)
                                                     
                                                     
                                                     suppressMessages(
                                                       cm_plot_train <- ggplot2::autoplot(cm_train, type = "heatmap") +
                                                         ggplot2::scale_fill_gradient(low="#D6EAF8", high = "#2E86C1")  +
                                                         ggplot2::theme(axis.text.x =
                                                                          ggplot2::element_text(angle=90,hjust=1, size = self$text_size)) +
                                                         ggplot2::theme(axis.text.y =
                                                                          ggplot2::element_text(size = self$text_size))
                                                     )
                                                     
                                                     suppressMessages(
                                                       cm_plot_val <- ggplot2::autoplot(cm_val, type = "heatmap") +
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
                                                       
                                                       suppressMessages(ggplot2::ggsave(filename = self$plot_path_train, plot = cm_plot_train, device = "png",
                                                                                        width = plot_size,
                                                                                        height = plot_size,
                                                                                        units = "cm"))
                                                       p_cm_train <- png::readPNG(self$plot_path_train)
                                                       
                                                       suppressMessages(ggplot2::ggsave(filename = self$plot_path_val, plot = cm_plot_val, device = "png",
                                                                                        width = plot_size,
                                                                                        height = plot_size,
                                                                                        units = "cm"))
                                                       p_cm_val <- png::readPNG(self$plot_path_val)
                                                       
                                                     } else {
                                                       
                                                       suppressMessages(ggplot2::ggsave(filename = self$plot_path_train, plot = cm_plot_train, device = "jpg",
                                                                                        width = plot_size,
                                                                                        height = plot_size,
                                                                                        units = "cm"))
                                                       p_cm_train <- jpeg::readJPEG(self$plot_path_train)
                                                       
                                                       suppressMessages(ggplot2::ggsave(filename = self$plot_path_val, plot = cm_plot_val, device = "jpg",
                                                                                        width = plot_size,
                                                                                        height = plot_size,
                                                                                        units = "cm"))
                                                       p_cm_train <- jpeg::readJPEG(self$plot_path_val)
                                                     }
                                                     
                                                     p_cm_train <- as.array(p_cm_train)
                                                     p_cm_train <- array(p_cm_train, dim = c(1, dim(p_cm_train)))
                                                     p_cm_val <- as.array(p_cm_val)
                                                     p_cm_val <- array(p_cm_val, dim = c(1, dim(p_cm_val)))
                                                     
                                                     num_images <- 1
                                                     train_images <- array(0, dim = c(num_images, dim(p_cm_train)[-1]))
                                                     train_images[1, , , ] <- p_cm_train
                                                     self$train_images <- train_images
                                                     
                                                     val_images <- array(0, dim = c(num_images, dim(p_cm_val)[-1]))
                                                     val_images[1, , , ] <- p_cm_val
                                                     self$val_images <- val_images
                                                     file.writer <- tensorflow::tf$summary$create_file_writer(file.path(self$path_tensorboard, self$run_name))
                                                     file.writer$set_as_default()
                                                     tensorflow::tf$summary$image(name = "confusion matrix train", data = self$train_images, step = as.integer(epoch-1))
                                                     tensorflow::tf$summary$image(name = "confusion matrix validation", data = self$val_images, step = as.integer(epoch-1))
                                                     file.writer$flush()
                                                   }
                                                 ))
  conf_matrix_cb_py_class(path_tensorboard = path_tensorboard,
                          run_name = run_name,
                          confMatLabels = confMatLabels,
                          cm_dir = cm_dir,
                          total_epochs = total_epochs)
}


get_callbacks <- function(default_arguments, path_tensorboard, run_name, train_type,
                          path_model, path, train_val_ratio, batch_size, epochs, format,
                          max_queue_size, lr_plateau_factor, patience, cooldown, path_checkpoint,
                          steps_per_epoch, step, shuffle_file_order, initial_epoch, vocabulary,
                          learning_rate, shuffle_input, vocabulary_label, solver,
                          file_limit, reverse_complement, wavenet_format,  cnn_format,
                          create_model_function = NULL, vocabulary_size, gen_cb, argumentList,
                          maxlen, labelGen, labelByFolder, vocabulary_label_size, tb_images,
                          target_middle, path_file_log, proportion_per_seq, n_gram,
                          train_val_split_csv, model = NULL,
                          skip_amb_nuc, max_samples, proportion_entries, path_log, output,
                          train_with_gen, random_sampling, reduce_lr_on_plateau,
                          save_weights_only, save_best_only, reset_states, early_stopping_time,
                          validation_only_after_training, gen.val, target_from_csv) {
  
  
  if (output$checkpoints) {
    # create folder for checkpoints using run_name
    # filenames contain epoch, validation loss and validation accuracy
    checkpoint_dir <- paste0(path_checkpoint, "/", run_name)
    dir.create(checkpoint_dir, showWarnings = FALSE)
    if (!is.list(model$output)) {
      filepath_checkpoints <- file.path(checkpoint_dir, "Ep.{epoch:03d}-val_loss{val_loss:.2f}-val_acc{val_acc:.3f}.hdf5")
    } else {
      filepath_checkpoints <- file.path(checkpoint_dir, "Ep.{epoch:03d}.hdf5")
      if (save_best_only) {
        warning("save_best_only not implemented for multi target. Setting save_best_only to FALSE")
        save_best_only <- FALSE
      }
    }
  }
  
  # Check if path_file_log is unique
  if (!is.null(path_file_log) && dir.exists(path_file_log)) {
    stop(paste0("path_file_log entry is already present. Please give this file a unique name."))
  }
  
  count_files <- !random_sampling
  callbacks <- list()
  callback_names <- NULL
  
  if (reduce_lr_on_plateau) {
    if (is.list(model$outputs)) {
      monitor <- "val_loss"
    } else {
      monitor <- "val_acc"
    }
    callbacks[[1]] <- reduce_lr_cb(patience = patience, cooldown = cooldown,
                                   lr_plateau_factor = lr_plateau_factor,
                                   monitor = monitor)
    callback_names <- c("reduce_lr", callback_names)
  }
  
  if (!is.null(path_log)) {
    callbacks <- c(callbacks, log_cb(path_log, run_name))
    callback_names <- c("log", callback_names)
  }
  
  if (!output$tensorboard) tb_images <- FALSE
  if (output$tensorboard) {
    
    # add balanced acc score
    
    # initialize metrics, temporary fix
    model <- manage_metrics(model)
    metrics <- model$metrics
    if (train_with_gen) {
      num_targets <- ifelse(train_type == "lm", length(vocabulary), length(vocabulary_label))
    } else {
      num_targets <- dim(dataset$Y)[2]
    }
    contains_macro_acc_metric <- FALSE
    for (i in 1:length(model$metrics)) {
      if (model$metrics[[i]]$name == "balanced_acc") contains_macro_acc_metric <- TRUE
    }
    
    if (!contains_macro_acc_metric) {
      if (tb_images) {
        if (!reticulate::py_has_attr(model, "cm_dir")) {
          cm_dir <- file.path(tempdir(), paste(sample(letters, 7), collapse = ""))
          dir.create(cm_dir)
          model$cm_dir <- cm_dir
        }
        metrics <- c(metrics, balanced_acc_wrapper(num_targets = as.integer(num_targets), cm_dir = model$cm_dir))
      }
    }
    
    # count files in path
    if (train_type == "label_rds" | train_type == "lm_rds") format <- "rds"
    if (train_with_gen) {
      num_train_files <- count_files(path = path, format = format, train_type = train_type, 
                                     target_from_csv = target_from_csv, 
                                     train_val_split_csv = train_val_split_csv)
    } else {
      num_train_files <- 1
    }
    
    complete_tb <- tensorboard_complete_cb(default_arguments = default_arguments, model = model, path_tensorboard = path_tensorboard, run_name = run_name, train_type = train_type,
                                           path_model = path_model, path = path, train_val_ratio = train_val_ratio, batch_size = batch_size, epochs = epochs,
                                           max_queue_size = max_queue_size, lr_plateau_factor = lr_plateau_factor, patience = patience, cooldown = cooldown,
                                           steps_per_epoch = steps_per_epoch, step = step, shuffle_file_order = shuffle_file_order, initial_epoch = initial_epoch, vocabulary = vocabulary,
                                           learning_rate = learning_rate, shuffle_input = shuffle_input, vocabulary_label = vocabulary_label, solver = solver,
                                           file_limit = file_limit, reverse_complement = reverse_complement, wavenet_format = wavenet_format,  cnn_format = cnn_format,
                                           create_model_function = NULL, vocabulary_size = vocabulary_size, gen_cb = gen_cb, argumentList = argumentList,
                                           maxlen = maxlen, labelGen = labelGen, labelByFolder = labelByFolder, vocabulary_label_size = vocabulary_label_size, tb_images = FALSE,
                                           target_middle = target_middle, num_train_files = num_train_files, path_file_log = path_file_log, proportion_per_seq = proportion_per_seq,
                                           skip_amb_nuc = skip_amb_nuc, max_samples = max_samples, proportion_entries = proportion_entries,
                                           train_with_gen = train_with_gen, count_files = !random_sampling)
    callbacks <- c(callbacks, complete_tb)
    callback_names <- c(callback_names, names(complete_tb))
  }
  
  if (output$checkpoints) {
    if (wavenet_format) {
      # can only save weights for wavenet
      save_weights_only <- TRUE
    }
    callbacks <- c(callbacks, checkpoint_cb(filepath = filepath_checkpoints, save_weights_only = save_weights_only,
                                            save_best_only = save_best_only))
    callback_names <- c(callback_names, "checkpoint")
  }
  
  if (reset_states) {
    callbacks <- c(callbacks, reset_states_cb(path_file_log = path_file_log, path_file_logVal = path_file_logVal))
    callback_names <- c(callback_names, "reset_states")
  }
  
  if (!is.null(early_stopping_time)) {
    callbacks <- c(callbacks, early_stopping_cb(early_stopping_patience = early_stopping_patience,
                                                early_stopping_time = early_stopping_time))
    callback_names <- c(callback_names, "early_stopping")
  }
  
  if (validation_only_after_training) {
    if (!train_with_gen) stop("Validation after training only implemented for generator")
    callbacks <- c(callbacks, validation_after_training_cb(gen.val = gen.val, validation_steps = validation_steps))
    callback_names <- c(callback_names, "validation_after_training")
  }
  
  if (tb_images) {
    if (is.list(model$output)) {
      warning("Tensorboard images (confusion matrix) not implemented for model with multiple outputs.
                 Setting tb_images to FALSE")
      tb_images <- FALSE
    }
    
    if (model$loss == "binary_crossentropy") {
      warning("Tensorboard images (confusion matrix) not implemented for sigmoid activation in last layer.
                 Setting tb_images to FALSE")
      tb_images <- FALSE
    }
  }
  
  if (tb_images) {
    
    confMatLabels <- vocabulary_label
    if (train_with_gen & train_type == "lm") {
      if (is.null(n_gram) || n_gram == 1) {
        confMatLabels <- vocabulary
      } else {
        l <- list()
        for (i in 1:n_gram) {
          l[[i]] <- vocabulary
        }
        confMatLabels <- expand.grid(l) %>% apply(1, paste0) %>% apply(2, paste, collapse = "") %>% sort()
      }
    }
    
    model %>% keras::compile(loss = model$loss,
                             optimizer = model$optimizer, metrics = metrics[-1])
    model <- manage_metrics(model)
    
    if (length(confMatLabels) > 16) {
      message("Cannot display confusion matrix with more than 16 labels.")
    } else {
      
      callbacks <- c(callbacks, conf_matrix_cb(path_tensorboard = path_tensorboard,
                                               run_name = run_name,
                                               confMatLabels = confMatLabels,
                                               cm_dir = model$cm_dir,
                                               total_epochs = epochs))
      callback_names <- c(callback_names, "conf_matrix")
    }
  }
  return(callbacks)
}
