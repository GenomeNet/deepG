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
                      "tensorboard.log", "fileLog", "train_type", "ambiguous_nuc", "added_label_path", "added_label_names",
                      "train_val_split_csv", "target_from_csv")
  } else {
    argsInQuotes <- c("model_path", "path", "path.val", "checkpoint_path", "run.name", "solver", "output_format",
                      "tensorboard.log", "fileLog", "train_type", "ambiguous_nuc", "format", "added_label_path", "added_label_names",
                      "train_val_split_csv", "target_from_csv")
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
                                    proportion_per_file, skip_amb_nuc, max_samples, proportion_entries) {
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

                                                        `__init__` = function(self, num_train_files, fileLog, tensorboard.log, run.name, labelVocabulary,
                                                                              path, train_type, start_index, proportion_per_file, max_samples, step,
                                                                              proportion_entries) {
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
                                                          self$proportion_entries <- proportion_entries
                                                          NULL
                                                        },

                                                        on_epoch_end = function(self, epoch, logs) {
                                                          if (is.null(self$proportion_entries)) self$proportion_entries <- 1
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
                                                                  self$description[i] <- paste0("Using step size ", self$step[i], ", proportion_entries ",
                                                                                                self$proportion_entries * 100, "% and ",
                                                                                                ifelse(is.null(self$proportion_per_file[i]), 1,
                                                                                                       self$proportion_per_file[i]) * 100, "% per sequence")
                                                                } else {
                                                                  self$description[i] <- paste0("Using step size ", self$step[i], ", ",
                                                                                                ifelse(is.null(self$proportion_per_file[i]), 1,
                                                                                                       self$proportion_per_file[i]) * 100, "% per sequence, maximum of ",
                                                                                                self$max_samples, " samples per file and proportion_entries ",
                                                                                                self$proportion_entries * 100, "%")
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
                                                              description <- paste0("Using step size ", step,
                                                                                    ", proportion_entries ", self$proportion_entries * 100, "% and ",
                                                                                    ifelse(is.null(self$proportion_per_file), 1,
                                                                                           self$proportion_per_file) * 100, "% per sequence")
                                                            } else {
                                                              description <- paste0("Using step size ", step, ", ",
                                                                                    ifelse(is.null(self$proportion_per_file), 1,
                                                                                           self$proportion_per_file) * 100, "% per sequence, maximum of ",
                                                                                    self$max_samples, " samples per file and proportion_entries ",
                                                                                    self$proportion_entries * 100, "%")

                                                            }
                                                            if (self$train_type == "label_rds") {
                                                              description <- paste0("Using step size ",
                                                                                    ifelse(is.null(self$proportion_per_file), 1,
                                                                                           self$proportion_per_file) * 100, "% per sequence and maximum of ",
                                                                                    self$max_samples, " samples per file.")
                                                            }
                                                            tensorflow::tf$summary$scalar(name = paste("training files seen (%)"), data = files_percentage, step = epoch,
                                                                                          description = description)
                                                          }

                                                          file.writer$flush()
                                                        }

                                                      ))


  l[[4]] <- proportion_training_files_cb(num_train_files = num_train_files, fileLog = fileLog, tensorboard.log = tensorboard.log, run.name = run.name,
                                         labelVocabulary = labelVocabulary, path = path, train_type = train_type, proportion_per_file = proportion_per_file,
                                         max_samples = max_samples, step = step, proportion_entries = proportion_entries)

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
conf_matrix_cb <- function(tensorboard.log, run.name, confMatLabels, cm_dir) {

  conf_matrix_cb_py_class <- reticulate::PyClass("conf_matrix_cb",
                                                 inherit = tensorflow::tf$keras$callbacks$Callback,
                                                 list(

                                                   `__init__` = function(self, cm_dir, tensorboard.log, run.name, confMatLabels, graphics = "png") {
                                                     self$cm_dir <- cm_dir
                                                     self$tensorboard.log <- tensorboard.log
                                                     self$run.name <- run.name
                                                     self$plot_path_train <- tempfile(pattern = "", fileext = paste0(".", graphics))
                                                     self$plot_path_val <- tempfile(pattern = "", fileext = paste0(".", graphics))
                                                     self$confMatLabels <- confMatLabels
                                                     self$epoch <- 0
                                                     self$train_images <- NULL
                                                     self$val_images <- NULL
                                                     self$graphics <- graphics
                                                     self$epoch <- 0
                                                     self$text_size <- NULL
                                                     if (length(confMatLabels) < 8) {
                                                       self$text_size <- (10 - (max(nchar(confMatLabels)) * 0.15)) * (0.95^length(confMatLabels))
                                                     }
                                                     self$cm_display_percentage <- TRUE
                                                     NULL
                                                   },

                                                   on_epoch_begin = function(self, epoch, logs) {
                                                     suppressMessages(library(yardstick))
                                                     if (epoch > 0) {

                                                       cm_train <- readRDS(file.path(self$cm_dir, paste0("cm_train_", epoch-1, ".rds")))
                                                       cm_val <- readRDS(file.path(self$cm_dir, paste0("cm_val_", epoch, ".rds")))
                                                       if (self$cm_display_percentage) {
                                                         cm_train <- cm_perc(cm_train, 2)
                                                         cm_val <- cm_perc(cm_val, 2)
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
                                                       file.writer <- tensorflow::tf$summary$create_file_writer(file.path(self$tensorboard.log, self$run.name))
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
                                                     cm_val <- readRDS(file.path(self$cm_dir, paste0("cm_val_", epoch, ".rds")))
                                                     if (self$cm_display_percentage) {
                                                       cm_train <- cm_perc(cm_train, 2)
                                                       cm_val <- cm_perc(cm_val, 2)
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
                                                     file.writer <- tensorflow::tf$summary$create_file_writer(file.path(self$tensorboard.log, self$run.name))
                                                     file.writer$set_as_default()
                                                     tensorflow::tf$summary$image(name = "confusion matrix train", data = self$train_images, step = as.integer(epoch-1))
                                                     tensorflow::tf$summary$image(name = "confusion matrix validation", data = self$val_images, step = as.integer(epoch-1))
                                                     file.writer$flush()
                                                   }
                                                 ))
  conf_matrix_cb_py_class(tensorboard.log = tensorboard.log,
                          run.name = run.name,
                          confMatLabels = confMatLabels,
                          cm_dir = cm_dir)
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

#' Balanced  accuracy metric.
#'
#' @export
balanced_acc_wrapper <- function(num_targets, cm_dir) {
  balanced_acc_stateful <- reticulate::PyClass("balanced_acc",
                                               inherit = tensorflow::tf$keras$metrics$Metric,
                                               list(

                                                 `__init__` = function(self, num_targets, cm_dir) {
                                                   super()$`__init__`(name = "balanced_acc")
                                                   self$num_targets <- num_targets
                                                   self$cm_dir <- cm_dir
                                                   self$count <- 0
                                                   self$cm <- self$add_weight(name = "cm_matrix", shape = c(num_targets, num_targets), initializer="zeros")
                                                   NULL
                                                 },

                                                 update_state = function(self, y_true, y_pred, sample_weight = NULL) {
                                                   self$cm$assign_add(self$compute_cm(y_true, y_pred))
                                                   NULL
                                                 },

                                                 result = function(self) {
                                                   balanced_acc <- self$compute_balanced_acc()
                                                   self$store_cm()
                                                   return(balanced_acc)
                                                 },

                                                 compute_cm = function(self, y_true, y_pred) {
                                                   labels <- tensorflow::tf$math$argmax(y_true, axis = 1L)
                                                   predictions <- tensorflow::tf$math$argmax(y_pred, axis = 1L)
                                                   current_cm <- tensorflow::tf$math$confusion_matrix(
                                                     labels = labels, predictions = predictions,
                                                     dtype = "float32", num_classes = self$num_targets)
                                                   current_cm <- tensorflow::tf$transpose(current_cm)
                                                   return(current_cm)
                                                 },

                                                 compute_balanced_acc = function(self) {
                                                   diag <- tensorflow::tf$linalg$diag_part(self$cm)
                                                   col_sums <- tensorflow::tf$math$reduce_sum(self$cm, axis=0L)
                                                   average_per_class <- tensorflow::tf$math$divide(diag, col_sums)
                                                   nan_index <- tensorflow::tf$math$logical_not(tensorflow::tf$math$is_nan(average_per_class))
                                                   average_per_class <- tensorflow::tf$boolean_mask(average_per_class, nan_index)
                                                   acc_sum <- tensorflow::tf$math$reduce_sum(average_per_class)
                                                   balanced_acc <- tensorflow::tf$math$divide(acc_sum, length(average_per_class))
                                                   return(balanced_acc)
                                                 },

                                                 reset_states = function(self) {
                                                   self$count <- self$count + 1
                                                   self$cm$assign_sub(self$cm)
                                                   NULL
                                                 },

                                                 store_cm = function(self) {
                                                   if (self$count %% 2 == 0) {
                                                     file_name <- file.path(self$cm_dir, paste0("cm_val_", floor(self$count/2), ".rds"))
                                                   } else {
                                                     file_name <- file.path(self$cm_dir, paste0("cm_train_", floor(self$count/2), ".rds"))
                                                   }
                                                   saveRDS(keras::k_eval(self$cm), file_name)
                                                   NULL
                                                 }

                                               ))
  return(balanced_acc_stateful(num_targets = num_targets, cm_dir = cm_dir))
}

#' F1 metric
#'
#' @export
f1_wrapper <- function(num_targets = 2) {
  f1_stateful <- reticulate::PyClass("f1",
                                     inherit = tensorflow::tf$keras$metrics$Metric,
                                     list(

                                       `__init__` = function(self, num_targets) {
                                         super()$`__init__`(name = "f1")
                                         self$num_targets <- num_targets
                                         self$f1_score <- 0
                                         self$cm <- self$add_weight(name = "cm_matrix", shape = c(num_targets, num_targets), initializer="zeros")
                                         NULL
                                       },

                                       update_state = function(self, y_true, y_pred, sample_weight = NULL) {
                                         self$cm$assign_add(self$compute_cm(y_true, y_pred))
                                         NULL
                                       },

                                       result = function(self) {
                                         self$f1_score <- self$compute_f1()
                                         return(self$f1_score)
                                       },

                                       compute_cm = function(self, y_true, y_pred) {
                                         labels <- tensorflow::tf$math$argmax(y_true, axis = 1L)
                                         predictions <- tensorflow::tf$math$argmax(y_pred, axis = 1L)
                                         current_cm <- tensorflow::tf$math$confusion_matrix(
                                           labels = labels, predictions = predictions,
                                           dtype = "float32", num_classes = self$num_targets)
                                         current_cm <- tensorflow::tf$transpose(current_cm)
                                         return(current_cm)
                                       },

                                       compute_f1 = function(self) {
                                         diag <- tensorflow::tf$linalg$diag_part(self$cm)
                                         precision <- diag/(tensorflow::tf$reduce_sum(self$cm, 0L) + tensorflow::tf$constant(1e-15))
                                         recall <- diag/(tensorflow::tf$reduce_sum(self$cm, 1L) + tensorflow::tf$constant(1e-15))
                                         f1 = (2 * precision * recall)/(precision + recall + tensorflow::tf$constant(1e-15))
                                         return(f1)
                                       },

                                       reset_states = function(self) {
                                         self$cm$assign_sub(self$cm)
                                         NULL
                                       }

                                     ))
  return(f1_stateful(num_targets = num_targets))
}

#' Mean AUC score (mean or median)
#'
#' @param model_output_size Number of neurons in output layer of model, for which metric will be applied to.
#' @param loss Loss function of model, for which metric will be applied to; must be "binary_crossentropy"
#' or "catergorical_crossentropy".
#' @export
auc_wrapper <- function(model_output_size,
                        loss = "binary_crossentropy") {

  stopifnot(loss %in% c("binary_crossentropy", "categorical_crossentropy"))
  if (loss == "categorical_crossentropy" & model_output_size != 2) {
    stop("Output size must be two, when loss is catergorical_crossentropy")
  }
  metric_name <- ifelse(loss == "binary_crossentropy" & model_output_size > 1,
                        "mean_AUC", "AUC")

  auc_stateful <- reticulate::PyClass("AUC",
                                      inherit = tensorflow::tf$keras$metrics$Metric,
                                      list(

                                        `__init__` = function(self, model_output_size, loss, metric_name) {
                                          super()$`__init__`(name = metric_name)
                                          self$model_output_size <- model_output_size
                                          self$loss <- loss
                                          if (loss == "binary_crossentropy") {
                                            self$auc_scores <- self$add_weight(name = "auc_vector", shape = c(1, model_output_size), initializer="zeros")
                                            self$auc_score <-  self$add_weight(name = "auc_score", initializer="zeros")
                                            self$auc_list <- vector("list", model_output_size)
                                            for (i in 1:model_output_size) {
                                              assign(paste0("m_", i), tensorflow::tf$keras$metrics$AUC())
                                            }
                                            #purrr::map(1:model_output_size, ~assign(paste0("m_", .x), tensorflow::tf$keras$metrics$AUC()))
                                            parse_text <- purrr::map(1:model_output_size, ~parse(text = paste0("m_", .x)))
                                            self$auc_list <- purrr::map(1:model_output_size, ~eval(parse_text[[.x]]))
                                            self$shape <- c(1L, as.integer(self$model_output_size))
                                          } else {
                                            self$auc_scores <- self$add_weight(name = "auc_vector", shape = c(1, 1), initializer="zeros")
                                            self$auc_score <-  self$add_weight(name = "auc_score", initializer="zeros")
                                            self$auc_list <- vector("list", 1)
                                            assign(paste0("m_", 1), tensorflow::tf$keras$metrics$AUC())
                                            parse_text <- purrr::map(1, ~parse(text = paste0("m_", .x)))
                                            self$auc_list <- purrr::map(1, ~eval(parse_text[[.x]]))
                                            self$shape <- c(1L, 1L)
                                          }
                                          NULL
                                        },

                                        update_state = function(self, y_true, y_pred, sample_weight = NULL) {
                                          self$compute_auc(y_true, y_pred)
                                          current_auc_list <- vector("list", length(self$auc_list))
                                          # for (i in 1:length(self$auc_list)) {
                                          #   current_auc_list[[i]] <- self$auc_list[[i-1]]$result()
                                          # }
                                          current_auc_list <-  purrr::map(1:length(self$auc_list),
                                                                          ~self$auc_list[[.x - 1]]$result())
                                          current_auc <- unlist(current_auc_list)
                                          current_auc <- tensorflow::tf$reshape(tensor = current_auc, shape = self$shape)
                                          self$auc_scores$assign(current_auc)
                                          NULL
                                        },

                                        result = function(self) {
                                          self$auc_score$assign(tensorflow::tf$math$reduce_mean(self$auc_scores))
                                          return(self$auc_score)
                                        },

                                        compute_auc = function(self, y_true, y_pred) {

                                          if (self$loss == "binary_crossentropy") {
                                            if (self$model_output_size > 1) {
                                              # for (i in 0:(length(self$auc_list) - 1)) {
                                              #   self$auc_list[[i]]$update_state(y_true[ , i+1], y_pred[ , i+1])
                                              # }
                                              purrr::map(0:(length(self$auc_list) - 1),
                                                         ~self$auc_list[[.x]]$update_state(y_true[ , .x+1], y_pred[ , .x+1]))
                                            } else {
                                              self$auc_list[[0]]$update_state(y_true, y_pred)
                                            }

                                          } else {
                                            y_true_temp <- y_true[ , 1]
                                            y_pred_temp <- tensorflow::tf$math$argmax(y_pred, axis = 1L)
                                            self$auc_list[[0]]$update_state(y_true_temp, y_pred_temp)
                                          }
                                          NULL
                                        },

                                        reset_states = function(self) {
                                          # for (i in 0:(length(self$auc_list) - 1)) {
                                          #   self$auc_list[[i]]$reset_states()
                                          # }
                                          purrr::map(0:(length(self$auc_list) - 1),
                                                     ~self$auc_list[[.x]]$reset_states())
                                          self$auc_scores$assign_sub(self$auc_scores)
                                          NULL
                                        }

                                      ))

  return(auc_stateful(model_output_size = model_output_size, loss = loss, metric_name = metric_name))
}
