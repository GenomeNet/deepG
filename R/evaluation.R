#' Evaluates a trained model on fasta/fastq or rds files
#'
#' Returns evaluation metric like confusion matrix, loss, AUC, AUPRC, MAE, MSE (depending and output layer).
#' Evaluates \code{batch_size} * \code{number_batches} samples.
#'
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_folder
#' @inheritParams generator_fasta_label_header_csv
#' @param path_input Input directory where fasta/fastq or rds files are located.
#' @param model A keras model.
#' @param batch_size Number of samples per batch.
#' @param step How often to take a sample.
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as specified in ambiguous_nuc.
#' @param vocabulary_label List of labels for targets of each output layer.
#' @param number_batches How many batches to evaluate.
#' @param format File format, "fasta", "fastq" or "rds".
#' @param mode Either "lm" for language model and "label_header", "label_csv" or "label_folder" for label classification.
#' @param evaluate_all_files Boolean, if TRUE will iterate over all files in \code{path_input} once. \code{number_batches} will be overwritten.
#' @param auc Whether to include AUC metric. Only possible for 2 targets if layer activation is "softmax".
#' @param auprc Whether to include AUPRC metric. Only possible for 2 targets if layer activation is "softmax".
#' @param exact_num_samples Exact number of samples to evaluate. If you want to evaluate a number of samples not devisible by batch_size. Useful if you want
#' to evaluate a data set exactly ones and know the number of samples already. Should be a vector if mode = "label_folder" (with same length as vocabulary_label)
#' and else an integer.
#' @param ... Further generator options.
#' @export
evaluate_model <- function(path_input,
                           model = NULL,
                           batch_size = 100,
                           step = 1,
                           padding = FALSE,
                           vocabulary = c("a", "c", "g", "t"),
                           vocabulary_label = list(c("a", "c", "g", "t")),
                           number_batches = 10,
                           format = "fasta",
                           target_middle = FALSE,
                           mode = "lm",
                           output_format = "target_right",
                           ambiguous_nuc = "zero",
                           evaluate_all_files = FALSE,
                           verbose = TRUE,
                           max_iter = 20000,
                           target_from_csv = NULL,
                           max_samples = NULL,
                           proportion_per_seq = NULL,
                           seed = 1234,
                           auc = FALSE,
                           auprc = FALSE,
                           exact_num_samples = NULL,
                           ...) {

  set.seed(seed)
  path_model <- NULL
  stopifnot(mode %in% c("lm", "label_header", "label_folder", "label_csv", "lm_rds", "label_rds"))
  stopifnot(format %in% c("fasta", "fastq", "rds"))
  stopifnot(is.null(proportion_per_seq) || proportion_per_seq <= 1)
  if (!is.null(exact_num_samples) & evaluate_all_files) {
    warning(paste("Will evaluate number of samples as specified in exact_num_samples argument. Setting evaluate_all_files to FALSE."))
    evaluate_all_files <- FALSE
  }
  eval_exact_num_samples <- !is.null(exact_num_samples) | evaluate_all_files
  activations <- get_output_activations(model)

  if (is.null(vocabulary_label)) vocabulary_label <- list(vocabulary)
  number_batches <- rep(ceiling(number_batches/length(path_input)), length(path_input))
  num_classes <- ifelse(mode == "label_folder", length(path_input), 1)
  num_out_layers <- length(activations)

  # extract maxlen from model
  num_in_layers <- length(model$inputs)
  if (num_in_layers == 1) {
    maxlen <- model$input$shape[[2]]
  } else {
    if (!target_middle) {
      maxlen <- model$input[[num_in_layers]]$shape[[2]]
    } else {
      maxlen <- model$input[[num_in_layers - 1]]$shape[[2]] + model$input[[num_in_layers]]$shape[[2]]
    }
  }

  if (evaluate_all_files & (format %in% c("fasta", "fastq"))) {

    number_batches <- NULL
    num_samples <- rep(0, length(path_input))

    for (i in 1:num_classes) {
      if (mode == "label_folder") {
        files <- list_fasta_files(path_input[[i]], format = format, file_filter = NULL)
      } else {
        files <- list_fasta_files(path_input, format = format, file_filter = NULL)
      }

      # remove files not in csv table
      if (mode == "label_csv") {
        csv_file <- read.csv2(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
        if (dim(csv_file)[2] == 1) {
          csv_file <- read.csv(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
        }
        index <- basename(files) %in% csv_file$file
        files <- files[index]
      }

      for (file in files) {
        if (format == "fasta") {
          fasta_file <- microseq::readFasta(file)
        } else {
          fasta_file <- microseq::readFastq(file)
        }

        # remove entries with wrong header
        if (mode == "label_header") {
          index <- fasta_file$Header %in% vocabulary_label
          fasta_file <- fasta_file[index, ]
        }

        seq_vector <- fasta_file$Sequence

        if (!is.null(proportion_per_seq)) {
          fasta_width <- nchar(seq_vector)
          sample_range <- floor(fasta_width - (proportion_per_seq * fasta_width))
          start <- mapply(sample_range, FUN = sample, size = 1)
          perc_length <- floor(fasta_width * proportion_per_seq)
          stop <- start + perc_length
          seq_vector <- mapply(seq_vector, FUN = substr, start = start, stop = stop)
        }

        if (mode == "lm") {
          if (!padding) {
            seq_vector <- seq_vector[nchar(seq_vector) >= (maxlen + 1)]
          } else {
            length_vector <- nchar(seq_vector)
            short_seq_index <- which(length_vector < (maxlen + 1))
            for (ssi in short_seq_index) {
              seq_vector[ssi] <- paste0(paste(rep("0", (maxlen + 1) - length_vector[ssi]), collapse = ""), seq_vector[ssi])
            }
          }
        } else {
          if (!padding) {
            seq_vector <- seq_vector[nchar(seq_vector) >= (maxlen)]
          } else {
            length_vector <- nchar(seq_vector)
            short_seq_index <- which(length_vector < (maxlen))
            for (ssi in short_seq_index) {
              seq_vector[ssi] <- paste0(paste(rep("0", (maxlen) - length_vector[ssi]), collapse = ""), seq_vector[ssi])
            }
          }
        }

        if (length(seq_vector) == 0) next
        new_samples <- get_start_ind(seq_vector = seq_vector,
                                     length_vector = nchar(seq_vector),
                                     maxlen = maxlen,
                                     step = step,
                                     train_mode = ifelse(mode == "lm", "lm", "label"),
                                     discard_amb_nuc = ifelse(ambiguous_nuc == "discard", TRUE, FALSE),
                                     vocabulary = vocabulary
        ) %>% length()

        if (is.null(max_samples)) {
          num_samples[i] <- num_samples[i] + new_samples
        } else {
          num_samples[i] <- num_samples[i] + min(new_samples, max_samples)
        }
      }
      number_batches[i] <- ceiling(num_samples[i]/batch_size)

    }
    if (mode == "label_folder") {
      message_string <- paste0("Evaluate ", num_samples, " samples for class ", vocabulary_label[[1]], ".\n")
    } else {
      message_string <- paste0("Evaluate ", sum(num_samples), " samples.")
    }
    message(message_string)
  }

  if (evaluate_all_files & format == "rds") {
    rds_files <- list_fasta_files(path_corpus = file_path,
                                  format = "rds",
                                  file_filter = NULL)
    num_samples <- 0
    for (file in rds_files) {
      rds_file <- readRDS(file)
      num_samples <- dim(rds_file[[1]])[1] + num_samples
    }
    number_batches <- ceiling(num_samples/batch_size)
    message_string <- paste0("Evaluate ", num_samples, " samples.")
    message(message_string)
  }

  if (!is.null(exact_num_samples)) {
    num_samples <- exact_num_samples
    number_batches <- ceiling(num_samples/batch_size)
  }

  overall_num_batches <- sum(number_batches)

  if (mode == "lm") {
    gen <- generator_fasta_lm(path_corpus = path_input,
                              format = format,
                              batch_size = batch_size,
                              maxlen = maxlen,
                              max_iter = max_iter,
                              vocabulary = vocabulary,
                              verbose = FALSE,
                              shuffle_file_order = FALSE,
                              step = step,
                              padding = padding,
                              shuffle_input = FALSE,
                              reverse_complement = FALSE,
                              output_format = output_format,
                              ambiguous_nuc = ambiguous_nuc,
                              proportion_per_seq = proportion_per_seq,
                              max_samples = max_samples,
                              seed = seed,
                              ...)
  }

  if (mode == "label_header" | mode == "label_csv") {
    gen <- generator_fasta_label_header_csv(path_corpus = path_input,
                                            format = format,
                                            batch_size = batch_size,
                                            maxlen = maxlen,
                                            max_iter = max_iter,
                                            vocabulary = vocabulary,
                                            verbose = FALSE,
                                            shuffle_file_order = FALSE,
                                            step = step,
                                            padding = padding,
                                            shuffle_input = FALSE,
                                            vocabulary_label = vocabulary_label[[1]],
                                            reverse_complement = FALSE,
                                            ambiguous_nuc = ambiguous_nuc,
                                            target_from_csv = target_from_csv,
                                            proportion_per_seq = proportion_per_seq,
                                            max_samples = max_samples,
                                            seed = seed, ...)
  }

  if (mode == "label_rds" | mode == "lm_rds") {
    gen <- generator_rds(rds_folder = path_input, batch_size = batch_size, path_file_log = NULL, ...)
  }

  batch_index <- 1
  start_time <- Sys.time()
  ten_percent_steps <- seq(overall_num_batches/10, overall_num_batches, length.out = 10)
  percentage_index <- 1
  count <- 1
  y_conf_list <- vector("list", overall_num_batches)
  y_list <- vector("list", overall_num_batches)

  for (k in 1:num_classes) {

    index <- NULL
    if (mode == "label_folder") {
      gen <- generator_fasta_label_folder(path_corpus = path_input[k],
                                          format = format,
                                          batch_size = batch_size,
                                          maxlen = maxlen,
                                          max_iter = max_iter,
                                          vocabulary = vocabulary,
                                          step = step,
                                          padding = padding,
                                          reverse_complement = FALSE,
                                          num_targets = length(path_input),
                                          ones_column = k,
                                          ambiguous_nuc = ambiguous_nuc,
                                          proportion_per_seq = proportion_per_seq,
                                          max_samples = max_samples,
                                          seed = seed, ...)
    }

    for (i in 1:number_batches[k]) {
      z <- gen()
      x <- z[[1]]
      y <- z[[2]]

      y_conf <- model(x)
      batch_index <- batch_index + 1

      # remove double predictions
      if (eval_exact_num_samples & (i == number_batches[k])) {
        double_index <- (i * batch_size) - num_samples[k]
        if (double_index > 0) {
          index <- 1:(nrow(y_conf) - double_index)
          y_conf <- y_conf[index, ]
          y <- y[index, ]
        }
      }

      y_conf_list[[count]] <- y_conf
      if (batch_size == 1 | (!is.null(index) && length(index == 1))) {
        y_list[[count]] <- matrix(y, ncol = ncol(y_conf))
      } else {
        y_list[[count]] <- y
      }
      count <- count + 1

      if (verbose & (batch_index == 10)) {
        time_passed <- as.double(difftime(Sys.time(), start_time, units = "hours"))
        time_estimation <- (overall_num_batches/10) * time_passed
        cat("Evaluation will take approximately", round(time_estimation, 3), "hours. Starting time:", format(Sys.time(), "%F %R."), " \n")

      }

      if (verbose & (batch_index > ten_percent_steps[percentage_index])) {
        cat("Progress: ", percentage_index * 10 ,"% \n")
        time_passed <- as.double(difftime(Sys.time(), start_time, units = "hours"))
        cat("Time passed: ", round(time_passed, 3), "hours \n")
        percentage_index <- percentage_index + 1
      }

    }
  }

  y_conf_list <- reshape_y_list(y_conf_list, num_out_layers = num_out_layers, tf_format = TRUE)
  y_list <- reshape_y_list(y_list, num_out_layers = num_out_layers, tf_format = FALSE)

  eval_list <- list()
  for (i in 1:num_out_layers) {

    if (activations[i] == "softmax") {
      eval_list[[i]] <- evaluate_softmax(y = y_list[[i]], y_conf = y_conf_list[[i]],
                                         auc = auc, auprc = auprc,
                                         label_names = vocabulary_label[[i]])
    }

    if (activations[i] == "sigmoid") {
      eval_list[[i]] <- evaluate_sigmoid(y = y_list[[i]], y_conf = y_conf_list[[i]],
                                         auc = auc, auprc = auprc,
                                         label_names = vocabulary_label[[i]])
    }

    if (activations[i] == "linear") {
      eval_list[[i]] <- evaluate_linear(y = y_list[[i]], y_conf = y_conf_list[[i]], label_names = vocabulary_label[[i]])
    }

  }

  return(eval_list)
}



reshape_y_list <- function(y, num_out_layers, tf_format = TRUE) {

  if (num_out_layers > 1) {
    y <- do.call(c, y)
  }

  reshaped_list <- vector("list", num_out_layers)

  for (i in 1:num_out_layers) {
    index <- seq(i, length(y), by = num_out_layers)
    if (tf_format) {
      reshaped_list[[i]] <- y[index] %>%
        tensorflow::tf$concat(axis = 0L) %>%
        keras::k_eval()
    } else {
      reshaped_list[[i]] <- do.call(rbind, y[index])
    }
  }
  return(reshaped_list)
}


get_output_layer_names <- function(model) {
  out_layers <- model$get_config()$output
  names_vec <- vector("character", length(out_layers))
  for (i in 1:length(out_layers)) {
    names_vec[i] <- out_layers[[i]][[1]]
  }
  return(names_vec)
}

get_output_activations <- function(model) {
  out_names <- get_output_layer_names(model)
  num_layers <- length(model$get_config()$layers)

  act_vec <- vector("character", length(out_names))
  count <- 1
  for (i in 1:num_layers) {
    layer_name <- model$get_config()$layers[[i]]$name
    if (layer_name %in% out_names) {
      act_name <- model$layers[[i]]$get_config()$activation
      if (is.null(act_name)) act_name <- "linear"
      act_vec[count] <- act_name
      count <- count + 1
    }
  }
  return(act_vec)
}

evaluate_softmax <- function(y, y_conf, auc = FALSE, auprc = FALSE, label_names = NULL) {

  if (ncol(y) != 2 & (auc | auprc)) {
    message("Can only compute AUC or AUPRC if output layer with softmax acticvation has two neurons.")
    auc <- FALSE
    auprc <- FALSE
  }

  y_pred <- apply(y_conf, 1, which.max)
  y_true <- apply(y, 1, FUN = which.max)

  df_true_pred <- data.frame(
    true = factor(y_true, levels = 1:(length(label_names)), labels = label_names),
    pred = factor(y_pred, levels = 1:(length(label_names)), labels = label_names)
  )

  loss_per_class <- list()
  for (i in 1:ncol(y)) {
    index <- y_true == i
    if (any(index)) {
      cce_loss_class <- tensorflow::tf$keras$losses$categorical_crossentropy(y[index, ], y_conf[index, ])
      loss_per_class[[i]] <- cce_loss_class$numpy()
    } else {
      loss_per_class[[i]] <- NA
    }
  }

  cm <- yardstick::conf_mat(df_true_pred, true, pred)
  confMat <- cm[[1]]

  acc <- sum(diag(confMat))/sum(confMat)
  loss <- mean(unlist(loss_per_class))

  for (i in 1:length(loss_per_class)) {
    loss_per_class[[i]] <- mean(unlist(loss_per_class[[i]]), na.rm = TRUE)
  }

  loss_per_class <- unlist(loss_per_class)
  m <- as.matrix(confMat)
  class_acc <- vector("numeric")
  for (i in 1:ncol(m)) {
    if (sum(m[ , i]) == 0) {
      class_acc[i] <- NA
    } else {
      class_acc[i] <- m[i, i]/sum(m[ , i])
    }
  }
  names(class_acc) <- label_names
  names(loss_per_class) <- label_names
  balanced_acc <- mean(class_acc)

  if (auc) {
    auc_list <- PRROC::roc.curve(
      scores.class0 = y_conf[ , 2],
      weights.class0 = y_true - 1)
  } else {
    auc_list <- NULL
  }

  if (auprc) {
    auprc_list <- PRROC::pr.curve(
      scores.class0 = y_conf[ , 2],
      weights.class0 = y_true - 1)
  } else {
    auprc_list <- NULL
  }

  return(list(confusion_matrix = confMat,
              accuracy = acc,
              categorical_crossentropy_loss = loss,
              #balanced_accuracy = balanced_acc,
              #loss_per_class = loss_per_class,
              #accuracy_per_class = class_acc,
              AUC = auc_list$auc,
              AUPRC = auprc_list$auc.integral))
}


evaluate_sigmoid <- function(y, y_conf, auc = FALSE, auprc = FALSE, label_names = NULL) {

  y_pred <- ifelse(y_conf > 0.5, 1, 0)

  loss_per_class <- list()
  for (i in 1:ncol(y)) {
    bce_loss_class <- tensorflow::tf$keras$losses$binary_crossentropy(y[ , i], y_conf[ , i])
    loss_per_class[[i]] <- bce_loss_class$numpy()
  }

  loss_per_class <- unlist(loss_per_class)
  names(loss_per_class) <- label_names
  loss <- mean(unlist(loss_per_class))

  class_acc <- vector("numeric", ncol(y))
  for (i in 1:ncol(y)) {
    num_true_pred <-  sum(y[ , i] == y_pred[ , i])
    class_acc[i] <- num_true_pred /nrow(y)
  }
  names(class_acc) <- label_names
  acc <- mean(class_acc)

  if (auc) {
    auc_list <- purrr::map(1:ncol(y_conf), ~PRROC::roc.curve(
      scores.class0 = y_conf[ , .x],
      weights.class0 = y[ , .x]))
    auc_vector <- vector("numeric", ncol(y))
    for (i in 1:length(auc_vector)) {
      auc_vector[i] <- auc_list[[i]]$auc
    }
  } else {
    auc_list <- NULL
  }

  if (auprc) {
    auprc_list <- purrr::map(1:ncol(y_conf), ~PRROC::pr.curve(
      scores.class0 = y_conf[ , .x],
      weights.class0 = y[ , .x]))
    auprc_vector <- vector("numeric", ncol(y))
    for (i in 1:length(auprc_vector)) {
      auprc_vector[i] <- auprc_list[[i]]$auc.integral
    }
  } else {
    auprc_list <- NULL
  }

  return(list(accuracy = acc,
              binary_crossentropy_loss = loss,
              #loss_per_class = loss_per_class,
              #accuracy_per_class = class_acc,
              AUC = mean(auc_vector),
              AUPRC = mean(auprc_vector)))

}

evaluate_linear <- function(y_true, y_conf, label_names = NULL) {

  loss_per_class_mse <- list()
  loss_per_class_mae <- list()
  for (i in 1:ncol(y)) {
    mse_loss_class <- tensorflow::tf$keras$losses$mean_squared_error(y_true[ ,i], y_conf[ , i])
    mae_loss_class <- tensorflow::tf$keras$losses$mean_absolute_error(y_true[ ,i], y_conf[ , i])
    loss_per_class_mse[[i]] <- mse_loss_class$numpy()
    loss_per_class_mae[[i]] <- mae_loss_class$numpy()
  }

  return(list(mse = mean(loss_per_class_mse),
              mae = mean(loss_per_class_mae)))

}
