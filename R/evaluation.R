#' Evaluates a trained model on fasta, fastq or rds files
#'
#' Returns evaluation metric like confusion matrix, loss, AUC, AUPRC, MAE, MSE (depending on output layer).
#'
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_folder
#' @inheritParams generator_fasta_label_header_csv
#' @param path_input Input directory where fasta, fastq or rds files are located.
#' @param model A keras model.
#' @param batch_size Number of samples per batch.
#' @param step How often to take a sample.
#' @param vocabulary Vector of allowed characters. Character outside vocabulary get encoded as specified in ambiguous_nuc.
#' @param vocabulary_label List of labels for targets of each output layer.
#' @param number_batches How many batches to evaluate.
#' @param format File format, `"fasta"`, `"fastq"` or `"rds"`.
#' @param mode Either `"lm"` for language model or `"label_header"`, `"label_csv"` or `"label_folder"` for label classification.
#' @param verbose Boolean.
#' @param target_middle Whether model is language model with separate input layers. 
#' @param evaluate_all_files Boolean, if `TRUE` will iterate over all files in \code{path_input} once. \code{number_batches} will be overwritten.
#' @param auc Whether to include AUC metric. If output layer activation is `"softmax"`, only possible for 2 targets. Computes the average if output layer has sigmoid
#' activation and multiple targets.
#' @param auprc Whether to include AUPRC metric. If output layer activation is `"softmax"`, only possible for 2 targets. Computes the average if output layer has sigmoid
#' activation and multiple targets.
#' @param path_pred_list Path to store list of predictions (output of output layers) and corresponding true labels as rds file. 
#' @param exact_num_samples Exact number of samples to evaluate. If you want to evaluate a number of samples not divisible by batch_size. Useful if you want
#' to evaluate a data set exactly ones and know the number of samples already. Should be a vector if `mode = "label_folder"` (with same length as `vocabulary_label`)
#' and else an integer.
#' @param activations List containing output formats for output layers (`softmax, sigmoid` or `linear`). If `ǸULL`, will be estimated from model.   
#' @param ... Further generator options. See \code{\link{get_generator}}.
#' @examples
#' # create dummy data
#' path_input <- tempfile()
#' dir.create(path_input)
#' create_dummy_data(file_path = path_input,
#'                   num_files = 3,
#'                   seq_length = 11, 
#'                   num_seq = 5,
#'                   vocabulary = c("a", "c", "g", "t"))
#' # create model
#' model <- create_model_lstm_cnn(layer_lstm = 8, layer_dense = 4, maxlen = 10, verbose = FALSE)
#' # evaluate
#' evaluate_model(path_input = path_input,
#'   model = model,
#'   step = 11,
#'   vocabulary = c("a", "c", "g", "t"),
#'   vocabulary_label = list(c("a", "c", "g", "t")),
#'   mode = "lm",
#'   output_format = "target_right",
#'   evaluate_all_files = TRUE,
#'   verbose = FALSE)
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
                           concat_seq = NULL,
                           seed = 1234,
                           auc = FALSE,
                           auprc = FALSE,
                           path_pred_list = NULL,
                           exact_num_samples = NULL,
                           activations = NULL,
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
  if (is.null(activations)) activations <- get_output_activations(model)
  
  if (is.null(vocabulary_label)) vocabulary_label <- list(vocabulary)
  if (!is.list(vocabulary_label)) vocabulary_label <- list(vocabulary_label)
  if (mode == "label_folder") {
    number_batches <- rep(ceiling(number_batches/length(path_input)), length(path_input))
  }
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
  maxlen <- 1000
  
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
        if (length(files) == 0) {
          stop("No files from path_input have label in target_from_csv file.")
        }
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
        
        if (!is.null(concat_seq)) {
          seq_vector <- paste(seq_vector, collapse = concat_seq)
        }
        
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
      x <- rds_file[[1]]
      while (is.list(x)) {
        x <- x[[1]]
      }
      num_samples <- dim(x)[1] + num_samples
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
                              concat_seq = concat_seq,
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
                                            concat_seq = concat_seq,
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
      gen <- generator_fasta_label_folder(path_corpus = path_input[[k]],
                                          format = format,
                                          batch_size = batch_size,
                                          maxlen = maxlen,
                                          max_iter = max_iter,
                                          vocabulary = vocabulary,
                                          step = step,
                                          padding = padding,
                                          concat_seq = concat_seq,
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
          
          if (is.list(y_conf)) {
            for (m in 1:length(y_conf)) {
              y_conf[[m]] <- y_conf[[m]][index, ]
              y[[m]] <- y[[m]][index, ]
            }
          } else {
            y_conf <- y_conf[index, ]
            y <- y[index, ]
          }
          
          # vector to matrix
          if (length(index) == 1) {
            if (is.list(y_conf)) {
              for (m in 1:length(y_conf)) {
                y_conf[[m]] <- array(as.array(y_conf[[m]]), dim = c(1, length(y_conf[[m]])))
                y[[m]] <- matrix(y[[m]], ncol = length(y[[m]]))
              }
            } else {
              y_conf <- array(as.array(y_conf), dim = c(1, length(y_conf)))
              y <- matrix(y, ncol = length(y))
            }
          }
          
        }
      }
      
      y_conf_list[[count]] <- y_conf
      if (batch_size == 1 | (!is.null(index) && length(index == 1))) {
        col_num <- ncol(y_conf)
        if (is.na(col_num)) col_num <- length(y_conf)
        y_list[[count]] <- matrix(y, ncol = col_num)
      } else {
        y_list[[count]] <- y
      }
      count <- count + 1
      
      if (verbose & (batch_index == 10)) {
        time_passed <- as.double(difftime(Sys.time(), start_time, units = "hours"))
        time_estimation <- (overall_num_batches/10) * time_passed
        cat("Evaluation will take approximately", round(time_estimation, 3), "hours. Starting time:", format(Sys.time(), "%F %R."), " \n")
        
      }
      
      if (verbose & (batch_index > ten_percent_steps[percentage_index]) & percentage_index < 10) {
        cat("Progress: ", percentage_index * 10 ,"% \n")
        time_passed <- as.double(difftime(Sys.time(), start_time, units = "hours"))
        cat("Time passed: ", round(time_passed, 3), "hours \n")
        percentage_index <- percentage_index + 1
      }
      
    }
  }
  
  if (verbose) {
    cat("Progress: 100 % \n")
    time_passed <- as.double(difftime(Sys.time(), start_time, units = "hours"))
    cat("Time passed: ", round(time_passed, 3), "hours \n")
  }
  
  y_conf_list <- reshape_y_list(y_conf_list, num_out_layers = num_out_layers, tf_format = TRUE)
  y_list <- reshape_y_list(y_list, num_out_layers = num_out_layers, tf_format = FALSE)
  
  if (!is.null(path_pred_list)) {
    saveRDS(list(pred = y_conf_list, true = y_list), path_pred_list)
  }
  
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
      eval_list[[i]] <- evaluate_linear(y_true = y_list[[i]], y_pred = y_conf_list[[i]], label_names = vocabulary_label[[i]])
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

#' Evaluate matrices of true targets and predictions from layer with softmax activation. 
#' 
#' Compute confusion matrix, accuracy, categorical crossentropy and (optionally) AUC or AUPRC, given predictions and
#' true targets. AUC and AUPRC only possible for 2 targets. 
#' 
#' @param y Matrix of true target.
#' @param y_conf Matrix of predictions.
#' @param auc Whether to include AUC metric. Only possible for 2 targets. 
#' @param auprc Whether to include AUPRC metric. Only possible for 2 targets. 
#' @param label_names Names of corresponding labels. Length must be equal to number of columns of \code{y}.
#' @examples
#' y <- matrix(c(1, 0, 0, 0, 1, 1), ncol = 2)
#' y_conf <- matrix(c(0.3, 0.5, 0.1, 0.7, 0.5, 0.9), ncol = 2)
#' evaluate_softmax(y, y_conf, auc = TRUE, auprc = TRUE, label_names = c("A", "B")) 
#' @export    
evaluate_softmax <- function(y, y_conf, auc = FALSE, auprc = FALSE, label_names = NULL) {
  
  if (ncol(y) != 2 & (auc | auprc)) {
    message("Can only compute AUC or AUPRC if output layer with softmax acticvation has two neurons.")
    auc <- FALSE
    auprc <- FALSE
  }
  
  y_pred <- apply(y_conf, 1, which.max)
  y_true <- apply(y, 1, FUN = which.max) - 1
  
  df_true_pred <- data.frame(
    true = factor(y_true + 1, levels = 1:(length(label_names)), labels = label_names),
    pred = factor(y_pred, levels = 1:(length(label_names)), labels = label_names)
  )
  
  loss_per_class <- list()
  for (i in 1:ncol(y)) {
    index <- (y_true + 1) == i
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
      weights.class0 = y_true)
  } else {
    auc_list <- NULL
  }
  
  if (auprc) {
    auprc_list <- PRROC::pr.curve(
      scores.class0 = y_conf[ , 2],
      weights.class0 = y_true)
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

#' Evaluate matrices of true targets and predictions from layer with sigmoid activation. 
#' 
#' Compute accuracy, binary crossentropy and (optionally) AUC or AUPRC, given predictions and
#' true targets. Outputs columnwise average.  
#' 
#' @inheritParams evaluate_model
#' @inheritParams evaluate_softmax
#' @param auc Whether to include AUC metric.
#' @param auprc Whether to include AUPRC metric. 
#' @examples
#' y <- matrix(sample(c(0, 1), 30, replace = TRUE), ncol = 3)
#' y_conf <- matrix(runif(n = 30), ncol = 3)
#' evaluate_sigmoid(y, y_conf, auc = TRUE, auprc = TRUE)
#' 
#' @export    
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
    
    na_count <- sum(is.na(auc_vector))
    if (na_count > 0) {
      message(paste(sum(na_count), ifelse(na_count > 1, "columns", "column"),
                    "removed from AUC evaluation since they contain only one label"))
    }
    AUC <- mean(auc_vector, na.rm = TRUE)
  } else {
    AUC <- NULL
  }
  
  if (auprc) {
    auprc_list <- purrr::map(1:ncol(y_conf), ~PRROC::pr.curve(
      scores.class0 = y_conf[ , .x],
      weights.class0 = y[ , .x]))
    auprc_vector <- vector("numeric", ncol(y))
    for (i in 1:length(auprc_vector)) {
      auprc_vector[i] <- auprc_list[[i]]$auc.integral
    }
    AUPRC <- mean(auprc_vector, na.rm = TRUE) 
  } else {
    AUPRC <- NULL
  }
  
  return(list(accuracy = acc,
              binary_crossentropy_loss = loss,
              #loss_per_class = loss_per_class,
              #accuracy_per_class = class_acc,
              AUC = AUC,
              AUPRC = AUPRC))
  
}

#' Evaluate matrices of true targets and predictions from layer with linear activation. 
#' 
#' Compute MAE and MSE, given predictions and
#' true targets. Outputs columnwise average.  
#' 
#' @inheritParams evaluate_model
#' @inheritParams evaluate_softmax
#' @param y_true Matrix of true labels.
#' @param y_pred Matrix of predictions.
#' @examples 
#' y_true <- matrix(rnorm(n = 12), ncol = 3)
#' y_pred <- matrix(rnorm(n = 12), ncol = 3)
#' evaluate_linear(y_true, y_pred)
#' @export    
evaluate_linear <- function(y_true, y_pred, label_names = NULL) {
  
  loss_per_class_mse <- list()
  loss_per_class_mae <- list()
  for (i in 1:ncol(y_true)) {
    mse_loss_class <- tensorflow::tf$keras$losses$mean_squared_error(y_true[ ,i], y_pred[ , i])
    mae_loss_class <- tensorflow::tf$keras$losses$mean_absolute_error(y_true[ ,i], y_pred[ , i])
    loss_per_class_mse[[i]] <- mse_loss_class$numpy()
    loss_per_class_mae[[i]] <- mae_loss_class$numpy()
  }
  
  return(list(mse = mean(unlist(loss_per_class_mse)),
              mae = mean(unlist(loss_per_class_mae))))
  
}


#' Plot ROC
#' 
#' Compute ROC and AUC from target and prediction matrix and plot ROC. Target/prediction matrix should 
#' have one column if output of layer with sigmoid activation and two columns for softmax activation. 
#' 
#' @inheritParams evaluate_softmax
#' @inheritParams evaluate_linear
#' @param path_roc_plot Where to store ROC plot.
#' @param return_plot Whether to return plot.
#' @examples
#' y_true <- matrix(c(1, 0, 0, 0, 1, 1), ncol = 1)
#' y_conf <- matrix(runif(n = nrow(y_true)), ncol = 1)
#' p <- plot_roc(y_true, y_conf, return_plot = TRUE)
#' p
#' @export    
plot_roc <- function(y_true, y_conf, path_roc_plot = NULL,
                     return_plot = TRUE) {
  
  if (!all(y_true == 0 | y_true == 1)) {
    stop("y_true should only contain 0 and 1 entries")
  }
  
  if (is.matrix(y_true) && ncol(y_true) > 2) {
    stop("y_true can contain 1 or 2 columns")
  }
  
  if (is.matrix(y_true) && ncol(y_true) == 2) {
    y_true <- y_true[ , 1] 
    y_conf <- y_conf[ , 2]
  }
  
  if (var(y_true) == 0) {
    stop("y_true contains just one label")
  }
  
  y_true <- as.vector(y_true)
  y_conf <- as.vector(y_conf)
  
  rocobj <-  pROC::roc(y_true, y_conf, quiet = TRUE)
  auc <- round(pROC::auc(y_true, y_conf, quiet = TRUE), 4)
  p <- pROC::ggroc(rocobj,  size = 1, color = "black")
  p <- p + ggplot2::theme_classic() + ggplot2::theme(aspect.ratio = 1) 
  p <- p + ggplot2::ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')'))
  p <- p + ggplot2::geom_abline(intercept = 1, linetype = 2, color = "grey50")
  p <- p + ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "grey50")
  p <- p + ggplot2::geom_hline(yintercept = 1,  linetype = 2, color = "grey50")
  
  if (!is.null(path_roc_plot)) {
    ggplot2::ggsave(path_roc_plot, p)
  }
  
  if (return_plot) {
    return(p)
  } else {
    return(NULL)
  }
  
}

# plot_roc_auprc <- function(y_true, y_conf, path_roc_plot = NULL, path_auprc_plot = NULL,
#                            return_plot = TRUE, layer_activation = "softmax") {
#   
#   if (layer_activation == "softmax") {
#     
#     if (!all(y_true == 0 | y_true == 1)) {
#       stop("y_true should only contain 0 and 1 entries")
#     }
#     
#     if (ncol(y_true) != 2 & (auc | auprc)) {
#       message("Can only compute AUC or AUPRC if output layer with softmax acticvation has two neurons.")
#     }
#     
#     auc_list <- PRROC::roc.curve(
#       scores.class0 = y_conf[ , 2],
#       weights.class0 = y_true[ , 2], curve = TRUE)
#     
#     
#     auprc_list <- PRROC::pr.curve(
#       scores.class0 = y_conf[ , 2],
#       weights.class0 = y_true[ , 2], curve = TRUE)
#     
#     #auc_plot <- NULL
#     #auprc_plot <- NULL  
#     
#   }
#   
#   if (layer_activation == "sigmoid") {
#     
#     auc_list <- purrr::map(1:ncol(y_conf), ~PRROC::roc.curve(
#       scores.class0 = y_conf[ , .x],
#       weights.class0 = y[ , .x], curve = TRUE))
#     auc_vector <- vector("numeric", ncol(y))
#     
#     
#     auprc_list <- purrr::map(1:ncol(y_conf), ~PRROC::pr.curve(
#       scores.class0 = y_conf[ , .x],
#       weights.class0 = y[ , .x], curve = TRUE))
#     auprc_vector <- vector("numeric", ncol(y))
#     
#   }
#   
#   if (!is.null(path_roc_plot)) {
#     
#   }
#   
#   if (!is.null(path_auprc_plot)) {
#     
#   }
#   
# }
