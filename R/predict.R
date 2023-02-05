#' Predict next nucleotides in sequence 
#' 
#' The output is a S4 class.
#'
#' @param sequence input sequence, length should be in sync with the model.
#' If length exceeds input.shape of model then only the right side of the
#' sequence will be used.
#' @param model trained model from the function \code{trainNetwork()}
#' @param vocabulary vocabulary of input sequence
#' @param verbose TRUE/FALSE
#' @examples 
#' \dontrun{
#' example.model <- keras::load_model_hdf5("example_model.hdf5")
#' sequence <- strrep("A", 100)
#' predictNextNucleotide(sequence, example.model)}
#' @export
predictNextNucleotide <- function(sequence,
                                  model,
                                  vocabulary =  c("l", "a", "c", "g", "t"),
                                  verbose = F){
  
  stopifnot(!missing(sequence))
  stopifnot(!missing(model))
  stopifnot(nchar(sequence) >= model$input_shape[2])
  
  substringright <- function(x, n){
    substr(x, nchar(x)- n + 1, nchar(x))
  }
  # sequence can be longer then model input shape
  # if so just use the last input_shape chars
  
  sentence <- tokenizers::tokenize_characters(
    stringr::str_to_lower(substringright(sequence, as.numeric(model$input_shape[2]))),
    strip_non_alphanum = FALSE, simplify = TRUE)
  
  x <- sapply(vocabulary, function(x){
    as.numeric(x == sentence)
  })
  x <- keras::array_reshape(x, c(1, dim(x)))
  
  if(verbose) {
    message("Prediction ...")}
  
  preds <- keras::predict_proba(model, x)
  next_index <- which.max(preds)
  next_char <- vocabulary[next_index]
  # return a S4 class
  return(new("prediction",
             next_char = next_char,
             probability = preds[next_index],
             index = next_index,
             alternative_probability = preds,
             solution = paste0(sequence, next_char)))
}


#' Replaces specific nucleotides in a sequence
#' 
#' @param sequence input sequence, length should be in sync with the model.
#' If length exceeds input.shape of model then only the right side of the
#' sequence will be used.
#' @param model trained model from the function \code{trainNetwork()}
#' @param char character in the sequence that will be replaced
#' @param vocabulary ordered vocabulary of input sequence
#' @examples 
#' \dontrun{
#' example.model <- keras::load_model_hdf5("example_model.hdf5")
#' replaceChar(sequence = sequence, model = example.model)}
#' @export
replaceChar <- function(sequence,
                        model,
                        char = "X",
                        vocabulary =  c("l", "a", "c", "g", "t")){
  
  stopifnot(!missing(sequence))
  stopifnot(!missing(model))
  stopifnot(nchar(sequence) >= model$input_shape[2])
  
  while (stringr::str_detect(sequence, char)) {
    # get the position
    next_position <- stringr::str_locate_all(pattern = 'X', sequence)[[1]][[1]]
    # seed text for model is the most-right chunk of text
    # with size of model$input_shape[[2]]
    seed <- substr(sequence,
                   next_position - model$input_shape[[2]] - 1,
                   next_position - 1)
    prediction <- predictNextNucleotide(seed, model, vocabulary)
    sequence <- paste0(prediction@solution,
                       substr(sequence, next_position + 1,
                              nchar(sequence)))
  }
  return(sequence)
}

#' Evaluates a trained model on .fasta/fastq files
#' 
#' Returns accuracies per batch and overall confusion matrix. Evaluates \code{batch.size} * \code{numberOfBatches} samples.
#' 
#' @inheritParams fastaFileGenerator
#' @inheritParams labelByFolderGenerator
#' @inheritParams fastaLabelGenerator
#' @param fasta.path Input directory where fasta/fastq files are located.
#' @param model A keras model. 
#' @param batch.size Number of samples per batch.
#' @param step How often to take a sample.
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as 0-vector.
#' @param label_vocabulary Labels for targets. Equal to vocabulary if not given.
#' @param numberOfBatches How many batches to evaluate.
#' @param filePath Where to store output, if missing output won't be written.
#' @param format File format, "fasta" or "fastq".
#' @param filename Name of output file.
#' @param mode Either "lm" for language model and "label_header", "label_csv" or "label_folder" for label classification.
#' @inheritParams fastaFileGenerator
#' @param ambiguous_nuc How to handle nucleotides outside vocabulary, either "zero", "discard" or "equal". If "zero", input gets encoded as zero vector; 
#' if "equal" input is 1/length(vocabulary) x length(vocabulary). If "discard" samples containing nucleotides outside vocabulary get discarded.        
#' @param evaluate_all_files Boolean, if TRUE will iterate over all files in \code{fasta.path} once. \code{numberOfBatches} will be overwritten. 
#' @param auc Whether to include auc metric. Only possible for 2 targets. 
#' @export
evaluateFasta <- function(fasta.path,
                          model = NULL,
                          batch.size = 100,
                          step = 1,
                          padding = FALSE,
                          vocabulary = c("a", "c", "g", "t"),
                          label_vocabulary = c("a", "c", "g", "t"),
                          numberOfBatches = 10,
                          filePath = NULL,
                          format = "fasta",
                          filename = "",
                          target_middle = FALSE,
                          mode = "lm",
                          output_format = "target_right",
                          ambiguous_nuc = "zero",
                          evaluate_all_files = FALSE,
                          verbose = TRUE,
                          max_iter = 20000,
                          target_from_csv = NULL,
                          max_samples = NULL,
                          proportion_per_file = NULL,
                          seed = 1234,
                          auc = FALSE,
                          ...) {
  
  set.seed(seed)
  model.path <- NULL 
  stopifnot(mode %in% c("lm", "label_header", "label_folder", "label_csv"))
  stopifnot(format %in% c("fasta", "fastq"))
  stopifnot(is.null(proportion_per_file) || proportion_per_file <= 1)
  
  if (is.null(label_vocabulary)) label_vocabulary <- vocabulary
  numberOfBatches <- rep(ceiling(numberOfBatches/length(fasta.path)), length(fasta.path))
  num_classes <- ifelse(mode == "label_folder", length(fasta.path), 1)
  
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
  
  if (evaluate_all_files) {
    numberOfBatches <- NULL
    num_samples <- rep(0, length(fasta.path))
    
    for (i in 1:num_classes) {
      if (mode == "label_folder") {
        files <- list_fasta_files(fasta.path[[i]], format = format, file_filter = NULL)
      } else {
        files <- list_fasta_files(fasta.path, format = format, file_filter = NULL)
      }
      for (file in files) {
        if (format == "fasta") {
          seq_vector <- microseq::readFasta(file)$Sequence
        } else {
          seq_vector <- microseq::readFastq(file)$Sequence
        }
        
        if (!is.null(proportion_per_file)) {
          fasta_width <- nchar(seq_vector)
          sample_range <- floor(fasta_width - (proportion_per_file * fasta_width))
          start <- mapply(sample_range, FUN = sample, size = 1)
          perc_length <- floor(fasta_width * proportion_per_file)
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
        new_samples <- getStartInd(seq_vector = seq_vector,
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
      numberOfBatches[i] <- ceiling(num_samples[i]/batch.size)
      
    }
    if (mode == "label_folder") {
      message_string <- paste0("Evaluate ", num_samples, " samples for class ", label_vocabulary, 
                               ". Setting numberOfBatches to ", numberOfBatches, ".")
    } else {
      message_string <- paste0("Evaluate ", sum(num_samples), " samples. Setting numberOfBatches to ", sum(numberOfBatches), ".")
    }
    message(message_string)  
  }
  
  overall_num_batches <- sum(numberOfBatches)
  
  if (mode == "lm") {
    gen <- fastaFileGenerator(corpus.dir = fasta.path,
                              format = format,
                              batch.size = batch.size,
                              maxlen = maxlen,
                              max_iter = max_iter,
                              vocabulary = vocabulary,
                              verbose = FALSE,
                              randomFiles = FALSE,
                              step = step,
                              padding = padding,
                              showWarnings = FALSE,
                              shuffleFastaEntries = FALSE,
                              reverseComplements = FALSE,
                              output_format = output_format,
                              ambiguous_nuc = ambiguous_nuc,
                              proportion_per_file = proportion_per_file,
                              max_samples = max_samples,
                              seed = seed,
                              ...)
  }
  
  if (mode == "label_header" | mode == "label_csv") {
    gen <- fastaLabelGenerator(corpus.dir = fasta.path,
                               format = format,
                               batch.size = batch.size,
                               maxlen = maxlen,
                               max_iter = max_iter,
                               vocabulary = vocabulary,
                               verbose = FALSE,
                               randomFiles = FALSE,
                               step = step,
                               padding = padding,
                               showWarnings = FALSE,
                               shuffleFastaEntries = FALSE,
                               labelVocabulary = label_vocabulary,
                               reverseComplements = FALSE,
                               ambiguous_nuc = ambiguous_nuc,
                               target_from_csv = target_from_csv,
                               proportion_per_file = proportion_per_file,
                               max_samples = max_samples,
                               seed = seed,
                               ...)
  }
  
  
  acc <- vector("numeric")
  loss_per_class <- vector("list", length = length(fasta.path))
  confMat <- matrix(0, nrow = length(label_vocabulary), ncol = length(label_vocabulary))
  batch_index <- 1
  start_time <- Sys.time()
  ten_percent_steps <- seq(overall_num_batches/10, overall_num_batches, length.out = 10) 
  percentage_index <- 1
  if (auc) {
    if (length(label_vocabulary) != 2) {
      auc <- FALSE
      warning("AUC score only possible for 2 targets")
    } else {
      auc_score <- tensorflow::tf$keras$metrics$AUC()
    }  
  } 
  
  for (k in 1:num_classes) {
    
    cce_loss <- vector("list", length = numberOfBatches[k])
    
    if (mode == "label_folder") {
      # Bug, order of classes = order fasta.path ?  
      gen <- labelByFolderGenerator(corpus.dir = fasta.path[k],
                                    format = format,
                                    batch.size = batch.size,
                                    maxlen = maxlen,
                                    max_iter = max_iter,
                                    vocabulary = vocabulary,
                                    step = step, 
                                    padding = padding,
                                    reverseComplements = FALSE,
                                    numTargets = length(fasta.path),
                                    onesColumn = k,
                                    ambiguous_nuc = ambiguous_nuc,
                                    proportion_per_file = proportion_per_file,
                                    max_samples = max_samples,
                                    seed = seed,
                                    ...)
    }
    
    for (i in 1:numberOfBatches[k]) {
      z <- gen()
      x <- z[[1]]
      y <- z[[2]]
      
      y_conf <- model(x)
      y_pred <- apply(y_conf, 1, which.max)
      y_true <- apply(y, 1, FUN = which.max)
      batch_index <- batch_index + 1
      
      # remove double predictions
      if (evaluate_all_files & (i == numberOfBatches[k])) {
        double_index <- (i * batch.size) - num_samples[k]
        if (double_index > 0) {
          index <- 1:(length(y_true) - double_index)
          y_true <- y_true[index]
          y_pred <- y_pred[index]
          y_conf <- y_conf[index, ]
          y <- y[index, ]
        }
      }
      
      if (auc) auc_score$update_state(y_true-1, unlist(y_conf[ , 2]))
      
      df_true_pred <- data.frame(
        true = factor(y_true, levels = 1:(length(label_vocabulary)), labels = label_vocabulary),
        pred = factor(y_pred, levels = 1:(length(label_vocabulary)), labels = label_vocabulary)
      )
      
      cce_loss_new <- tensorflow::tf$keras$losses$categorical_crossentropy(y, y_conf)
      cce_loss[[i]] <- cce_loss_new$numpy()
      
      cm <- yardstick::conf_mat(df_true_pred, true, pred)
      confMat <- confMat + cm[[1]]
      
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
    loss_per_class[[k]] <- unlist(cce_loss)
  }
  
  acc <- sum(diag(confMat))/sum(confMat)
  loss <- mean(unlist(loss_per_class))
  
  for (i in 1:length(loss_per_class)) {
    loss_per_class[[i]] <- mean(unlist(loss_per_class[[i]]))
  }
  loss_per_class <- unlist(loss_per_class)
  m <- as.matrix(confMat)
  class_acc <- vector("numeric")
  for (i in 1:ncol(m)) {
    class_acc[i] <- m[i, i]/sum(m[ , i])
  }
  names(class_acc) <- label_vocabulary
  macro_acc <- mean(class_acc)
  
  if (mode == "label_folder") {
    
    names(loss_per_class) <- label_vocabulary
    
    output_list <- list(
      accuracy = acc,
      confusion_matrix = confMat,
      loss = loss,
      macro_acc = macro_acc,
      class_acc = class_acc,
      loss_per_class = loss_per_class) 
    
    if (auc) output_list <- c(output_list, AUC = auc_score$result()$numpy())
    
  } else {
    
    output_list <- list(
      accuracy = acc,
      confusion_matrix = confMat,
      loss = loss)
    if (auc) output_list <- c(output_list, AUC = auc_score$result()$numpy())
    
  }
  
  if (!is.null(filePath)) {
    saveRDS(output_list, file = paste0(filePath, "/", filename, "_outputList.rds"))
  }
  
  return(output_list)
}

#' Computes mean and median auc scores for model with 1 or multiple binary classifiers in last layer.
#'
#' @inheritParams evaluateFasta
#' @param file_path Path to folder containing rds files.
#' @param return_auc_vector Whether to return AUC vector containing scores for each class.  
#' @param return_auprc_vector Whether to return AUPRC vector containing scores for each class.  
#' @param evaluate_all_files Whether to evaluate all samples in \code{file_path} excactly once. 
#' @export 
auc_roc_metric <- function(file_path,
                           model = NULL,
                           batch_size = 100,
                           numberOfBatches = 10,
                           #roc_plot_path = NULL,
                           return_auc_vector = FALSE,
                           evaluate_all_files = FALSE, 
                           return_auprc_vector = FALSE, 
                           seed = 1234) {
  
  format <- "rds"
  plot_roc <- FALSE #!is.null(roc_plot_path)
  
  if (evaluate_all_files) {
    rds_files <- list_fasta_files(corpus.dir = file_path,
                                  format = "rds",
                                  file_filter = NULL)
    num_samples <- 0
    for (file in rds_files) {
      rds_file <- readRDS(file)
      num_samples <- dim(rds_file[[1]])[1] + num_samples
    }
    numberOfBatches <- ceiling(num_samples/batch_size)
    message_string <- paste0("Evaluate ", num_samples, " samples. Setting numberOfBatches to ", numberOfBatches, ".")
    message(message_string)
  }
  
  if (length(model$outputs) > 1) {
    stop("model should have only one output layer")
  }
  
  num_layers <- length(model$get_config()$layers)
  layer_name <- model$get_config()$layers[[num_layers]]$name
  activation_string <- as.character(model$get_layer(layer_name)$activation)
  if (!(stringr::str_detect(activation_string, "sigmoid"))) {
    stop("model should have sigmoid activation at last layer")
  }
  
  set.seed(seed)
  gen <- gen_rds(rds_folder = file_path, batch_size = batch_size, fileLog = NULL)
  y_conf_list <- vector("list", numberOfBatches)
  y_true_list <- vector("list", numberOfBatches)
  for (batch in 1:numberOfBatches) {
    z <- gen()
    x <- z[[1]]
    y <- z[[2]]
    y_conf <- predict(model, x, verbose = 0)
    
    if (batch_size == 1) {
      y_conf_list[[batch]] <- matrix(y_conf, nrow = 1) %>% as.data.frame() 
      y_true_list[[batch]] <- matrix(y, nrow = 1) %>% as.data.frame() 
    } else {
      y_conf_list[[batch]] <- y_conf %>% as.data.frame() 
      y_true_list[[batch]] <- y  %>% as.data.frame() 
    }
  }
  
  df_conf <- data.table::rbindlist(y_conf_list) %>% as.data.frame()
  df_true <- data.table::rbindlist(y_true_list) %>% as.data.frame()
  
  if (evaluate_all_files) {
    df_conf <- df_conf[1:num_samples, ]
    df_true <- df_true[1:num_samples, ]
  }
  
  zero_var_col <- apply(df_true, 2, stats::var) == 0
  if (sum(zero_var_col) > 0) {
    warning_message <- paste(sum(zero_var_col), "columns contain just one label and will be removed from evaluation")
    df_conf <- df_conf[ , !zero_var_col]
    df_true <- df_true[ , !zero_var_col]
    warning(warning_message)
  }
  
  ## auc with pROC package
  #auc_list <- purrr::map(1:ncol(df_conf), ~pROC::roc(df_true[ , .x], df_conf[ , .x], quiet = TRUE))
  
  # auc and auprc with PRROC package
  auc_list <- purrr::map(1:ncol(df_conf), ~PRROC::roc.curve(
    scores.class0 = df_conf[ , .x],
    weights.class0 = df_true[ , .x]))
  auc_vector <- vector("numeric", ncol(df_true))
  for (i in 1:length(auc_vector)) {
    auc_vector[i] <- auc_list[[i]]$auc
  }
  
  auprc_list <- purrr::map(1:ncol(df_conf), ~PRROC::pr.curve(
    scores.class0 = df_conf[ , .x],
    weights.class0 = df_true[ , .x]))
  auprc_vector <- vector("numeric", ncol(df_true))
  for (i in 1:length(auprc_vector)) {
    auprc_vector[i] <- auprc_list[[i]]$auc.integral
  }
  
  output_auc <- list(mean_auc = mean(auc_vector), median_auc = median(auc_vector), 
                     summary = summary(auc_vector), standard_deviation = sd(auc_vector))
  
  if (return_auc_vector) {
    output_auc$auc_vector <- auc_vector
  }
  
  output_auprc <- list(mean_auprc = mean(auprc_vector), median_auprc = median(auprc_vector), 
                       summary = summary(auprc_vector), standard_deviation = sd(auprc_vector))
  
  if (return_auprc_vector) {
    output_auprc$auprc_vector <- auprc_vector
  }
  
  output <- list(AUC = output_auc, AUPRC = output_auprc)
  
  if (plot_roc) {
    roc_plot <- pROC::ggroc(data = auc_list, linetype = 1, size = 0.2) +
      theme_minimal() + theme(legend.position = "none") + 
      geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="blue", linetype="dashed") + 
      scale_colour_manual(values = rep("black", length(auc_list)))
    ggsave(plot = roc_plot, filename =  roc_plot_path)
  } 
  
  return(output) 
}

#' Compute F1-score from confusion matrix.
#' 
#' @param conf_mat A confusion matrix. Columns should 
#' correspond to true label and rows to predictions.
#' @export
f1_from_conf_matrix <- function(conf_mat) {
  conf_mat <- as.matrix(conf_mat)
  if (nrow(conf_mat) != nrow(conf_mat)) {
    stop("Confusion matrix needs to have same number of rows and columns.")
  }
  precision_per_class <- diag(conf_mat)/rowSums(conf_mat)
  precision_per_class[is.na(precision_per_class)] <- 0
  recall_per_class <- diag(conf_mat)/colSums(conf_mat)  
  recall_per_class[is.na(recall_per_class)] <- 0
  f1_per_class <- vector("numeric", ncol(conf_mat))
  for (i in 1:ncol(conf_mat)) {
    f1_per_class[i] <- (2 * precision_per_class[i] * recall_per_class[i])/
                           (precision_per_class[i] + recall_per_class[i])
  }
  f1_per_class[is.na(f1_per_class)] <- 0
  balanced_f1 <- sum(colSums(conf_mat) * f1_per_class)/sum(conf_mat)
  return(list(f1_per_class = f1_per_class, balanced_f1 = balanced_f1))
}
