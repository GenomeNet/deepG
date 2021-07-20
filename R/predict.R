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
#' @param model.path Path to pretrained model.
#' @param model A keras model. 
#' @param batch.size Number of samples per batch.
#' @param step How often to take a sample.
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as 0-vector.
#' @param label_vocabulary Labels for targets. Equal to vocabulary if not given.
#' @param numberOfBatches How many batches to evaluate.
#' @param filePath Where to store output, if missing output won't be written.
#' @param format File format, "fasta" or "fastq".
#' @param filename Name of output file.
#' @param plot Returns density plot of accuracies if TRUE.
#' @param mode Either "lm" for language model and "label_header", "label_csv" or "label_folder" for label classification.
#' @param acc_per_batch Whether to return vector with accuracies for every batch.
#' @inheritParams fastaFileGenerator
#' @param ambiguous_nuc How to handle nucleotides outside vocabulary, either "zero", "discard" or "equal". If "zero", input gets encoded as zero vector; 
#' if "equal" input is 1/length(vocabulary) x length(vocabulary). If "discard" samples containing nucleotides outside vocabulary get discarded.        
#' @param evaluate_all_files Boolean, if TRUE will iterate over all files in \code{fasta.path} once. \code{numberOfBatches} will be overwritten. 
#' @export
evaluateFasta <- function(fasta.path,
                          #model.path = NULL,
                          model = NULL,
                          batch.size = 100,
                          step = 1,
                          vocabulary = c("a", "c", "g", "t"),
                          label_vocabulary = c("a", "c", "g", "t"),
                          numberOfBatches = 10,
                          filePath = NULL,
                          format = "fasta",
                          filename = "",
                          target_middle = FALSE,
                          plot = FALSE, 
                          mode = "lm",
                          acc_per_batch = FALSE, 
                          output_format = "target_right",
                          ambiguous_nuc = "zero",
                          evaluate_all_files = FALSE,
                          verbose = TRUE,
                          max_iter = 20000,
                          target_from_csv = NULL) {
  
  model.path <- NULL 
  stopifnot(mode %in% c("lm", "label_header", "label_folder", "label_csv"))
  
  if (is.null(label_vocabulary)) label_vocabulary <- vocabulary
  
  if (!is.null(model.path)) {
    model <- keras::load_model_hdf5(model.path)
  }
  
  if (length(fasta.path) > 1) {
    numberOfBatches <- rep(ceiling(numberOfBatches/length(fasta.path)), length(fasta.path))
  }  
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
      for (file in files) {
        if (format == "fasta") {
          seq_vector <- microseq::readFasta(file)$Sequence
          if (mode == "lm") {
            seq_vector <- seq_vector[nchar(seq_vector) >= (maxlen + 1)]
          } else {
            seq_vector <- seq_vector[nchar(seq_vector) >=  maxlen]
          }
        } else {
          seq_vector <- microseq::readFastq(file)$Sequence
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
        num_samples[i] <- num_samples[i] + new_samples
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
                              showWarnings = FALSE,
                              shuffleFastaEntries = FALSE,
                              reverseComplements = FALSE,
                              output_format = output_format,
                              ambiguous_nuc = ambiguous_nuc)
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
                               showWarnings = FALSE,
                               shuffleFastaEntries = FALSE,
                               labelVocabulary = label_vocabulary,
                               reverseComplements = FALSE,
                               ambiguous_nuc = ambiguous_nuc,
                               target_from_csv = target_from_csv)
  }
  
  
  acc <- vector("numeric")
  loss_per_class <- vector("list", length = length(fasta.path))
  confMat <- matrix(0, nrow = length(label_vocabulary), ncol = length(label_vocabulary))
  batch_index <- 1
  start_time <- Sys.time()
  ten_percent_steps <- 1/(10:1) *  overall_num_batches
  percentage_index <- 1
  
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
                                    reverseComplements = FALSE,
                                    numTargets = length(fasta.path),
                                    onesColumn = k,
                                    ambiguous_nuc = ambiguous_nuc,
                                    padding = FALSE)
    }
    
    for (i in 1:numberOfBatches[k]) {
      z <- gen()
      x <- z[[1]]
      y <- z[[2]]
      
      y_conf <- predict(model, x, verbose = 0)
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
      
      df_true_pred <- data.frame(
        true = factor(y_true, levels = 1:(length(label_vocabulary)), labels = label_vocabulary),
        pred = factor(y_pred, levels = 1:(length(label_vocabulary)), labels = label_vocabulary)
      )
      
      cce_loss_new <- tensorflow::tf$keras$losses$categorical_crossentropy(y, y_conf)
      cce_loss[[i]] <- cce_loss_new$numpy()
      
      if (acc_per_batch) {
        acc[i] <- sum(y_pred == y_true)/length(y_true)
      }
      
      cm <- yardstick::conf_mat(df_true_pred, true, pred)
      confMat <- confMat + cm[[1]]
      
      if (verbose & (batch_index == 10)) {
        time_passed <- as.double(difftime(Sys.time(), start_time, units = "hours"))
        time_estimation <- (overall_num_batches/10) * time_passed
        cat("Evaluation will take approximately", round(time_estimation, 3), "hours \n")
      }
      
      if (verbose & (batch_index > ten_percent_steps[percentage_index])) {
        cat("Progress: ", percentage_index * 10 ,"% \n")  
        percentage_index <- percentage_index + 1 
      }
      
    }
    loss_per_class[[k]] <- unlist(cce_loss)
  }
  
  if (plot & acc_per_batch) {
    df <- data.frame(accuracies = acc)
    acc_plot <- ggplot2::ggplot(df) + ggplot2::geom_density(ggplot2::aes(x = accuracies), color = "blue", alpha = 0.1, fill = "blue")
    print(acc_plot)
  }
  
  if (!acc_per_batch) {
    acc <- sum(diag(confMat))/sum(confMat)
  } 
  
  loss <- mean(unlist(loss_per_class))
  for (i in 1:length(loss_per_class)) {
    loss_per_class[[i]] <- mean(loss_per_class[[i]])
  }
  loss_per_class <- unlist(loss_per_class)
  m <- as.matrix(confMat)
  class_acc <- vector("numeric")
  for (i in 1:ncol(m)) {
    class_acc[i] <- m[i, i]/sum(m[ , i])
  }
  names(class_acc) <- label_vocabulary
  macro_acc <- mean(class_acc)
  
  # if (!is.null(filePath)) {
  #   saveRDS(acc, file = paste0(filePath, "/", filename, "_Acc.rds"))
  #   saveRDS(confMat, file = paste0(filePath, "/", filename, "_ConfMat.rds"))
  #   saveRDS(loss, file = paste0(filePath, "/", filename, "_Loss.rds"))
  #   saveRDS(macro_acc, file = paste0(filePath, "/", filename, "_MacroAcc.rds"))
  #   # if (plot){
  #   #   ggplot2::ggsave(acc_plot, filename = paste0(filePath, "/", filename, "_AccPlot.pdf" ))
  #   # }
  # }
  
  if (mode == "label_folder") {
    
    names(loss_per_class) <- label_vocabulary
    
    output_list <- list(
      accuracy = acc,
      confusion_matrix = confMat,
      loss = loss,
      macro_acc = macro_acc,
      class_acc = class_acc,
      loss_per_class = loss_per_class
    ) 
  } else {
    
    output_list <- list(
      accuracy = acc,
      confusion_matrix = confMat,
      loss = loss
    ) 
  }
  
  if (!is.null(filePath)) {
    saveRDS(output_list, file = paste0(filePath, "/", filename, "_outputList.rds"))
  }
  
  return(output_list)
}
