#' Make prediction for nucleotide sequence or entries in fasta/fastq file
#'
#' @description Removes layers (optional) from pretrained model and calculates states of fasta/fastq file or nucleotide sequence.
#' Writes states to h5/csv file. To acces the content of h5 output use \code{load_predictions} function.
#' If "one_seq" computes prediction for sequence argument or fasta/fastq file.
#' Combines fasta entries in file to one sequence. This means predictor sequences can contain elements from more than one fasta entry.
#' If "by_entry" will output a separate file for each fasta/fastq entry.
#' Names of output files are: output_dir + "Nr" + i + filename + output_type, where i is the number of the fasta entry.
#' If "by_entry_one_file", will store prediction for all fasta entry in one h5 file.
#' If "one_pred_per_entry" will make one prediction for each entry by either picking random sample for long sequences
#' or pad sequence for short sequences.
#'
#' @param output_format Either "one_seq", "by_entry", "by_entry_one_file", "one_pred_per_entry".
#' @param output_type "h5" or "csv". If output_format is "by_entries_one_file", "one_pred_per_entry" can only be "h5".
#' @param return_states Return predictions as data frame. Only supported for output_format "one_seq".
#' @inheritParams predict_model_one_seq
#' @inheritParams predict_model_by_entry
#' @inheritParams predict_model_by_entry_one_file
#' @inheritParams predict_model_one_pred_per_entry
#' @examples
#' # make prediction for single sequence and write to h5 file
#' model <- create_model_lstm_cnn(maxlen = 20, layer_lstm = 8, layer_dense = 2, verbose = FALSE)
#' vocabulary <- c("a", "c", "g", "t")
#' sequence <- paste(sample(vocabulary, 200, replace = TRUE), collapse = "")
#' output_file <- tempfile(fileext = ".h5")
#' predict_model(output_format = "one_seq", model = model, step = 10,
#'              sequence = sequence, filename = output_file, mode = "label")
#' 
#' # make prediction for fasta file with multiple entries, write output to separate h5 files
#' fasta_path <- tempfile(fileext = ".fasta")
#' create_dummy_data(file_path = fasta_path, num_files = 1,
#'                  num_seq = 5, seq_length = 100,
#'                  write_to_file_path = TRUE)
#' model <- create_model_lstm_cnn(maxlen = 20, layer_lstm = 8, layer_dense = 2, verbose = FALSE)
#' output_dir <- tempfile()
#' dir.create(output_dir)
#' predict_model(output_format = "by_entry", model = model, step = 10, verbose = FALSE,
#'                output_dir = output_dir, mode = "label", path_input = fasta_path)
#' list.files(output_dir)
#' @export
predict_model <- function(output_format = "one_seq", model = NULL, layer_name = NULL, sequence = NULL, path_input = NULL,
                          round_digits = 2, filename = "states.h5", step = 1, vocabulary = c("a", "c", "g", "t"),
                          batch_size = 256, verbose = TRUE, return_states = FALSE, padding = FALSE, output_type = "h5",
                          path_model = NULL, mode = "lm", lm_format = "target_right", output_dir = NULL,
                          format = "fasta", include_seq = TRUE, reverse_complement_encoding = FALSE) {
  
  stopifnot(output_format %in% c("one_seq", "by_entry", "by_entry_one_file", "one_pred_per_entry"))
  if (output_format %in% c("by_entry_one_file", "one_pred_per_entry") & output_type == "csv") {
    message("by_entry_one_file or one_pred_per_entry only implemented for h5 output.
            Setting output_type to h5")
    output_type <- "h5"
  }
  
  # # add file ending if not given
  # if (!stringr::str_detect(filename, pattern = paste0("\\.", output_type, "$"))) {
  #   filename <- paste0(filename, ".", output_type)
  # }
  
  if (output_format == "one_seq") {
    output_list <- predict_model_one_seq(path_model = path_model, layer_name = layer_name, sequence = sequence, path_input = path_input,
                                         round_digits = round_digits, filename = filename, step = step, vocabulary = vocabulary,
                                         batch_size = batch_size, verbose = verbose, return_states = return_states, padding = padding,
                                         output_type = output_type, model = model, mode = mode, lm_format = lm_format,
                                         format = format, include_seq = include_seq, reverse_complement_encoding = reverse_complement_encoding)
    return(output_list)
  }
  
  if (output_format == "by_entry") {
    predict_model_by_entry(path_model = path_model, layer_name = layer_name, path_input = path_input,
                           round_digits = round_digits, filename = filename, step = step, vocabulary = vocabulary,
                           batch_size = batch_size, verbose = verbose, padding = padding,
                           output_type = output_type, model = model, mode = mode, lm_format = lm_format,
                           output_dir = output_dir, format = format, include_seq = include_seq,
                           reverse_complement_encoding = reverse_complement_encoding)
  }
  
  if (output_format == "by_entry_one_file") {
    predict_model_by_entry_one_file(path_model = path_model, layer_name = layer_name, path_input = path_input,
                                    round_digits = round_digits, filename = filename, step = step, vocabulary = vocabulary,
                                    batch_size = batch_size, verbose = verbose, padding = padding,
                                    model = model, mode = mode, lm_format = lm_format, format = format,
                                    include_seq = include_seq, reverse_complement_encoding = reverse_complement_encoding)
  }
  
  if (output_format == "one_pred_per_entry") {
    if (mode == "lm") {
      stop("one_pred_per_entry only implemented for label classification")
    }
    predict_model_one_pred_per_entry(path_model = path_model, layer_name = layer_name, path_input = path_input,
                                     round_digits = round_digits, filename = filename, vocabulary = vocabulary,
                                     batch_size = batch_size, verbose = verbose, model = model, format = format,
                                     reverse_complement_encoding = reverse_complement_encoding)
  }
  
}


#' Write output of specific model layer to h5 or csv file.
#'
#' Removes layers (optional) from pretrained model and calculates states of fasta file, writes states to h5/csv file.
#' Function combines fasta entries in file to one sequence. This means predictor sequences can contain elements from more than one fasta entry.
#' h5 file also contains sequence and positions of targets corresponding to states.
#'
#' @param path_model Path to a pretrained model.
#' @param layer_name Name of layer to get output from. If NULL, will use the last layer.
#' @param path_input Path to fasta file.
#' @param sequence Character string, ignores path_input if argument given.
#' @param round_digits Number of decimal places.
#' @param batch_size Number of samples to evaluate at once. Does not change output, only relevant for speed and memory.
#' @param step Frequency of sampling steps.
#' @param filename Filename to store states in. No file output if argument is NULL.
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as 0-vector.
#' @param return_states Logical scalar, return states matrix.
#' @param verbose Whether to print model before and after removing layers.
#' @param padding Logical scalar, generate states for first maxlen nucleotides by
#' padding beginning of sequence with 0-vectors.
#' @param output_type Either "h5" or "csv".
#' @param model A keras model. If model and path_model are not NULL, model will be used for inference.
#' @param mode Either "lm" for language model or "label" for label classification.
#' @param format Either "fasta" or "fastq".
#' @param include_seq Whether to include input sequence in h5 file.
#' @keywords internal
predict_model_one_seq <- function(path_model = NULL, layer_name = NULL, sequence = NULL, path_input = NULL, round_digits = 2,
                                  filename = "states.h5", step = 1, vocabulary = c("a", "c", "g", "t"), batch_size = 256, verbose = TRUE,
                                  return_states = FALSE, padding = FALSE, output_type = "h5", model = NULL, mode = "lm", lm_format = "target_right",
                                  format = "fasta", include_seq = TRUE, reverse_complement_encoding = FALSE) {
  
  stopifnot(mode %in% c("lm", "label"))
  stopifnot(output_type %in% c("h5", "csv"))
  file_output <- !is.null(filename)
  if (!file_output) {
    if (!return_states) stop("If filename is NULL, return_states must be TRUE; otherwise function produces no output.")
    filename <- tempfile(fileext = paste0(".", output_type))
  }
  stopifnot(batch_size > 0)
  stopifnot(!file.exists(paste0(filename, ".", output_type)) & !file.exists(filename))
  if (reverse_complement_encoding) {
    test_len <- length(vocabulary) != 4
    if (test_len || all(sort(stringr::str_to_lower(vocabulary)) != c("a", "c", "g", "t"))) {
      stop("reverse_complement_encoding only implemented for A,C,G,T vocabulary")
    }
  }
  tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), vocabulary)
  
  # load model and sequence/file
  if (is.null(model)) {
    model <- keras::load_model_hdf5(path_model, compile = FALSE)
  }
  
  if (is.null(layer_name)) {
    layer_name <- model$output_names
    message(paste("layer_name not specified. Using layer", layer_name))
  }
  
  if (!is.null(sequence) && (!missing(sequence) & sequence != "")) {
    seq <- sequence %>% stringr::str_to_lower()
  } else {
    if (format == "fasta") {
      fasta.file <- microseq::readFasta(path_input)
    }
    if (format == "fastq") {
      fasta.file <- microseq::readFastq(path_input)
    }
    if (nrow(fasta.file) > 1) {
      text_1 <- paste("Your file has", nrow(fasta.file), "entries. 'one_seq'  output_format will concatenate them to a single sequence.\n")
      text_2 <- "Use 'by_entry' or 'by_entry_one_file' output_format to evaluate them separately."
      message(paste0(text_1, text_2))
    }
    seq <- paste(fasta.file$Sequence, collapse = "") %>% stringr::str_to_lower()
  }
  
  # extract maxlen
  target_middle <- ifelse(mode == "lm" && (lm_format %in% c("target_middle_lstm", "target_middle_cnn")), TRUE, FALSE)
  if (!target_middle) {
    if (reverse_complement_encoding) {
      model$input[[1]]$shape[[2]]
    } else {
      maxlen <- model$input$shape[[2]]
    }
  } else {
    maxlen_1 <- model$input[[1]]$shape[[2]]
    maxlen_2 <- model$input[[2]]$shape[[2]]
    maxlen <- maxlen_1 + maxlen_2
  }
  
  if (padding) {
    seq <- paste0(paste(rep("0", maxlen), collapse = ""), seq)
  }
  
  if (nchar(seq) < maxlen) {
    stop("Input sequence is shorter than maxlen. Set padding = TRUE to pad sequence to bigger size.")
  }
  
  # start of samples
  if (mode == "lm") {
    numberPotentialSamples <- ceiling((nchar(seq) - maxlen)/step)
    start_indices <- (0:(numberPotentialSamples - 1) * step) + 1
  } else {
    numberPotentialSamples <- ceiling((nchar(seq) - maxlen + 1)/step)
    start_indices <- (0:(numberPotentialSamples - 1) * step) + 1
  }
  
  check_layer_name(model, layer_name)
  model <- tensorflow::tf$keras$Model(model$input, model$get_layer(layer_name)$output)
  if (verbose) {
    cat("Computing output for model at layer", layer_name,  "\n")
    print(model)
  }
  
  # extract number of neurons in last layer
  if (length(model$output$shape$dims) == 3) {
    if (!("lstm" %in% stringr::str_to_lower(model$output_names))) {
      stop("Output dimension of layer is > 1, format not supported yet")
    }
    layer.size <- model$output$shape[[3]]
  } else {
    layer.size <- model$output$shape[[2]]
  }
  
  # tokenize sequence
  tokSeq <- stringr::str_to_lower(seq)
  tokSeq <- keras::texts_to_sequences(tokenizer, seq)[[1]] - 1
  
  if (mode == "lm") {
    if (padding) {
      if (target_middle) {
        pos_arg <- start_indices - maxlen_2 + 1
      } else {
        pos_arg <- start_indices
      }
    } else {
      if (target_middle) {
        pos_arg <- start_indices + maxlen_1
      } else {
        pos_arg <- start_indices + maxlen
      }
    }
  } else {
    if (padding) {
      pos_arg <- start_indices - 1
    } else {
      pos_arg <- start_indices + maxlen - 1
    }
  }
  
  if (include_seq) {
    if (!padding) {
      output_seq <- seq
    } else {
      output_seq <- substr(seq, maxlen + 1, nchar(seq))
    }
  }
  
  if (output_type == "h5") {
    # create h5 file to store states
    h5_file <- hdf5r::H5File$new(filename, mode = "w")
    h5_file[["multi_entries"]] <- FALSE
    #if (!missing(path_input) & !(is.null(path_input))) h5_file[["fasta_file"]] <- path_input
    if (mode == "lm") {
      h5_file[["target_position"]] <- pos_arg
    } else {
      h5_file[["sample_end_position"]] <- pos_arg
    }
    if (include_seq) h5_file[["sequence"]] <- output_seq
    h5_file[["states"]] <- array(0, dim = c(0, layer.size))
    writer <- h5_file[["states"]]
  }
  
  number_batches <- ceiling(length(start_indices)/batch_size)
  
  row <- 1
  for (i in 1:number_batches) {
    # last entry might be shorter than batch_size
    if (i == number_batches) {
      subsetStartInd <- start_indices[((i-1) * batch_size + 1) : length(start_indices)]
      if (mode == "lm") {
        subsetSeq <- tokSeq[subsetStartInd[1]: (subsetStartInd[length(subsetStartInd)] + maxlen)]
      } else {
        subsetSeq <- tokSeq[subsetStartInd[1]: (subsetStartInd[length(subsetStartInd)] + maxlen - 1)]
      }
    } else {
      subsetStartInd <- start_indices[((i-1) * batch_size + 1) : (i*batch_size)]
      subsetSeq <- tokSeq[subsetStartInd[1] : (subsetStartInd[length(subsetStartInd)] + maxlen)]
    }
    
    if (mode == "lm") {
      input <- seq_encoding_lm(sequence = subsetSeq, maxlen = maxlen, vocabulary = vocabulary, start_ind = subsetStartInd,
                               output_format = lm_format, adjust_start_ind = TRUE)[[1]]
    } else {
      input <- seq_encoding_label(sequence = subsetSeq, maxlen = maxlen, vocabulary = vocabulary, start_ind = subsetStartInd,
                                  adjust_start_ind = TRUE)
    }
    
    if (target_middle) {
      current_batch_size <- dim(input[[1]])[1]
    } else {
      current_batch_size <- dim(input)[1]
    }
    
    if (reverse_complement_encoding) {
      x_1 <- input
      x_2 <- array(x_1[ , (dim(x_1)[2]):1, 4:1], dim = dim(x_1))
      input <- list(x_1, x_2)
    }
    
    activations <- keras::predict_on_batch(object = model, x = input)
    activations <- as.array(activations)
    # some layers give predictions for every char in sequence,
    # since return_sequences = TRUE, filter last prediction
    if (length(dim(activations)) == 3) {
      activations <- activations[ , maxlen, ]
    }
    # last entry might be shorter
    if (i == number_batches) {
      numRemainingSamples <- length(start_indices) - row + 1
      activation_batch <- matrix(round(activations[1:numRemainingSamples, ], digits = round_digits), ncol = ncol(activations))
      if (output_type == "h5") {
        writer[row:(row + numRemainingSamples - 1), ] <- activation_batch
      } else {
        # add column with target position for lm or end position for label classification
        if (mode == "lm") {
          if (padding) {
            target_position <- subsetStartInd
          } else {
            target_position <- subsetStartInd + maxlen
          }
        } else {
          if (padding) {
            target_position <- subsetStartInd - 1
          } else {
            target_position <- subsetStartInd + maxlen - 1
          }
        }
        activation_batch <- cbind(activation_batch, target_position)
        if (i == 1) {
          if (mode == "lm") {
            col.names <- c(as.character(1:layer.size), "target position")
          } else {
            col.names <- c(as.character(1:layer.size), "sample_end_position")
          }
          write.table(x = activation_batch, file = filename, col.names = col.names, row.names = FALSE, append = FALSE)
        } else {
          write.table(x = activation_batch, file = filename, col.names = FALSE, row.names = FALSE, append = TRUE)
        }
      }
    } else {
      activation_batch <- round(activations, digits = round_digits)
      if (output_type == "h5") {
        writer[row:(row + current_batch_size - 1), ] <-  activation_batch
      } else {
        # add column with target position for lm or end position for label classification
        if (mode == "lm") {
          if (padding) {
            if (target_middle) {
              target_position <- subsetStartInd - maxlen_2 + 1
            } else {
              target_position <- subsetStartInd
            }
          } else {
            if (target_middle) {
              target_position <- subsetStartInd + maxlen_1
            } else {
              target_position <- subsetStartInd + maxlen
            }
          }
        } else {
          if (padding) {
            target_position <- subsetStartInd - 1
          } else {
            target_position <- subsetStartInd + maxlen - 1
          }
        }
        activation_batch <- cbind(activation_batch, target_position)
        if (i == 1) {
          if (mode == "lm") {
            col.names <- c(as.character(1:layer.size), "target position")
          } else {
            col.names <- c(as.character(1:layer.size), "sample_end_position")
          }
          write.table(x = activation_batch, file = filename, col.names = col.names, row.names = FALSE, append = FALSE)
        } else {
          write.table(x = activation_batch, file = filename, col.names = FALSE, row.names = FALSE, append = TRUE)
        }
      }
      row <- row + current_batch_size
    }
  }
  if (return_states & (output_type == "h5")) states <- writer[ , ]
  if (return_states & (output_type == "csv")) states <- read.table(filename, header = TRUE)
  if (output_type == "h5") h5_file$close_all()
  if (!file_output) file.remove(filename)
  
  if (return_states) {
    output_list <- list()
    output_list$states <- states
    if (mode == "lm") output_list$target_position <- pos_arg
    if (mode == "label") output_list$sample_end_position <- pos_arg
    if (include_seq) output_list$sequence <- output_seq
    return(output_list)
  } else {
    return(NULL)
  }
}

#' Write states to h5 file
#'
#' @description Removes layers (optional) from pretrained model and calculates states of fasta file, writes a separate
#' h5 file for every fasta entry in fasta file. h5 files also contain the nucleotide sequence and positions of targets corresponding to states.
#' Names of output files are: file_path + "Nr" + i + filename + output_type, where i is the number of the fasta entry.
#'
#' @param filename Filename to store states, function adds "_nr_" + "i" after name, where i is entry number.
#' @param output_dir Path to folder, where to write output.
#' @keywords internal
predict_model_by_entry <- function(path_model = NULL, layer_name = NULL, path_input, round_digits = 2,
                                   filename = "states.h5", output_dir, step = 1, vocabulary = c("a", "c", "g", "t"),
                                   batch_size = 256, padding = FALSE, output_type = "h5", model = NULL, mode = "lm",
                                   lm_format = "target_right", format = "fasta", reverse_complement_encoding = FALSE,
                                   verbose = FALSE, include_seq = FALSE) {
  
  stopifnot(mode %in% c("lm", "label"))
  if (endsWith(filename, paste0(".", output_type))) {
    filename <- stringr::str_remove(filename, paste0(".", output_type, "$"))
    filename <- basename(filename)
  }
  
  if (is.null(model)) {
    model <- keras::load_model_hdf5(path_model, compile = FALSE)
    model <- keras::load_model_weights_hdf5(model, path_model)
  }
  
  # extract maxlen
  target_middle <- ifelse(mode == "lm" && (lm_format %in% c("target_middle_lstm", "target_middle_cnn")), TRUE, FALSE)
  if (!target_middle) {
    if (reverse_complement_encoding) {
      model$input[[1]]$shape[[2]]
    } else {
      maxlen <- model$input$shape[[2]]
    }
  } else {
    maxlen_1 <- model$input[[1]]$shape[[2]]
    maxlen_2 <- model$input[[2]]$shape[[2]]
    maxlen <- maxlen_1 + maxlen_2
  }
  
  # load fasta file
  if (format == "fasta") {
    fasta.file <- microseq::readFasta(path_input)
  }
  if (format == "fastq") {
    fasta.file <- microseq::readFastq(path_input)
  }
  df <- fasta.file[ , c("Sequence", "Header")]
  names(df) <- c("seq", "header")
  rownames(df) <- NULL
  
  for (i in 1:nrow(df)) {
    
    if (i > 1) verbose <- FALSE
    # skip entry if too short
    if ((nchar(df[i, "seq"]) < maxlen + 1) & !padding) next
    predict_model_one_seq(path_model = path_model, layer_name = layer_name, sequence = df[i, "seq"],
                          round_digits = round_digits, path_input = path_input,
                          filename = paste0(output_dir, "/", filename, "_nr_", as.character(i), ".", output_type),
                          step = step, vocabulary = vocabulary, batch_size = batch_size,
                          verbose = verbose, padding = padding, output_type = output_type, mode = mode,
                          lm_format = lm_format, model = model, include_seq = include_seq,
                          reverse_complement_encoding = reverse_complement_encoding)
  }
}

#' Write states to h5 file
#'
#' @description Removes layers (optional) from pretrained model and calculates states of fasta file,
#' writes separate states matrix in one .h5 file for every fasta entry.
#' h5 file also contains the nucleotide sequences and positions of targets corresponding to states.
#' @keywords internal
predict_model_by_entry_one_file <- function(path_model, path_input, round_digits = 2, filename = "states.h5",
                                            step = 1,  vocabulary = c("a", "c", "g", "t"), batch_size = 256, layer_name = NULL,
                                            padding = FALSE, verbose = TRUE, model = NULL, mode = "lm", lm_format = lm_format,
                                            format = "fasta", include_seq = TRUE, reverse_complement_encoding = FALSE) {
  
  stopifnot(mode %in% c("lm", "label"))
  # if (endsWith(filename, paste0(".", output_type))) {
  #   filename <- stringr::str_remove(filename, paste0(".", output_type, "$"))
  # }
  
  if (is.null(model)) {
    model <- keras::load_model_hdf5(path_model, compile = FALSE)
    model <- keras::load_model_weights_hdf5(model, path_model)
  }
  
  target_middle <- ifelse(mode == "lm" && (lm_format %in% c("target_middle_lstm", "target_middle_cnn")), TRUE, FALSE)
  # extract maxlen
  if (!target_middle) {
    if (reverse_complement_encoding) {
      model$input[[1]]$shape[[2]]
    } else {
      maxlen <- model$input$shape[[2]]
    }
  } else {
    maxlen_1 <- model$input[[1]]$shape[[2]]
    maxlen_2 <- model$input[[2]]$shape[[2]]
    maxlen <- maxlen_1 + maxlen_2
  }
  
  # extract number of neurons in last layer
  if (length(model$output$shape$dims) == 3) {
    if (!("lstm" %in% stringr::str_to_lower(model$output_names))) {
      stop("Output dimension of layer is > 1, format not supported yet")
    }
    layer.size <- model$output$shape[[3]]
  } else {
    layer.size <- model$output$shape[[2]]
  }
  
  # load fasta file
  if (format == "fasta") {
    fasta.file <- microseq::readFasta(path_input)
  }
  if (format == "fastq") {
    fasta.file <- microseq::readFastq(path_input)
  }
  df <- fasta.file[ , c("Sequence", "Header")]
  names(df) <- c("seq", "header")
  rownames(df) <- NULL
  
  if (verbose) {
    # check if names are unique
    if (length(df$header) != length(unique(df$header))) {
      message("Header names are not unique, adding '_header_x' to names (x being the header number)")
      df$header <- paste0(df$header, paste0("_header_", 1:length(df$header)))
    }
  }
  
  # create h5 file to store states
  
  h5_file <- hdf5r::H5File$new(filename, mode = "w")
  h5_file[["multi_entries"]] <- TRUE
  states.grp <- h5_file$create_group("states")
  if (mode == "lm") {
    target_position.grp <- h5_file$create_group("target_position")
  } else {
    sample_end_position.grp <- h5_file$create_group("sample_end_position")
  }
  if (include_seq) seq.grp <- h5_file$create_group("sequence")
  
  for (i in 1:nrow(df)) {
    
    seq_name <- df$header[i]
    temp_file <- tempfile(fileext = ".h5")
    if (i > 1) verbose <- FALSE
    # skip entry if too short
    if ((nchar(df[i, "seq"]) < maxlen + 1) & !padding) next
    output_list <- predict_model_one_seq(path_model = path_model, layer_name = layer_name, sequence = df$seq[i], path_input = path_input,
                                         round_digits = round_digits, filename = temp_file, step = step, vocabulary = vocabulary,
                                         batch_size = batch_size, verbose = verbose, return_states = TRUE, padding = padding,
                                         output_type = "h5", model = model, mode = mode, lm_format = lm_format,
                                         format = format, include_seq = include_seq, reverse_complement_encoding = reverse_complement_encoding)
    
    
    states.grp[[seq_name]] <- output_list$states
    if (mode == "lm") {
      target_position.grp[[seq_name]] <- output_list$target_position
    } else {
      sample_end_position.grp[[seq_name]] <- output_list$sample_end_position
    }
    if (include_seq) seq.grp[[seq_name]] <- output_list$sequence
  }
  h5_file$close_all()
}

#' Get states for label classification model.
#'
#' Computes output at specified model layer. Forces every fasta entry to have length maxlen by either padding sequences shorter than maxlen or taking random subsequence for
#' longer sequences.
#'
#' @inheritParams predict_model_one_seq
#' @keywords internal

predict_model_one_pred_per_entry <- function(model = NULL, layer_name = NULL, path_input, round_digits = 2, format = "fasta",
                                             filename = "states.h5", vocabulary = c("a", "c", "g", "t"), batch_size = 256, verbose = TRUE,
                                             return_states = FALSE, path_model = NULL, reverse_complement_encoding = FALSE, include_seq = FALSE) {
  
  file_type <- "h5"
  stopifnot(batch_size > 0)
  stopifnot(!file.exists(paste0(filename, ".", file_type)) & !file.exists(filename))
  # token for ambiguous nucleotides
  for (i in letters) {
    if (!(i %in% stringr::str_to_lower(vocabulary))) {
      amb_nuc_token <- i
      break
    }
  }
  tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "N"), vocabulary)
  
  if (is.null(layer_name)) {
    layer_name <- model$output_names
    message(paste("layer_name not specified. Using layer", layer_name))
  }
  
  # load model and sequence/file
  if (is.null(model)) {
    model <- keras::load_model_hdf5(path_model)
  }
  
  # extract maxlen
  if (reverse_complement_encoding) {
    model$input[[1]]$shape[[2]]
  } else {
    maxlen <- model$input$shape[[2]]
  }
  
  if (format == "fasta") {
    fasta.file <- microseq::readFasta(path_input)
  }
  if (format == "fastq") {
    fasta.file <- microseq::readFastq(path_input)
  }
  
  num_samples <- nrow(fasta.file)
  
  nucSeq <- as.character(fasta.file$Sequence)
  seq_length <- nchar(fasta.file$Sequence)
  
  for (i in 1:length(nucSeq)) {
    # take random subsequence
    if (seq_length[i] > maxlen) {
      start <- sample(1 : (seq_length[i] - maxlen + 1) , size = 1)
      nucSeq[i] <- substr(nucSeq[i], start = start, stop = start + maxlen - 1)
    }
    # pad sequence
    if (seq_length[i] < maxlen) {
      nucSeq[i] <- paste0(paste(rep("N", maxlen - seq_length[i]), collapse = ""), nucSeq[i])
    }
  }
  
  model <- tensorflow::tf$keras$Model(model$input, model$get_layer(layer_name)$output)
  cat("Computing output for model at layer", layer_name,  "\n")
  print(model)
  
  # extract number of neurons in last layer
  if (length(model$output$shape$dims) == 3) {
    if (!("lstm" %in% stringr::str_to_lower(model$output_names))) {
      stop("Output dimension of layer is > 1, format not supported yet")
    }
    layer.size <- model$output$shape[[3]]
  } else {
    layer.size <- model$output$shape[[2]]
  }
  
  if (file_type == "h5") {
    # create h5 file to store states
    h5_file <- hdf5r::H5File$new(filename, mode = "w")
    if (!missing(path_input)) h5_file[["fasta_file"]] <- path_input
    
    h5_file[["header_names"]] <- fasta.file$Header
    if (include_seq) h5_file[["sequences"]] <- nucSeq
    h5_file[["states"]] <- array(0, dim = c(0, layer.size))
    h5_file[["multi_entries"]] <- FALSE
    writer <- h5_file[["states"]]
  }
  
  rm(fasta.file)
  #nucSeq <- paste(nucSeq, collapse = "") %>% stringr::str_to_lower()
  number_batches <- ceiling(num_samples/batch_size)
  if (verbose) cat("Evaluating", number_batches, ifelse(number_batches > 1, "batches", "batch"), "\n")
  row <- 1
  string_start_index <- 1 
  ten_percent_steps <- seq(number_batches/10, number_batches, length.out = 10)
  percentage_index <- 1
  
  if (number_batches > 1) {
    for (i in 1:(number_batches - 1)) {
      string_end_index <-string_start_index + batch_size - 1 
      char_seq <- nucSeq[string_start_index : string_end_index] %>% paste(collapse = "") 
      if (i == 1) start_ind <- seq(1, nchar(char_seq), maxlen)
      one_hot_batch <- seq_encoding_label(sequence = NULL, maxlen = maxlen, vocabulary = vocabulary,
                                          start_ind = start_ind, ambiguous_nuc = "zero", nuc_dist = NULL,
                                          use_quality = FALSE, quality_vector = NULL, use_coverage = FALSE, max_cov = NULL,
                                          cov_vector = NULL, n_gram = NULL, n_gram_stride = 1, char_sequence = char_seq,
                                          tokenizer = tokenizer, adjust_start_ind = TRUE) 
      if (reverse_complement_encoding) one_hot_batch <- list(one_hot_batch, reverse_complement_tensor(one_hot_batch))
      activations <- keras::predict_on_batch(model, one_hot_batch)
      writer[row : (row + batch_size - 1), ] <- activations
      row <- row + batch_size
      string_start_index <- string_end_index + 1
      
      if (verbose & (i > ten_percent_steps[percentage_index]) & percentage_index < 10) {
        cat("Progress: ", percentage_index * 10 ,"% \n")
        percentage_index <- percentage_index + 1
      }
      
    }
  }
  
  # last batch might be shorter
  char_seq <- nucSeq[string_start_index : length(nucSeq)] %>% paste(collapse = "") 
  one_hot_batch <- seq_encoding_label(sequence = NULL, maxlen = maxlen, vocabulary = vocabulary,
                                      start_ind = seq(1, nchar(char_seq), maxlen), ambiguous_nuc = "zero", nuc_dist = NULL,
                                      use_quality = FALSE, quality_vector = NULL, use_coverage = FALSE, max_cov = NULL,
                                      cov_vector = NULL, n_gram = NULL, n_gram_stride = 1, char_sequence = char_seq,
                                      tokenizer = tokenizer, adjust_start_ind = TRUE) 
  if (reverse_complement_encoding) one_hot_batch <- list(one_hot_batch, reverse_complement_tensor(one_hot_batch))
  activations <- keras::predict_on_batch(model, one_hot_batch)
  writer[row : num_samples, ] <- activations[1 : length(row:num_samples), ]
  
  if (verbose) cat("Progress: 100 % \n")
  
  if (return_states & (file_type == "h5")) states <- writer[ , ]
  if (file_type == "h5") h5_file$close_all()
  if (return_states) return(states)
}


#' Read states from h5 file
#'
#' Reads rows from h5 file created by  \code{\link{predict_model}} function,
#' rows correspond to position in sequence and columns to neurons.
#'
#' @param h5_path Path to h5 file.
#' @param rows Range of rows to read.  If NULL read everything.
#' @param get_target_positions Return position of corresponding targets if TRUE.
#' @param get_seq Return nucleotide sequence if TRUE.
#' @examples
#' # make prediction for single sequence and write to h5 file
#' model <- create_model_lstm_cnn(maxlen = 20, layer_lstm = 8, layer_dense = 2, verbose = FALSE)
#' vocabulary <- c("a", "c", "g", "t")
#' sequence <- paste(sample(vocabulary, 200, replace = TRUE), collapse = "")
#' output_file <- tempfile(fileext = ".h5")
#' predict_model(output_format = "one_seq", model = model, step = 10,
#'               sequence = sequence, filename = output_file, mode = "label")
#' load_prediction(h5_path = output_file)
#' @export
load_prediction <- function(h5_path, rows = NULL, verbose = TRUE,
                            get_target_positions = FALSE, get_seq = FALSE) {
  
  if (is.null(rows)) complete <- TRUE
  h5_file <- hdf5r::H5File$new(h5_path, mode = "r")
  
  multi_entries <- ifelse(h5_file[["multi_entries"]][], TRUE, FALSE)
  if (!multi_entries) {
    number_entries <- 1
  } else {
    entry_names <- names(h5_file[["states"]])
    number_entries <- length(entry_names)
    output_list <- list()
  }
  train_mode <- ifelse("target_position" %in% names(h5_file), "lm", "label")
  
  if (get_target_positions & !any(c("sample_end_position", "target_position") %in% names(h5_file))) {
    get_target_positions <- FALSE
    message("File does not contain target positions.")
  }
  
  if (!multi_entries) {
    
    read_states <- h5_file[["states"]]
    
    if (get_target_positions) {
      if (train_mode == "label") {
        read_targetPos <- h5_file[["sample_end_position"]]
      } else {
        read_targetPos <- h5_file[["target_position"]]
      }
    }
    
    if (verbose) {
      cat("states matrix has", dim(read_states[ , ])[1], "rows and",  dim(read_states[ , ])[2], "columns \n")
    }
    if (complete) {
      states <- read_states[ , ]
      if (get_target_positions) {
        targetPos <- read_targetPos[ ]
      }
    } else {
      states <- read_states[rows, ]
      if (get_target_positions) {
        targetPos <- read_targetPos[rows]
      }
    }
    
    if (is.null(dim(states))) {
      states <- matrix(states, nrow = 1)
    }
    
    contains_seq <- FALSE
    if (get_seq) {
      if ("sequence" %in% names(h5_file)) {
        contains_seq <- TRUE
        sequence <- h5_file[["sequence"]][]
      } else {
        contains_seq <- FALSE
        message("File does not contain sequence.")
      }
    }
    
    h5_file$close_all()
    output_list <- list(states = states)
    if (get_target_positions) {
      if (train_mode == "label") {
        output_list$sample_end_position <- targetPos
      } else {
        output_list$target_position <- targetPos
      }
    }
    
    if (get_seq && contains_seq) {
      output_list$sequence <- sequence
    }
    
    return(output_list)
    
    # multi entries
  } else {
    
    if (verbose) {
      cat("file contains", number_entries, "entries \n")
    }
    
    if (get_target_positions) {
      target_name <- ifelse("target_position" %in% names(h5_file), "target_position", "sample_end_position")
    }
    
    if (get_seq & !("sequence" %in% names(h5_file))) {
      message("File does not contain sequence.")
      get_seq <- FALSE
    }
    
    for (i in 1:number_entries) {
      
      entry_name <- entry_names[i]
      states <- h5_file[["states"]][[entry_name]][ , ]
      if (is.null(dim(states))) {
        states <- matrix(states, nrow = 1)
      }
      
      if (get_seq) {
        sequence <- h5_file[["sequence"]][[entry_name]][ ]
      }
      
      if (get_target_positions) {
        targetPos <- h5_file[[target_name]][[entry_name]][ ]
      }
      
      if (!complete) {
        states <- states[rows, ]
        
        if (is.null(dim(states))) {
          states <- matrix(states, nrow = 1)
        }
        
        if (get_target_positions) {
          targetPos <- hdf5r::h5attr(read_states, target_name)[rows]
        }
      }
      
      if (get_target_positions) {
        if (train_mode == "label") {
          l <- list(states = states, sample_end_position = targetPos)
        } else {
          l <- list(states = states, target_position = targetPos)
        }
      } else {
        l <- list(states = states)
      }
      
      if (get_seq) {
        l[["sequence"]] <- sequence
      }
      output_list[[entry_name]] <- l
    }
    h5_file$close_all()
    names(output_list) <- entry_names
    return(output_list)
  }
}

#' Create summary of predictions 
#'
#' Create summary data frame for confidence predictions over 1 or several state files or a data frame.
#' Columns in file or data frame should be confidence preditions for one class,
#' i.e. each rows should sum to 1 and have nonnegative entries. 
#' Output data frame contains average confidence scores, max score and percentage of votes for each class.
#'
#' @param states_path Folder containing state files or a single file with same ending as file_type.
#' @param label_names Names of predicted classes.
#' @param file_type "h5" or "csv".
#' @param df A states data frame. Ignores states_dir if not NULL.
#' @examples 
#' m <- matrix(runif(20*3), ncol = 3)
#' m <- m/rowSums(m)
#' label_names <- paste0("class_", 1:3)
#' df <- as.data.frame(m)
#' summarize_states(label_names = label_names, df = df)
#' @export
summarize_states <- function(states_path = NULL, label_names = NULL, file_type = "h5", df = NULL) {
  
  if (!is.null(df)) {
    states_path <- NULL
  }
  
  if (is.null(states_path)) {
    state_files <- 1
  } else {
    if (endsWith(states_path, file_type)) {
      state_files <- states_path
    } else {
      state_files <- list.files(states_path, full.names = TRUE)
    }
  }
  
  if (!is.null(label_names)) {
    num_labels <- length(label_names)
  }
  
  summary_list <- list()
  
  for (state_file in state_files) {
    
    if (is.null(df)) {
      if (file_type == "h5") {
        df <- load_prediction(h5_path = state_file, get_target_positions = FALSE, verbose = FALSE)
        df <- as.data.frame(df$states)
      }
      if (file_type == "csv") {
        df <- read.csv(state_file)
        if (ncol(df) != num_labels) {
          df <- read.csv2(state_file)
        }
      }
    }
    
    if (state_file == state_files[1] & is.null(label_names)) {
      label_names <- paste0("X_", 1:ncol(df))
      num_labels <- length(label_names)
    }
    
    stopifnot(ncol(df) == num_labels)
    
    names(df) <- c(label_names)
    
    mean_df <- data.frame(matrix(0, nrow = 1, ncol = num_labels))
    names(mean_df) <- paste0("mean_conf_", label_names)
    max_df <- data.frame(matrix(0, nrow = 1, ncol = num_labels))
    names(max_df) <- paste0("max_conf_", label_names)
    
    for (label in label_names) {
      mean_df[[paste0("mean_conf_", label)]] <-  mean(df[[label]])
      max_df[[paste0("max_conf_", label)]] <-  max(df[[label]])
    }
    
    vote_distribution <- apply(df[label_names], 1, which.max)
    vote_perc <- table(factor(vote_distribution, levels = 1:length(label_names)))/length(vote_distribution)
    votes_df <- data.frame(matrix(vote_perc, nrow = 1, ncol = num_labels))
    names(votes_df) <- paste0("vote_perc_", label_names)
    
    mean_prediction <- label_names[which.max(unlist(mean_df))]
    max_prediction <- label_names[which.max(unlist(max_df))]
    vote_prediction <- label_names[which.max(vote_perc)]
    
    if (is.null(states_path)) {
      file_name <- NA
    } else {
      file_name <- basename(state_file)
    }
    
    summary_list[[state_file]] <- data.frame(file_name, mean_df, max_df, votes_df,
                                             mean_prediction, max_prediction, vote_prediction,
                                             num_prediction = nrow(df))
    
  }
  
  summary_df <- data.table::rbindlist(summary_list)
  return(summary_df)
}
