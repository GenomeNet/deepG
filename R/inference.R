#' Get cell states of semi-redundant chunks
#'
#' @param model.path path to keras model in hdf5 format
#' @param x semi-redundant chunks (one-hot)
#' @param maxlen time steps to unroll for
#' @param batch.size how many samples are trained in parallel
#' @param run.name name of output files without ending
#' @param type save-type, will save as hdf5 if type is set to 'hdf5' (default .csv)
#' @param verbose TRUE/FALSE
#' @export
getStates <- function(model.path,
                      x,
                      maxlen = 30,
                      batch.size = 100,
                      run.name = "output",
                      type = "csv",
                      verbose = F) {
  
  stopifnot(maxlen > 0)
  stopifnot(batch.size > 0)
  
  model <- keras::load_model_hdf5(model.path)
  # Remove the last 2 layers
  keras::pop_layer(model)
  keras::pop_layer(model)
  states <- predict(model, x, batch_size = batch.size)
  # we dont have predictions in the beginning so create some empty cell response (set it zero)
  states.begining <- states[1:maxlen,] * 0
  states.final <- rbind(states.begining, states)
  # save states as hdf5
  if (verbose) message("Saving states ...")
  if (type == "hdf5") {
    file <- hdf5r::H5File$new(paste0(run.name, "_states.hdf5"), mode = "a")
    file.grp <- hdf5r::file.h5$create_group("states1")
    file.grp <- states.final
    hdf5r::h5close(file)
  } else {
    write.table(states.final,
                file = paste0(run.name, "_states.csv"),
                sep = ";", quote = F, col.names = F,
                row.names = F)
    
  }
  return(states.final)
}

#' Get cell states from a fasta file
#'
#' @param model keras model
#' @param fasta.path path to fasta file
#' @param maxlen time steps to unroll for
#' @param batch.size how many subsequences are predicted in parallel
#' @param verbose print output
#' @export
getStatesFromFasta <- function(model = NULL,
                               fasta.path = "example_files/fasta/a.fasta",
                               maxlen = 80,
                               batch.size = 100,
                               verbose = F) {
  if (verbose)
    message("Preprocess ...")  
  # prepare fasta
  preprocessed <- deepG::preprocessFasta(fasta.path,
                                         maxlen = maxlen,
                                         vocabulary = c("\n", "a", "c", "g", "t"))
  batch.num <- 1
  batch.start <- 1
  batch.end <- batch.start + batch.size
  states <- NULL
  states <- list()
  while (batch.start < nrow(preprocessed$X)) {
    if ((batch.start + batch.size) > nrow(preprocessed$X)) {
      if (verbose)
        message("Reduce batch.size temporarily")
      batch.end <- nrow(preprocessed$X)
      # reduced batch size
    }
    if (verbose)
      message(
        paste(
          "Generating batch number",
          batch.num,
          batch.start,
          "-",
          batch.end
        ))
    x.batch <-
      preprocessed$X[batch.start:batch.end, , ] # dim should be (batch_size, length, words)
    states[[batch.num]] <- keras::predict_on_batch(model, x.batch)
    # update batch index
    batch.num <- batch.num + 1 
    batch.start <- batch.end + 1
    batch.end <- batch.start + batch.size
  }
  states.matrix <- do.call(rbind, states)
  return(states.matrix)
}

#' Write states to h5 or csv file
#' 
#' \code{writeStatesToH5} Removes layers (optional) from pretrained model and calculates states of fasta file, writes states to h5/csv file.
#' Function combines fasta entries in file to one sequence. This means predictor sequences can contain elements from more than one fasta entry. 
#' h5 file also contains sequence and positions of targets corresponding to states.   
#'     
#' @param model.path Path to a pretrained model.
#' @param layer.depth Depth of layer to evaluate. If NULL last layer is used. Can only be used for sequential model. 
#' Use \code{layer_name} for non-sequential model-
#' @param layer_name Name of layer to get output from.
#' @param fasta.path Path to fasta file.
#' @param sequence Character string, ignores fasta.path if argument given.
#' @param round_digits Number of decimal places. 
#' @param batch.size Number of samples to evaluate at once. Does not change output, only relevant for speed and memory.  
#' @param step Frequency of sampling steps.
#' @param filename Filename to store states in.
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as 0-vector.
#' @param returnStates Logical scalar, return states matrix.
#' @param verbose Whether to print model before and after removing layers.
#' @param padding Logical scalar, generate states for first maxlen nucleotides by
#' padding beginning of sequence with 0-vectors.
#' @param file_type Either "h5" or "csv".
#' @param model A keras model. If model and model.path are not NULL, model will be used for inference. 
#' @param mode Either "lm" for language model or "label" for label classification.
#' @param target_middle Whether target to predict is in middle of sequence.
#' @param format Either "fasta" or "fastq".
#' @export
writeStates <- function(model.path = NULL, layer.depth = NULL, layer_name = NULL, sequence = NULL, fasta.path = NULL, round_digits = 2,
                        filename = "states.h5", step = 1, vocabulary = c("a", "c", "g", "t"), batch.size = 256, verbose = TRUE,
                        returnStates = FALSE, padding = FALSE, file_type = "h5", model = NULL, mode = "lm", target_middle = FALSE,
                        format = "fasta") {
  
  stopifnot(mode %in% c("lm", "label"))
  stopifnot(file_type %in% c("h5", "csv"))
  stopifnot(batch.size > 0)
  stopifnot(is.null(layer.depth) | layer.depth > 0)
  stopifnot(!file.exists(paste0(filename, ".", file_type)) & !file.exists(filename))
  tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), vocabulary)
  
  # load model and sequence/file
  if (is.null(model)) {
    model <- keras::load_model_hdf5(model.path, compile = FALSE)
    model <- keras::load_model_weights_hdf5(model, model.path)
  }
  if (!is.null(sequence) && (!missing(sequence) & sequence != "")) {
    seq <- sequence %>% stringr::str_to_lower()
  } else {
    if (format == "fasta") {
      fasta.file <- microseq::readFasta(fasta.path)
    }
    if (format == "fastq") {
      fasta.file <- microseq::readFastq(fasta.path)
    }
    seq <- paste(fasta.file$Sequence, collapse = "") %>% stringr::str_to_lower()
  }
  
  # extract maxlen
  if (!target_middle) {
    maxlen <- model$input$shape[[2]] 
  } else {
    maxlen_1 <- model$input[[1]]$shape[[2]] 
    maxlen_2 <- model$input[[2]]$shape[[2]]
    maxlen <- maxlen_1 + maxlen_2
  }
  
  if (padding) {
    seq <- paste0(paste(rep("0", maxlen), collapse = ""), seq)
  }
  
  # start of samples
  if (mode == "lm") {
    numberPotentialSamples <- ceiling((nchar(seq) - maxlen)/step)   
    start_indices <- (0:(numberPotentialSamples - 1) * step) + 1
  } else {
    numberPotentialSamples <- ceiling((nchar(seq) - maxlen + 1)/step)   
    start_indices <- (0:(numberPotentialSamples - 1) * step) + 1
  }
  
  if (!is.null(layer_name)) {
    check_layer_name(model, layer_name)
    model <- tensorflow::tf$keras$Model(model$input, model$get_layer(layer_name)$output)
    cat("Computing output for model at layer", layer_name,  "\n")
    print(model)
  } else {
    
    if (!is.null(layer.depth)) {
      num_of_layers <- length(model$layers)
      if (is.null(layer.depth)) layer.depth <- num_of_layers
      stopifnot(num_of_layers > layer.depth)
      if (verbose) {
        model
        cat("Original model has", num_of_layers, "layers \n")
      }
      
      # remove layers
      if (num_of_layers - layer.depth !=0) {
        for (i in 1:(num_of_layers - layer.depth)) {
          keras::pop_layer(model)
        }
      }
    }
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
  
  # add file ending if not given
  if (!stringr::str_detect(filename, pattern = paste0("\\.", file_type, "$"))) {
    filename <- paste0(filename, ".", file_type)
  }
  
  if (file_type == "h5") {
    # create h5 file to store states
    h5_file <- hdf5r::H5File$new(filename, mode = "w") 
    if (!missing(fasta.path) & !(is.null(fasta.path))) h5_file[["fasta_file"]] <- fasta.path
    
    if (mode == "lm") {  
      if (padding) {
        if (target_middle) {
          h5_file[["target_positions"]] <- start_indices - maxlen_2 + 1
        } else {
          h5_file[["target_positions"]] <- start_indices
        }
      } else {
        if (target_middle) {
          h5_file[["target_positions"]] <- start_indices + maxlen_1 
        } else {
          h5_file[["target_positions"]] <- start_indices + maxlen
        }
      }
    } else {
      if (padding) {
        h5_file[["sample_end_position"]] <- start_indices - 1
      } else {
        h5_file[["sample_end_position"]] <- start_indices + maxlen - 1
      }
    }
    
    if (!padding) {
      h5_file[["sequence"]] <- seq
    } else {
      h5_file[["sequence"]] <- substr(seq, maxlen + 1, nchar(seq))
    }
    h5_file[["states"]] <- array(0, dim = c(0, layer.size)) 
    writer <- h5_file[["states"]]
  }
  
  numberOfBatches <- ceiling(length(start_indices)/batch.size)
  
  row <- 1
  for (i in 1:numberOfBatches) {
    # last entry might be shorter than batch.size
    if (i == numberOfBatches) {
      subsetStartInd <- start_indices[((i-1) * batch.size + 1) : length(start_indices)]
      if (mode == "lm") {
        subsetSeq <- tokSeq[subsetStartInd[1]: (subsetStartInd[length(subsetStartInd)] + maxlen)]
      } else {
        subsetSeq <- tokSeq[subsetStartInd[1]: (subsetStartInd[length(subsetStartInd)] + maxlen - 1)]
      }
    } else {
      subsetStartInd <- start_indices[((i-1) * batch.size + 1) : (i*batch.size)]
      subsetSeq <- tokSeq[subsetStartInd[1] : (subsetStartInd[length(subsetStartInd)] + maxlen)]
    }
    if (mode == "lm") {
      input <- sequenceToArray(sequence = subsetSeq, maxlen = maxlen, vocabulary = vocabulary, startInd = subsetStartInd, 
                               target_middle = target_middle)[[1]]
    } else {
      input <- sequenceToArrayLabel(sequence = subsetSeq, maxlen = maxlen, vocabulary = vocabulary, startInd = subsetStartInd)
    }
    if (target_middle) {
      current_batch_size <- dim(input[[1]])[1]
    } else {
      current_batch_size <- dim(input)[1]
    }
    activations <- keras::predict_on_batch(object = model, x = input)
    activations <- as.array(activations)
    # some layers give predictions for every char in sequence, 
    # since return_sequences = TRUE, filter last prediction 
    if (length(dim(activations)) == 3) {
      activations <- activations[ , maxlen, ]
    }
    # last entry might be shorter
    if (i == numberOfBatches) {
      numRemainingSamples <- length(start_indices) - row + 1
      activation_batch <- matrix(round(activations[1:numRemainingSamples, ], digits = round_digits), ncol = ncol(activations))
      if (file_type == "h5") { 
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
      if (file_type == "h5") {
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
  if (returnStates & (file_type == "h5")) states <- writer[ , ]
  if (returnStates & (file_type == "csv")) states <- read.table(filename, header = TRUE)
  if (file_type == "h5") h5_file$close_all()
  if (returnStates) return(states)
}

#' Write states to h5 file
#' 
#' @description \code{writeStatesByFastaEntries} Removes layers (optional) from pretrained model and calculates states of fasta file, writes a separate 
#' h5 file for every fasta entry in fasta file. h5 files also contain the nucleotide sequence and positions of targets corresponding to states.      
#' Names of output files are: file_path + "Nr" + i + filename + file_type, where i is the number of the fasta entry. 
#' 
#' @param model.path Path to a pretrained model.
#' @param layer.depth Depth of layer to evaluate. If NULL last layer is used. 
#' @param layer_name Name of layer to get output from.
#' @param fasta.path Path to fasta file.
#' @param round_digits Number of decimal places. 
#' @param file_name Filename to store states, function adds "Nr" + "i" before name, where i is entry number. 
#' @param file_path Path to folder, where to write output.
#' @param step Frequency of sampling steps.
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as 0-vector.
#' @param batch.size Number of samples to evaluate at once. Does not change output, only relevant for speed and memory.  
#' @param padding Logical scalar, generate states for first maxlen nucleotides by
#' padding beginning of sequence with 0-vectors.
#' @param file_type Either "h5" or "csv".
#' @param model A keras model. If model and model.path are not NULL, model will be used for inference.
#' @param mode Either "lm" for language model or "label" for label classification.
#' @param target_middle Whether target to predict is in middle of sequence.
#' @param format Either "fasta" or "fastq".
#' @export
writeStatesByFastaEntries <- function(model.path = NULL, layer.depth = NULL, layer_name = NULL, fasta.path, round_digits = 2, 
                                      file_name = "states.h5", file_path, step = 1, vocabulary = c("a", "c", "g", "t"),
                                      batch.size = 256, padding = FALSE, file_type = "h5", model = NULL, mode = "lm",
                                      target_middle = FALSE, format = "fasta") {
  
  stopifnot(mode %in% c("lm", "label"))
  
  if (is.null(model)) {
    model <- keras::load_model_hdf5(model.path, compile = FALSE)
    model <- keras::load_model_weights_hdf5(model, model.path)
  }
  
  # extract maxlen
  if (!target_middle) {
    maxlen <- model$input$shape[[2]] 
  } else {
    maxlen_1 <- model$input[[1]]$shape[[2]] 
    maxlen_2 <- model$input[[2]]$shape[[2]]
    maxlen <- maxlen_1 + maxlen_2
  }
  
  # load fasta file
  if (format == "fasta") {
    fasta.file <- microseq::readFasta(fasta.path)
  }
  if (format == "fastq") {
    fasta.file <- microseq::readFastq(fasta.path)
  }
  df <- fasta.file[ , c("Sequence", "Header")]
  names(df) <- c("seq", "header") 
  rownames(df) <- NULL
  
  for (i in 1:nrow(df)) {
    verbose <- FALSE
    if (i == 1) verbose <- TRUE 
    # skip entry if too short 
    if ((nchar(df[i, "seq"]) < maxlen + 1) & !padding) next
    writeStates(model.path = model.path, layer.depth = layer.depth, layer_name = layer_name, sequence = df[i, "seq"], 
                round_digits = round_digits, fasta.path = fasta.path,
                filename = paste0(file_path, "/Nr", as.character(i), "_", file_name),
                step = step, vocabulary = vocabulary, batch.size = batch.size,
                verbose = verbose, padding = padding, file_type = file_type, mode = mode,
                target_middle = target_middle, model = model)  
  }
}

#' Write states to h5 file
#' 
#' @description \code{writeStatesByFastaEntries} Removes layers (optional) from pretrained model and calculates states of fasta file,
#' writes separate states matrix in one .h5 file for every fasta entry. 
#' h5 file also contains the nucleotide sequences and positions of targets corresponding to states.         
#' 
#' To acces the content of h5 file: 
#' h5_path <- "/path/to/file"
#' h5_file <- hdf5r::H5File$new(h5_path, mode = "r")
#' a <- h5_file[["states"]]
#' # shows header names
#' # names(a)
#' b <- a[["someHeaderName"]]
#' #shows state matrix 
#' head(b[,])
#' # get correspondig sequence
#' hdf5r::h5attr(b, "sequence")
#' #show corresponding positions 
#' #hdf5r::h5attr(b, "target_positions") # for language model
#' #hdf5r::h5attr(b, "sample_end_position") # for label classification
#' h5_file$close_all()
#'   
#' 
#' @param model.path Path to a pretrained model.
#' @param layer.depth Depth of layer to evaluate. If NULL last layer is used. 
#' @param fasta.path Path to fasta file.
#' @param round_digits Number of decimal places. 
#' @param step Frequency of sampling steps.
#' @param h5.filename Filename of h5 file to store states.
#' @param vocabulary Vector of allowed characters, character outside vocabulary get encoded as 0-vector.
#' @param batch.size Number of samples to evaluate at once. Does not change output, only relevant for speed and memory.  
#' @param layer_name Name of layer to get output from.
#' @param padding Logical scalar, generate states for first maxlen nucleotides by
#' padding beginning of sequence with 0-vectors. 
#' @param verbose Whether to print model before and after removing layers.
#' @param model A keras model. If model and model.path are not NULL, model will be used for inference. 
#' @param mode Either "lm" for language model or "label" for label classification.
#' @param target_middle Whether target to predict is in middle of sequence.
#' @param format Either "fasta" or "fastq".
#' @export
statesByFastaOneFile <- function(model.path, layer.depth = NULL, fasta.path, round_digits = 2, h5.filename = "states.h5",
                                 step = 1,  vocabulary = c("a", "c", "g", "t"), batch.size = 256, layer_name = NULL,
                                 padding = FALSE, verbose = TRUE, model = NULL, mode = "lm", target_middle = FALSE, 
                                 format = "fasta") {
  
  stopifnot(mode %in% c("lm", "label"))
  stopifnot(batch.size > 0)
  stopifnot(is.null(layer.depth) | layer.depth > 0)
  stopifnot(!file.exists(paste0(h5.filename,".h5")) & !file.exists(h5.filename))
  
  if (is.null(model)) {
    model <- keras::load_model_hdf5(model.path, compile = FALSE)
    model <- keras::load_model_weights_hdf5(model, model.path)
  }
  
  # extract maxlen
  if (!target_middle) {
    maxlen <- model$input$shape[[2]] 
  } else {
    maxlen_1 <- model$input[[1]]$shape[[2]] 
    maxlen_2 <- model$input[[2]]$shape[[2]]
    maxlen <- maxlen_1 + maxlen_2
  }
  
  if (!is.null(layer_name)) {
    check_layer_name(model, layer_name)
    model <- tensorflow::tf$keras$Model(model$input, model$get_layer(layer_name)$output)
    cat("Computing output for model at layer", layer_name,  "\n")
    print(model)
  } else {
    
    if (!is.null(layer.depth)) {
      num_of_layers <- length(model$layers)
      if (is.null(layer.depth)) layer.depth <- num_of_layers
      stopifnot(num_of_layers > layer.depth)
      if (verbose) {
        model
        cat("Original model has", num_of_layers, "layers \n")
      }
      
      # remove layers
      if (num_of_layers - layer.depth !=0) {
        for (i in 1:(num_of_layers - layer.depth)) {
          keras::pop_layer(model)
        }
      }
    }
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
  
  # read fasta file and create separate row in data frame per entry
  if (format == "fasta") {
    fasta.file <- microseq::readFasta(fasta.path)
  }
  if (format == "fastq") {
    fasta.file <- microseq::readFastq(fasta.path)
  }
  df <- as.data.frame(fasta.file$Sequence)
  names(df) <- "seq"
  df$length <- nchar(fasta.file$Sequence)
  df$name <- fasta.file$Header
  rownames(df) <- NULL
  df$seq <- stringr::str_to_lower(df$seq) 
  df$name <- trimws(df$name)
  
  if (padding) {
    for (i in 1:nrow(df)) {
      df[i, "seq"] <- paste0(paste(rep("0", maxlen), collapse = ""), df[i, "seq"])
      df[i, "length"] <- df[i, "length"] + maxlen
    }
  }
  
  # tokenize sequences
  tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), vocabulary)
  tokSeqComplete <- keras::texts_to_sequences(tokenizer, paste(df$seq, collapse = ""))[[1]] - 1
  
  # create h5 file to store states
  if (!stringr::str_detect(h5.filename, pattern = "\\.h5$")) {
    h5.filename <- paste0(h5.filename,".h5")
  }
  file.h5 <- hdf5r::H5File$new(h5.filename, mode = "w") 
  if (!missing(fasta.path)) file.h5[["fasta_path"]] <- fasta.path
  states.grp <- file.h5$create_group("states")
  
  # loop over fasta entries
  for (j in 1:nrow(df)) {
    seq <- df[j, "seq"]
    seq_length <- df[j, "length"]
    seq_name <- df[j, "name"]
    if (seq_name == "") seq_name <- paste0("dummyHeader_", as.character(j))
    if (j == 1) {
      start <- 1
    } else {
      start <- cumsum(df$length[j - 1])
    }
    end <- start + df$length[j] - 1
    tokSeq <- tokSeqComplete[start:end] 
    
    # start of samples
    if (mode == "lm") {
      numberPotentialSamples <- ceiling((nchar(seq) - maxlen)/step)   
      start_indices <- (0:(numberPotentialSamples - 1) * step) + 1
    } else {
      numberPotentialSamples <- ceiling((nchar(seq) - maxlen + 1)/step)   
      start_indices <- (0:(numberPotentialSamples - 1) * step) + 1
    }
    
    states.grp[[seq_name]] <- array(0, dim = c(0, layer.size)) 
    writer <- states.grp[[seq_name]]
    
    numberOfBatches <- ceiling(length(start_indices)/batch.size)
    
    row <- 1
    for (i in 1:numberOfBatches) {
      
      # last entry might be shorter than batch.size
      if (i == numberOfBatches) {
        subsetStartInd <- start_indices[((i-1) * batch.size + 1) : length(start_indices)]
        if (mode == "lm") {
          subsetSeq <- tokSeq[subsetStartInd[1]: (subsetStartInd[length(subsetStartInd)] + maxlen)]
        } else {
          subsetSeq <- tokSeq[subsetStartInd[1]: (subsetStartInd[length(subsetStartInd)] + maxlen - 1)]
        }
      } else {
        subsetStartInd <- start_indices[((i-1) * batch.size + 1) : (i*batch.size)]
        subsetSeq <- tokSeq[subsetStartInd[1] : (subsetStartInd[length(subsetStartInd)] + maxlen)]
      }
      if (mode == "lm") {
        input <- sequenceToArray(sequence = subsetSeq, maxlen = maxlen, vocabulary = vocabulary, startInd = subsetStartInd,
                                 target_middle = target_middle)[[1]]
      } else {
        input <- sequenceToArrayLabel(sequence = subsetSeq, maxlen = maxlen, vocabulary = vocabulary, startInd = subsetStartInd)
      }
      current_batch_size <- dim(input)[1]
      activations <- keras::predict_on_batch(object = model, x = input)
      activations <- as.array(activations)
      # some layers give predictions for every char in sequence, 
      # since return_sequences = TRUE, filter last prediction 
      if (length(dim(activations)) == 3) {
        activations <- activations[ , maxlen, ]
      }
      # last entry might be shorter
      if (i == numberOfBatches) {
        numRemainingSamples <- length(start_indices) - row + 1
        writer[row:(row + numRemainingSamples - 1), ] <- round(activations[1:numRemainingSamples, ],
                                                               digits = round_digits)
      } else {
        writer[row:(row + current_batch_size - 1), ] <- round(activations, digits = round_digits)
        row <- row + current_batch_size
      }
    }
    
    if (mode == "lm") {
      if (padding) {
        if (target_middle) {
          hdf5r::h5attr(states.grp[[seq_name]], "target_positions") <- start_indices - maxlen_2  
        } else {
          hdf5r::h5attr(states.grp[[seq_name]], "target_positions") <- start_indices  
        }
      } else {
        if (target_middle) {
          hdf5r::h5attr(states.grp[[seq_name]], "target_positions") <- start_indices + maxlen_1
        } else {
          hdf5r::h5attr(states.grp[[seq_name]], "target_positions") <- start_indices + maxlen
        }
      }
    } else {
      if (padding) {
        hdf5r::h5attr(states.grp[[seq_name]], "sample_end_position") <- start_indices - 1  
      } else {
        hdf5r::h5attr(states.grp[[seq_name]], "sample_end_position") <- start_indices + maxlen - 1
      }
    }
    
    if (!padding) {
      hdf5r::h5attr(states.grp[[seq_name]], "sequence") <- seq
    } else {
      hdf5r::h5attr(states.grp[[seq_name]], "sequence") <- substr(seq, maxlen + 1, nchar(seq))
    }
  }
  file.h5$close_all()
}

#' Inference on set of fasta files
#' 
#' @description Wrapper function for inference functions \code{\link{{writeStates}}, \code{\link{{writeStatesByFastaEntries}}
#' or \code{\link{{statesByFastaOneFile}}. These functions compute states for a single fasta file, \code{statesWrapper} extends functions 
#' to compute states on a set of fasta files.    
#' 
#' Output files in \{output_path} will be named by adding "_states" + ".h5/.csv" to name of fasta file. If \code{function_name} is "writeStatesByFastaEntries",
#' number of the fasta entry will be added additionally for every fasta file. 
#'  
#' @param function_name Name of base function, either "writeStates", "writeStatesByFastaEntries" or "statesByFastaOneFile".
#' @param function_args List of arguments for function in \code{function_name}. If no new argument is given, default values will be used.
#' @param fasta.path Path to folder containing fasta files.
#' @param output_path Path to folder of output files.
#' @examples 
#' \dontrun{
#' statesWrapper(function_name = "writeStates", function_args = list(model.path = "/path/to/model.h5", padding = TRUE, layer_depth = 1), 
#' fasta.path = "/path/to/fasta/folder", output_path = "/where/to/write/output")
#' }
#' @export
statesWrapper <- function(function_name = "writeStates", function_args = list(), fasta.path, output_path) {
  
  fasta_files <- list.files(
    path = xfun::normalize_path(fasta.path),
    pattern = paste0("\\.", "fasta", "$"),
    full.names = TRUE)
  stopifnot(length(fasta_files) > 0)
  
  fasta_names <- list.files(
    path = xfun::normalize_path(fasta.path),
    pattern = paste0("\\.", "fasta", "$"),
    full.names = FALSE)
  fasta_names <- stringr::str_remove(fasta_names, ".fasta$")
  
  # load model and remove layers before looping through files
  if (is.null(function_args$model)) {
    model <- keras::load_model_hdf5(function_args$model.path, compile = FALSE)
    model <- keras::load_model_weights_hdf5(model, function_args$model.path)
    function_args$model <- model
  } else {
    model <- function_args$model
  }
  
  num_of_layers <- length(model$layers)
  if (is.null(function_args$layer.depth)) function_args$layer.depth <- num_of_layers
  stopifnot(num_of_layers >= function_args$layer.depth)
  
  # remove layers
  if (is.null(function_args$layer_name)) {
    if (num_of_layers - function_args$layer.depth !=0) {
      for (i in 1:(num_of_layers - function_args$layer.depth)) {
        keras::pop_layer(model)
      }
    }
  } else {
    model <- tensorflow::tf$keras$Model(model$input, model$get_layer(function_args$layer_name)$output)
  }
  
  # layers have been removed before calling base function, thus set layer.depth to NULL (don't remove additional layers)  
  function_args$layer.depth <- NULL
  function_args$model <- model
  
  new_arguments <- names(function_args)
  default_arguments <- formals(function_name)
  
  # overwrite default arguments 
  for (arg in new_arguments) {
    default_arguments[arg] <- function_args[arg] 
  }
  
  if (function_name == "writeStates") {
    formals(writeStates) <- default_arguments
  }
  if (function_name == "writeStatesByFastaEntries") {
    formals(writeStatesByFastaEntries) <- default_arguments
  } 
  if (function_name == "statesByFastaOneFile") {
    formals(statesByFastaOneFile) <- default_arguments
  } 
  
  for (i in 1:length(fasta_files)) {
    
    if (function_name == "writeStates") {
      filename <- paste0(output_path, "/", fasta_names[i], "_states")
      writeStates(fasta.path = fasta_files[i], filename = filename, returnStates = FALSE, verbose = FALSE)
    } 
    
    if (function_name == "writeStatesByFastaEntries") {
      filename <- paste0(fasta_names[i], "_states")
      writeStatesByFastaEntries(file_name = filename, file_path = output_path, fasta.path = fasta_files[i], verbose = FALSE)
    } 
    
    if (function_name == "statesByFastaOneFile") {
      filename <- paste0(output_path, "/", fasta_names[i], "_states")
      statesByFastaOneFile(fasta.path = fasta_files[i], h5.filename = filename, verbose = FALSE)
    } 
  }
}

#' Get states for label classification model    
#'
#' Computes output at certain model layer. Forces every fasta entry to have length maxlen by either padding sequences shorter than maxlen or taking random subsequence for
#' longer sequences.  
#' 
#' @inheritParams writeStates
#' @export
predReads <- function(model.path = NULL, layer.depth = NULL, layer_name = NULL, fasta.path, round_digits = 2, format = "fasta",
                      filename = "states.h5", vocabulary = c("a", "c", "g", "t"), batch.size = 256, verbose = TRUE,
                      returnStates = FALSE, model = NULL) {
  
  file_type <- "h5"
  stopifnot(batch.size > 0)
  stopifnot(is.null(layer.depth) | layer.depth > 0)
  stopifnot(!file.exists(paste0(filename, ".", file_type)) & !file.exists(filename))
  tokenizer <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "N"), vocabulary)
  
  # load model and sequence/file
  if (is.null(model)) {
    model <- keras::load_model_hdf5(model.path)
  }
  # get maxlen
  maxlen <- model$input$shape[[2]] 
  
  if (format == "fasta") {
    fasta.file <- microseq::readFasta(fasta.path)
  }
  if (format == "fastq") {
    fasta.file <- microseq::readFastq(fasta.path)
  }
  nucSeq <- as.character(fasta.file$Sequence)
  seq_length <- nchar(fasta.file$Sequence)
  # TODO: check if loop is slow   
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
  seq <- paste(nucSeq, collapse = "") %>% stringr::str_to_lower()
  
  if (!is.null(layer_name)) {
    model <- tensorflow::tf$keras$Model(model$input, model$get_layer(layer_name)$output)
    cat("Computing output for model at layer", layer_name,  "\n")
    print(model)
  } else {
    
    if (!is.null(layer.depth)) {
      num_of_layers <- length(model$layers)
      if (is.null(layer.depth)) layer.depth <- num_of_layers
      stopifnot(num_of_layers >= layer.depth)
      if (verbose) {
        model
        cat("Original model has", num_of_layers, "layers \n")
      }
      
      # remove layers
      if (num_of_layers - layer.depth !=0) {
        for (i in 1:(num_of_layers - layer.depth)) {
          keras::pop_layer(model)
        }
      }
    }
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
  
  # add file ending if not given
  if (!stringr::str_detect(filename, pattern = paste0("\\.", file_type, "$"))) {
    filename <- paste0(filename, ".", file_type)
  }
  
  if (file_type == "h5") {
    # create h5 file to store states
    h5_file <- hdf5r::H5File$new(filename, mode = "w") 
    if (!missing(fasta.path)) h5_file[["fasta_file"]] <- fasta.path
    
    h5_file[["header_names"]] <- fasta.file$Header  
    h5_file[["sequences"]] <- nucSeq
    h5_file[["states"]] <- array(0, dim = c(0, layer.size)) 
    writer <- h5_file[["states"]]
  }
  
  # create temporary fasta file for generator
  temp_file <- tempfile(pattern = "", fileext = paste0(".", format))
  fasta.file <- fasta.file[1, ]
  fasta.file$Sequence <- seq
  fasta.file$Header <- "one_seq"
  
  if (format == "fasta") {
    microseq::writeFasta(fdta = fasta.file, out.file = temp_file)
  }
  if (format == "fastq") {
    microseq::writeFastq(fdta = fasta.file, out.file = temp_file)
  }
  gen <- fastaLabelGenerator(corpus.dir = temp_file, 
                             format = "fasta",
                             batch.size = batch.size,
                             maxlen = maxlen,
                             vocabulary = vocabulary,
                             step = maxlen, 
                             labelVocabulary = c("one_seq"),
                             reverseComplements = FALSE)
  
  numberOfBatches <- ceiling(length(nucSeq)/batch.size)
  
  row <- 1
  if (numberOfBatches > 1) {
    for (i in 1:(numberOfBatches - 1)) {
      activations <- keras::predict_generator(model, generator = gen, steps = 1)
      writer[row : (row + batch.size - 1), ] <- activations
      row <- row + batch.size
    } 
  }
  # last batch might be shorter
  activations <- keras::predict_generator(model, generator = gen, steps = 1)
  writer[row : length(nucSeq), ] <- activations[1 : length(row:length(nucSeq)), ]
  
  file.remove(temp_file)
  if (returnStates & (file_type == "h5")) states <- writer[ , ]
  if (file_type == "h5") h5_file$close_all()
  if (returnStates) return(states)
}


#' Read states from h5
#' 
#' Reads rows from h5 file created by  \code{\link{{writeStatesByFastaEntries}} or  \code{\link{{writeStates}},
#' rows correspond to position in sequence and columns to neurons.
#' 
#' @param h5_path Path to h5 file.
#' @param rows Range of rows to read.
#' @param complete Return all entries if TRUE.
#' @param getTargetPositions Return position of corresponding targets if TRUE.
#' @export
readRowsFromH5 <- function(h5_path, rows, verbose = TRUE, complete = FALSE, getTargetPositions = FALSE) {
  h5_file <- hdf5r::H5File$new(h5_path, mode = "r")
  read_states <- h5_file[["states"]]
  
  if (getTargetPositions) {
    train_mode <- ifelse("sample_end_position" %in% names(h5_file), "label", "lm")
    if (train_mode == "label") {
      read_targetPos <- h5_file[["sample_end_position"]]
    } else {
      read_targetPos <- h5_file[["target_positions"]]
    }
  }
  
  if (verbose) {
    cat("states matrix has", dim(read_states[ , ])[1], "rows and ",  dim(read_states[ , ])[2], "columns \n")
  }
  if (complete) {
    states <- read_states[ , ]
    if (getTargetPositions) {
      targetPos <- read_targetPos[ ]
    }
  } else {
    states <- read_states[rows, ]
    if (getTargetPositions) {
      targetPos <- read_targetPos[rows]
    }
  }
  h5_file$close_all()
  if (getTargetPositions) {
    if (train_mode == "label") {
      return(list(states = states, sample_end_position = targetPos))
    } else {
      return(list(states = states, targetPos = targetPos))
    }
  } else {
    return(states)
  }
}

