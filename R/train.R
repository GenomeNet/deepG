#' @title Train neural network on genomic data
#'
#' @description
#' Train a neural network on genomic data. Data can be fasta/fastq files, rds files or a prepared data set.
#' If the data is given as collection of fasta, fastq or rds files, function will create a data generator that extracts training and validation batches
#' from files. Function includes several options to determine the sampling strategy of the generator and preprocessing of the data.  
#' Training progress can be visualized in tensorboard. Model weights can be stored during training using checkpoints.        
#' 
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_folder
#' @inheritParams generator_fasta_label_header_csv
#' @param train_type Either `"lm"`, `"lm_rds"` for language model; `"label_header"`, `"label_folder"`, `"label_csv"`, `"label_rds"` for classification or `"dummy_gen"`.
#' \itemize{
#' \item Language model is trained to predict character(s) in a sequence. \cr
#' \item `"label_header"`/`"label_folder"`/`"label_csv"` are trained to predict a corresponding class given a sequence as input.
#' \item If `"label_header"`, class will be read from fasta headers.
#' \item If `"label_folder"`, class will be read from folder, i.e. all files in one folder must belong to the same class. 
#' \item If `"label_csv"`, targets are read from a csv file. This file should have one column named "file". The targets then correspond to entries in that row (except "file"
#' column). Example: if we are currently working with a file called "a.fasta" and corresponding label is "label_1", there should be a row in our csv file  
#' 
#'  |  file       | label_1 | label_2 | 
#'  |   ---       |   ---   |  ---    |   
#'  | "a.fasta"   |    1    |    0    |
#'
#' \item If `"label_rds"`, generator will iterate over set of .rds files containing each a list of input and target tensors. Not implemented for model
#' with multiple inputs. 
#' \item If `"lm_rds"`, generator will iterate over set of .rds files and will split tensor according to `target_len` argument
#' (targets are last `target_len` nucleotides of each sequence). 
#' \item  If `"dummy_gen"`, generator creates random data once and repeatedly feeds these to model.
#' \item  If `"masked_lm"`, generator maskes some parts of the input. See `masked_lm` argument for details.
#' }
#' @param model A keras model.
#' @param path Path to training data. If \code{train_type} is \code{label_folder}, should be a vector or list
#' where each entry corresponds to a class (list elements can be directories and/or individual files). If \code{train_type} is not \code{label_folder}, 
#' can be a single directory or file or a list of directories and/or files.
#' @param path_val Path to validation data. See `path` argument for details.
#' @param dataset List of training data holding training samples in RAM instead of using generator. Should be list with two entries called `"X"` and `"Y"`.
#' @param dataset_val List of validation data. Should have two entries called `"X"` and `"Y"`.
#' @param path_checkpoint Path to checkpoints folder or `NULL`. If `NULL`, checkpoints don't get stored.
#' @param path_log Path to directory to write training scores. File name is `run_name` + `".csv"`. No output if `NULL`.
#' @param train_val_ratio For generator defines the fraction of batches that will be used for validation (compared to size of training data), i.e. one validation iteration
#' processes \code{batch_size} \eqn{*} \code{steps_per_epoch} \eqn{*} \code{train_val_ratio} samples. If you use dataset instead of generator and \code{dataset_val} is `NULL`, splits \code{dataset}
#' into train/validation data.
#' @param run_name Name of the run. Name will be used to identify output from callbacks. If `NULL`, will use date as run name. 
#' If name already present, will add `"_2"` to name or `"_{x+1}"` if name ends with `_x`, where `x` is some integer. 
#' @param batch_size Number of samples used for one network update.
#' @param epochs Number of iterations.
#' @param max_queue_size Maximum size for the generator queue.
#' @param reduce_lr_on_plateau Whether to use learning rate scheduler.
#' @param lr_plateau_factor Factor of decreasing learning rate when plateau is reached.
#' @param patience Number of epochs waiting for decrease in validation loss before reducing learning rate.
#' @param cooldown Number of epochs without changing learning rate.
#' @param steps_per_epoch Number of training batches per epoch.
#' @param step Frequency of sampling steps.
#' @param shuffle_file_order Boolean, whether to go through files sequentially or shuffle beforehand.
#' @param vocabulary Vector of allowed characters. Characters outside vocabulary get encoded as specified in \code{ambiguous_nuc}.
#' @param initial_epoch Epoch at which to start training. Note that network
#' will run for (\code{epochs} - \code{initial_epochs}) rounds and not \code{epochs} rounds.
#' @param path_tensorboard Path to tensorboard directory or `NULL`. If `NULL`, training not tracked on tensorboard.
#' @param save_best_only Only save model that improved on best validation loss score.
#' @param save_weights_only Whether to save weights only.
#' @param seed Sets seed for reproducible results.
#' @param shuffle_input Whether to shuffle entries in file.
#' @param tb_images Whether to show custom images (confusion matrix) in tensorboard "IMAGES" tab.
#' @param format File format, `"fasta"`, `"fastq"`, `"rds"` or `"fasta.tar.gz"`, `"fastq.tar.gz"` for `tar.gz` files. 
#' @param path_file_log Write name of files used for training to csv file if path is specified.
#' @param vocabulary_label Character vector of possible targets. Targets outside \code{vocabulary_label} will get discarded if
#' \code{train_type = "label_header"}.
#' @param file_limit Integer or `NULL`. If integer, use only specified number of randomly sampled files for training. Ignored if greater than number of files in \code{path}.
#' @param reverse_complement_encoding Whether to use both original sequence and reverse complement as two input sequences.
#' @param output_format Determines shape of output tensor for language model.
#' Either `"target_right"`, `"target_middle_lstm"`, `"target_middle_cnn"` or `"wavenet"`.
#' Assume a sequence `"AACCGTA"`. Output correspond as follows
#' \itemize{
#' \item `"target_right": X = "AACCGT", Y = "A"`
#' \item `"target_middle_lstm": X = (X_1 = "AAC", X_2 = "ATG"), Y = "C"` (note reversed order of X_2)
#' \item `"target_middle_cnn": X = "AACGTA", Y = "C"` 
#' \item `"wavenet": X = "AACCGT", Y = "ACCGTA"`
#' }
#' @param reset_states Whether to reset hidden states of RNN layer at every new input file and before/after validation.
#' @param use_quality_score Whether to use fastq quality scores. If `TRUE` input is not one-hot-encoding but corresponds to probabilities.
#' For example (0.97, 0.01, 0.01, 0.01) instead of (1, 0, 0, 0).
#' @param padding Whether to pad sequences too short for one sample with zeros.
#' @param early_stopping_time Time in seconds after which to stop training.
#' @param validation_only_after_training Whether to skip validation during training and only do one validation iteration after training.
#' @param skip_amb_nuc Threshold of ambiguous nucleotides to accept in fasta entry. Complete entry will get discarded otherwise.
#' @param class_weight List of weights for output. Order should correspond to \code{vocabulary_label}.
#' You can use \code{\link{get_class_weight}} function to estimate class weights:
#' 
#' \code{class_weights <- get_class_weights(path = path, train_type = train_type)}
#' 
#' If \code{train_type = "label_csv"} you need to add path to csv file:
#' 
#' \code{class_weights <- get_class_weights(path = path, train_type = train_type, csv_path = target_from_csv)}
#' @param print_scores Whether to print train/validation scores during training.
#' @param train_val_split_csv A csv file specifying train/validation split. csv file should contain one column named `"file"` and one column named
#' `"type"`. The `"file"` column contains names of fasta/fastq files and `"type"` column specifies if file is used for training or validation.
#' Entries in `"type"` must be named `"train"` or `"val"`, otherwise file will not be used for either. `path` and `path_val` arguments should be the same.
#' Not implemented for `train_type = "label_folder"`.
#' @param set_learning When you want to assign one label to set of samples. Only implemented for `train_type = "label_folder"`.
#' Input is a list with the following parameters 
#' \itemize{
#' \item `samples_per_target`: how many samples to use for one target.
#' \item `maxlen`: length of one sample.
#' \item `reshape_mode`: `"time_dist", "multi_input"` or `"concat"`. 
#' \itemize{
#' \item
#'  If `reshape_mode` is `"multi_input"`, generator will produce `samples_per_target` separate inputs, each of length `maxlen` (model should have
#' `samples_per_target` input layers).
#' \item If reshape_mode is `"time_dist"`, generator will produce a 4D input array. The dimensions correspond to
#' `(batch_size, samples_per_target, maxlen, length(vocabulary))`.
#' \item If `reshape_mode` is `"concat"`, generator will concatenate `samples_per_target` sequences
#' of length `maxlen` to one long sequence.
#' }
#' \item If `reshape_mode` is `"concat"`, there is an additional `buffer_len`
#' argument. If `buffer_len` is an integer, the subsequences are interspaced with `buffer_len` rows. The input length is
#' (`maxlen` \eqn{*} `samples_per_target`) + `buffer_len` \eqn{*} (`samples_per_target` - 1).
#' }
#' @param random_sampling Whether samples should be taken from random positions when using `max_samples` argument. If `FALSE` random 
#' samples are taken from a consecutive subsequence.
#' @param n_gram_stride Step size for n-gram encoding. For AACCGGTT with `n_gram = 4` and `n_gram_stride = 2`, generator encodes
#' `(AACC), (CCGG), (GGTT)`; for `n_gram_stride = 4` generator encodes `(AACC), (GGTT)`.
#' @param callback_list Add additional callbacks to `keras::fit` call.  
#' @param model_card List of arguments for training parameters of training run. Must contain at least an entry `path_model_card`, i.e. the 
#' directory where parameters are stored. List can contain additional (optional) arguments, for example 
#' `model_card = list(path_model_card = "/path/to/logs", description = "transfer learning with BERT model on virus data", ...)`  
#' @examples
#' # create dummy data
#' path_train_1 <- tempfile()
#' path_train_2 <- tempfile()
#' path_val_1 <- tempfile()
#' path_val_2 <- tempfile()
#' 
#' for (current_path in c(path_train_1, path_train_2,
#'                        path_val_1, path_val_2)) {
#'   dir.create(current_path)
#'   create_dummy_data(file_path = current_path,
#'                     num_files = 3,
#'                     seq_length = 10,
#'                     num_seq = 5,
#'                     vocabulary = c("a", "c", "g", "t"))
#' }
#' 
#' # create model
#' model <- create_model_lstm_cnn(layer_lstm = 8, layer_dense = 2, maxlen = 5)
#' 
#' # train model
#' hist <- train_model(train_type = "label_folder",
#'                     model = model,
#'                     path = c(path_train_1, path_train_2),
#'                     path_val = c(path_val_1, path_val_2),
#'                     batch_size = 8,
#'                     epochs = 3,
#'                     steps_per_epoch = 6,
#'                     step = 5,
#'                     format = "fasta",
#'                     vocabulary_label = c("label_1", "label_2"))
#'  
#' @export
train_model <- function(model = NULL,
                        dataset = NULL,
                        dataset_val = NULL,
                        # training args
                        train_val_ratio = 0.2,
                        run_name = "run_1",
                        initial_epoch = 0,
                        class_weight = NULL,
                        print_scores = TRUE,
                        epochs = 10,
                        max_queue_size = 100,
                        steps_per_epoch = 1000,
                        # callbacks
                        path_checkpoint = NULL,
                        path_tensorboard = NULL,
                        path_log = NULL,
                        save_best_only = TRUE,
                        save_weights_only = FALSE,
                        tb_images = FALSE,
                        path_file_log = NULL,
                        reset_states = FALSE,
                        early_stopping_time = NULL,
                        validation_only_after_training = FALSE,
                        train_val_split_csv = NULL,
                        reduce_lr_on_plateau = TRUE,
                        lr_plateau_factor = 0.9,
                        patience = 20,
                        cooldown = 1,
                        model_card = NULL,
                        callback_list = NULL,
                        # generator args
                        train_type = "label_folder",
                        path = NULL,
                        path_val = NULL,
                        batch_size = 64,
                        step = NULL,
                        shuffle_file_order = TRUE,
                        vocabulary = c("a", "c", "g", "t"),
                        format = "fasta",
                        ambiguous_nuc = "zero",
                        seed = c(1234, 4321),
                        file_limit = NULL,
                        use_coverage = NULL,
                        set_learning = NULL,
                        proportion_entries = NULL,
                        sample_by_file_size = FALSE,
                        n_gram = NULL,
                        n_gram_stride = 1,
                        masked_lm = NULL,
                        random_sampling = FALSE,
                        add_noise = NULL,
                        return_int = FALSE,
                        maxlen = NULL,
                        reverse_complement = FALSE,
                        reverse_complement_encoding = FALSE,
                        output_format = "target_right",
                        proportion_per_seq = NULL,
                        read_data = FALSE,
                        use_quality_score = FALSE,
                        padding = FALSE,
                        concat_seq = NULL,
                        target_len = 1,
                        skip_amb_nuc = NULL,
                        max_samples = NULL,
                        added_label_path = NULL,
                        add_input_as_seq = NULL,
                        target_from_csv = NULL,
                        target_split = NULL,
                        shuffle_input = TRUE,
                        vocabulary_label = NULL,
                        delete_used_files = FALSE) {
  
  # initialize metrics, temporary fix
  model <- manage_metrics(model)
  
  run_name <- get_run_name(run_name, path_tensorboard, path_checkpoint, path_log,
                           path_model_card = model_card$path_model_card,
                           auto_extend = TRUE)
  train_with_gen <- is.null(dataset)
  output <- list(tensorboard = FALSE, checkpoints = FALSE)
  if (!is.null(path_tensorboard)) output$tensorboard <- TRUE
  if (!is.null(path_checkpoint)) output$checkpoints <- TRUE
  wavenet_format <- FALSE ; target_middle <- FALSE ; cnn_format <- FALSE
  if (train_type != "label_csv") target_from_csv <- NULL
  
  if (train_with_gen) {
    stopifnot(train_type %in% c("lm", "label_header", "label_folder", "label_csv", "label_rds", "lm_rds", "dummy_gen", "masked_lm"))
    stopifnot(ambiguous_nuc %in% c("zero", "equal", "discard", "empirical"))
    stopifnot(length(vocabulary) == length(unique(vocabulary)))
    stopifnot(length(vocabulary_label) == length(unique(vocabulary_label)))
    labelByFolder <- FALSE
    labelGen <- ifelse(train_type == "lm", FALSE, TRUE)
    
    if (train_type == "label_header") target_from_csv <- NULL
    if (train_type == "label_csv") {
      #train_type <- "label_header"
      if (is.null(target_from_csv)) {
        stop('You need to add a path to csv file for target_from_csv when using train_type = "label_csv"')
      }
      if (!is.null(vocabulary_label)) {
        message("Reading vocabulary_label from csv header")
        if (!is.data.frame(target_from_csv)) {
          output_label_csv <- read.csv2(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
          if (dim(output_label_csv)[2] == 1) {
            output_label_csv <- read.csv(target_from_csv, header = TRUE, stringsAsFactors = FALSE)
          }
        } else {
          output_label_csv <- target_from_csv
        }
        vocabulary_label <- names(output_label_csv)
        vocabulary_label <- vocabulary_label[vocabulary_label != "file"]
      }
    }
    
    if (!is.null(skip_amb_nuc)) {
      if((skip_amb_nuc > 1) | (skip_amb_nuc <0)) {
        stop("skip_amb_nuc should be between 0 and 1 or NULL")
      }
    }
    
    if (!is.null(proportion_per_seq)) {
      if(any(proportion_per_seq > 1) | any(proportion_per_seq  < 0)) {
        stop("proportion_per_seq should be between 0 and 1 or NULL")
      }
    }
    
    # TODO: adjust for multi output model
    # if (!is.null(class_weight) && (length(class_weight) != length(vocabulary_label))) {
    #   stop("class_weight and vocabulary_label must have same length")
    # }
    
    if (!is.null(concat_seq)) {
      if (!is.null(use_coverage)) stop("Coverage encoding not implemented for concat_seq")
    }
    
    # train train_val_ratio via csv file
    if (!is.null(train_val_split_csv)) {
      
      train_val_file <- read.csv2(train_val_split_csv, header = TRUE, stringsAsFactors = FALSE)
      
      if (is.null(path)) {
        path <- train_val_file %>% dplyr::filter(type %in% c("train", "val", "validation")) %>% 
          dplyr::select(file) %>% as.list()
      }
      
      if (train_type == "label_folder") {
        stop('train_val_split_csv not implemented for train_type = "label_folder"')
      }
      if (is.null(path_val)) {
        path_val <- path
      } else {
        if (!all(unlist(path_val) %in% unlist(path))) {
          warning("Train/validation split done via file in train_val_split_csv. Only using files from path argument.")
        }
        path_val <- path
      }
      
      if (dim(train_val_file)[2] == 1) {
        train_val_file <- read.csv(train_val_split_csv, header = TRUE, stringsAsFactors = FALSE)
      }
      train_val_file <- dplyr::distinct(train_val_file)
      
      if (!all(c("file", "type") %in% names(train_val_file))) {
        stop("Column names of train_val_split_csv file must be 'file' and 'type'")
      }
      
      if (length(train_val_file$file) != length(unique(train_val_file$file))) {
        stop("In train_val_split_csv all entires in 'file' column must be unique")
      }
      
      train_files <- train_val_file %>% dplyr::filter(type == "train")
      train_files <- as.character(train_files$file)
      val_files <- train_val_file %>% dplyr::filter(type == "val" | type == "validation")
      val_files <- as.character(val_files$file)
    } else {
      train_files <- NULL
      val_files <- NULL
    }
    
    if (train_type == "lm") {
      stopifnot(output_format %in% c("target_right", "target_middle_lstm", "target_middle_cnn", "wavenet"))
      if (output_format == "target_middle_lstm") target_middle <- TRUE
      if (output_format == "target_middle_cnn") cnn_format <- TRUE
      if (output_format == "wavenet") wavenet_format <- TRUE
    }
    
    if (train_type == "label_header" & is.null(target_from_csv)) {
      stopifnot(!is.null(vocabulary_label))
    }
    
    if (train_type == "label_folder") {
      labelByFolder <- TRUE
      stopifnot(!is.null(vocabulary_label))
      stopifnot(length(path) == length(vocabulary_label))
    }
    
  }
  
  model_weights <- model$get_weights()
  
  # function arguments
  argumentList <- as.list(match.call(expand.dots=FALSE))
  #argumentList <- c(as.list(environment()), list(...)) log default args too
  argumentList <- argumentList[names(argumentList) != ""]
  
  # extract maxlen from model
  if (is.null(maxlen)) {
    maxlen <- get_maxlen(model, set_learning, target_middle, read_data)
  }
  
  if (is.null(step)) step <- maxlen
  vocabulary_label_size <- length(vocabulary_label)
  vocabulary_size <- length(vocabulary)
  
  if (is.null(dataset) && labelByFolder) {
    if (length(path) == 1) warning("Training with just one label")
  }
  
  # add empty hparam dict if non exists
  if (!reticulate::py_has_attr(model, "hparam")) {
    model$hparam <- reticulate::dict()
  }
  
  # tempory file to log training data
  removeLog <- FALSE
  if (is.null(path_file_log)) {
    removeLog <- TRUE
    path_file_log <- tempfile(pattern = "", fileext = ".csv")
  } else {
    if (!endsWith(path_file_log, ".csv")) path_file_log <- paste0(path_file_log, ".csv")
    #path_file_logVal <- tempfile(pattern = "", fileext = ".csv")
  }
  if (reset_states) {
    path_file_logVal <- tempfile(pattern = "", fileext = ".csv")
  } else {
    path_file_logVal <- NULL
  }
  
  # if no dataset is supplied, external fasta generator will generate batches
  if (train_with_gen) {
    #message("Starting fasta generator...")
    
    gen <- get_generator(path = path, batch_size = batch_size, model = model,
                         maxlen = maxlen, step = step, shuffle_file_order = shuffle_file_order,
                         vocabulary = vocabulary, seed = seed[1], proportion_entries = proportion_entries,
                         shuffle_input = shuffle_input, format = format, 
                         path_file_log = path_file_log, reverse_complement = reverse_complement, n_gram_stride = n_gram_stride,
                         output_format = output_format, ambiguous_nuc = ambiguous_nuc,
                         proportion_per_seq = proportion_per_seq, skip_amb_nuc = skip_amb_nuc,
                         use_quality_score = use_quality_score, padding = padding, n_gram = n_gram,
                         added_label_path = added_label_path, add_input_as_seq = add_input_as_seq,
                         max_samples = max_samples, concat_seq = concat_seq, target_len = target_len,
                         file_filter = train_files, use_coverage = use_coverage, random_sampling = random_sampling,
                         train_type = train_type, set_learning = set_learning, file_limit = file_limit,
                         reverse_complement_encoding = reverse_complement_encoding, read_data = read_data,
                         sample_by_file_size = sample_by_file_size, add_noise = add_noise, target_split = target_split,
                         target_from_csv = target_from_csv, masked_lm = masked_lm, return_int = return_int,
                         path_file_logVal = path_file_logVal, delete_used_files = delete_used_files,
                         vocabulary_label = vocabulary_label, new_batch_size = new_batch_size, val = FALSE)
    
    if (!is.null(path_val)) {
      
      gen.val <- get_generator(path = path_val, batch_size = batch_size, model = model,
                               maxlen = maxlen, step = step, shuffle_file_order = shuffle_file_order,
                               vocabulary = vocabulary, seed = seed[2], proportion_entries = proportion_entries,
                               shuffle_input = shuffle_input, format = format, delete_used_files = FALSE,
                               path_file_log = path_file_logVal, reverse_complement = reverse_complement, n_gram_stride = n_gram_stride,
                               output_format = output_format, ambiguous_nuc = ambiguous_nuc,
                               proportion_per_seq = proportion_per_seq, skip_amb_nuc = skip_amb_nuc,
                               use_quality_score = use_quality_score, padding = padding, n_gram = n_gram,
                               added_label_path = added_label_path, add_input_as_seq = add_input_as_seq,
                               max_samples = max_samples, concat_seq = concat_seq, target_len = target_len,
                               file_filter = val_files, use_coverage = use_coverage, random_sampling = random_sampling,
                               train_type = train_type, set_learning = set_learning, file_limit = file_limit,
                               reverse_complement_encoding = reverse_complement_encoding, read_data = read_data,
                               sample_by_file_size = sample_by_file_size, add_noise = add_noise, target_split = target_split,
                               target_from_csv = target_from_csv, masked_lm = masked_lm, return_int = return_int,
                               path_file_logVal = path_file_logVal, vocabulary_label = vocabulary_label,
                               new_batch_size = new_batch_size, val = TRUE)
    } else {
      gen.val <- NULL
    }
    
  }
  
  # skip validation callback
  if (validation_only_after_training | is.null(train_val_ratio) || train_val_ratio == 0) {
    validation_data <- NULL
  } else {
    if (train_with_gen) {
      if (is.null(path_val)) {
        validation_data <- NULL
      } else {
        validation_data <- gen.val
      } 
    }
  }
  
  if (is.null(validation_data)) {
    validation_steps <- NULL
  } else {
    validation_steps <- ceiling(steps_per_epoch * train_val_ratio)
  }
  
  callbacks <- get_callbacks(default_arguments = default_arguments, model = model, path_tensorboard = path_tensorboard, run_name = run_name, train_type = train_type,
                             path_model = path_model, path = path, train_val_ratio = train_val_ratio, batch_size = batch_size, epochs = epochs,
                             max_queue_size = max_queue_size, lr_plateau_factor = lr_plateau_factor, patience = patience, cooldown = cooldown, format = format,
                             steps_per_epoch = steps_per_epoch, step = step, shuffle_file_order = shuffle_file_order, initial_epoch = initial_epoch, vocabulary = vocabulary,
                             learning_rate = learning_rate, shuffle_input = shuffle_input, vocabulary_label = vocabulary_label, solver = solver,
                             file_limit = file_limit, reverse_complement = reverse_complement, wavenet_format = wavenet_format,  cnn_format = cnn_format,
                             train_val_split_csv = train_val_split_csv, n_gram = n_gram, path_file_logVal = path_file_logVal,
                             create_model_function = NULL, vocabulary_size = vocabulary_size, gen_cb = gen_cb, argumentList = argumentList, output = output,
                             maxlen = maxlen, labelGen = labelGen, labelByFolder = labelByFolder, vocabulary_label_size = vocabulary_label_size, tb_images = tb_images,
                             target_middle = target_middle, path_file_log = path_file_log, proportion_per_seq = proportion_per_seq,
                             skip_amb_nuc = skip_amb_nuc, max_samples = max_samples, proportion_entries = proportion_entries, path_log = path_log,
                             train_with_gen = train_with_gen, random_sampling = random_sampling, reduce_lr_on_plateau = reduce_lr_on_plateau,
                             save_weights_only = save_weights_only, path_checkpoint = path_checkpoint, save_best_only = save_best_only, gen.val = gen.val,
                             target_from_csv = target_from_csv, reset_states = reset_states, early_stopping_time = early_stopping_time,
                             validation_only_after_training = validation_only_after_training, model_card = model_card)
  
  # training
  if (train_with_gen) {
    
    model <- keras::set_weights(model, model_weights)
    history <-
      model %>% keras::fit(
        x = gen,
        validation_data = validation_data,
        validation_steps = validation_steps,
        steps_per_epoch = steps_per_epoch,
        max_queue_size = max_queue_size,
        epochs = epochs,
        initial_epoch = initial_epoch,
        callbacks = c(callbacks, callback_list),
        class_weight = class_weight,
        verbose = print_scores)
    
    if (validation_only_after_training) {
      history$val_loss <- model$val_loss
      history$val_acc <- model$val_acc
      model$val_loss <- NULL
      model$val_acc <- NULL
    }
    
  } else {
    
    model <- keras::set_weights(model, model_weights)
    if (!is.null(dataset_val)) {
      validation_data <- list(dataset_val[[1]], dataset_val[[2]])
    } else {
      validation_data <- NULL
    }
    
    history <- keras::fit(
      object = model,
      x = dataset[[1]],
      y = dataset[[2]],
      batch_size = batch_size,
      validation_split = train_val_ratio,
      validation_data = validation_data,
      callbacks = c(callbacks, callback_list),
      epochs = epochs,
      class_weight = class_weight,
      verbose = print_scores)
  }
  
  if (removeLog & file.exists(path_file_log)) {
    file.remove(path_file_log)
  }
  
  message("Training done.")
  
  return(history)
}

#' Generate run_name if none is given or is already present.
#' 
#' If no run name is given, will use date as run name. If run name is already present will add _2 to name or 
#' _x+1 if name ends with _x and x is integer. 
#'
#' @param auto_extend If run_name is already present, add "_2" to name. If name already ends with "_x" replace x with x+1.
#' @keywords internal
get_run_name <- function(run_name = NULL, path_tensorboard, path_checkpoint, path_log, path_model_card, auto_extend = FALSE) {
  
  if (is.null(run_name)) {
    run_name_new <- as.character(Sys.time()) %>% stringr::str_replace_all(" ", "_")
  }
  
  tb_names <- ""
  cp_names <- ""
  log_names <- ""
  mc_names <- ""
  name_present_tb <- FALSE
  name_present_cp <- FALSE
  name_present_log <- FALSE
  name_present_mc <- FALSE
  
  if (!is.null(path_tensorboard)) {
    tb_names <- list.files(path_tensorboard)
    name_present_tb <- (run_name %in% tb_names) # & any(stringr::str_detect(tb_names, run_name))
  }
  if (!is.null(path_checkpoint)) {
    cp_names <- list.files(path_checkpoint)
    name_present_cp <- (run_name %in% cp_names) # & any(stringr::str_detect(cp_names, run_name))  
  }
  if (!is.null(path_log)) {
    log_names <- list.files(path_log)
    name_present_log <- (run_name %in% log_names) # & any(stringr::str_detect(log_names, run_name)) 
  }
  if (!is.null(path_model_card)) {
    mc_names <- list.files(path_model_card)
    name_present_mc <- (run_name %in% mc_names) # & any(stringr::str_detect(log_names, run_name)) 
  }
  
  name_present <- name_present_tb | name_present_cp | name_present_log | name_present_mc
  
  if (name_present & auto_extend) {
    
    ends_with_int <- stringr::str_detect(run_name, "_\\d+$")
    if (ends_with_int) {
      int_ending <- stringr::str_extract(run_name, "\\d+$") %>% as.integer()
      run_name_new <- paste0(stringr::str_remove(run_name, "\\d+$"), int_ending + 1)
    } else {
      run_name_new <- paste0(run_name, "_2")
    }
    
    int_ending <- stringr::str_subset(c(tb_names, cp_names, log_names, mc_names),
                                      paste0("^", stringr::str_remove(run_name, "_\\d+$"))) %>% unique()
    int_ending <- stringr::str_subset(int_ending, "_\\d+$")
    if (length(int_ending) > 0) {
      max_int_ending <- stringr::str_extract(int_ending, "_\\d+$") %>% stringr::str_remove("_") %>% as.integer() %>% max()
      if (!ends_with_int) {
        run_name_new <- paste0(run_name, "_", max_int_ending + 1)
      } else {
        run_name_new <- paste0(stringr::str_remove(run_name, "\\d+$"), max_int_ending + 1)
      }
    }
    
    if (length(int_ending) > 0) {
      name_order <- stringr::str_extract(int_ending, "\\d+$") %>% as.integer() %>% order()
      prev_names <- unique(c(run_name, int_ending[name_order]))
      if (ends_with_int) {
        name_order <- stringr::str_extract(prev_names, "\\d+$") %>% as.integer() %>% order()
        prev_names <- prev_names[name_order]
      }
      
      if (length(prev_names) > 8) {
        old_names_start <- paste(prev_names[1:2], collapse = ", ")
        old_names_end <- paste(prev_names[(length(prev_names)-1) : length(prev_names)], collapse = ", ")
        #old_names <- paste(old_names_start, ",...,", old_names_end) # outputs range of previously used names
        old_names <- run_name
      } else {
        old_names <- paste(prev_names, collapse = ", ")
      }
      message(paste("run_name", old_names, "already present, setting run_name to", run_name_new))
    } else {
      message(paste("run_name", run_name, "already present, setting run_name to", run_name_new))
    }
  }
  
  if (name_present & !auto_extend) {
    stop("run_name already present, please give your run a unique name")
  }
  
  if (!name_present) {
    return(run_name)
  }
  
  return(run_name_new)
}

#' @title Train CPC inspired model
#'   
#' @description
#' Train a CPC (Oord et al.) inspired neural network on genomic data.
#' 
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_folder
#' @inheritParams generator_fasta_label_header_csv
#' @param train_type Either `"cpc"`, `"Self-GenomeNet"`. 
#' @param encoder A keras encoder for the cpc function. 
#' @param context A keras context model for the cpc function.
#' @param path Path to training data. If \code{train_type} is \code{label_folder}, should be a vector or list
#' where each entry corresponds to a class (list elements can be directories and/or individual files). If \code{train_type} is not \code{label_folder}, 
#' can be a single directory or file or a list of directories and/or files.
#' @param path_val Path to validation data. See `path` argument for details.
#' @param path_checkpoint Path to checkpoints folder or `NULL`. If `NULL`, checkpoints don't get stored.
#' @param path_tensorboard Path to tensorboard directory or `NULL`. If `NULL`, training not tracked on tensorboard.
#' @param train_val_ratio For generator defines the fraction of batches that will be used for validation (compared to size of training data), i.e. one validation iteration
#' processes \code{batch_size} \eqn{*} \code{steps_per_epoch} \eqn{*} \code{train_val_ratio} samples. If you use dataset instead of generator and \code{dataset_val} is `NULL`, splits \code{dataset}
#' into train/validation data.
#' @param run_name Name of the run. Name will be used to identify output from callbacks.
#' @param batch_size Number of samples used for one network update.
#' @param epochs Number of iterations.
#' @param steps_per_epoch Number of training batches per epoch.
#' @param shuffle_file_order Boolean, whether to go through files sequentially or shuffle beforehand.
#' @param initial_epoch Epoch at which to start training. Note that network
#' will run for (\code{epochs} - \code{initial_epochs}) rounds and not \code{epochs} rounds.
#' @param seed Sets seed for reproducible results.
#' @param file_limit Integer or `NULL`. If integer, use only specified number of randomly sampled files for training. Ignored if greater than number of files in \code{path}.
#' @param patchlen The length of a patch when splitting the input sequence.
#' @param nopatches The number of patches when splitting the input sequence. 
#' @param step Frequency of sampling steps.
#' @param preloadGeneratorpath 
#' @param stride The overlap between two patches when splitting the input sequence.
#' @param pretrained_model A pretrained keras model, for which training will be continued
#' @param learningrate A Tensor, floating point value. If a schedule is defines, this value gives the initial learning rate. Defaults to 0.001.
#' @param learningrate_schedule A schedule for a non-constant learning rate over the training. Either "cosine_annealing", "step_decay", or "exp_decay".
#' @param k Value of k for sparse top k categorical accuracy. Defaults to 5.
#' @param stepsmin In CPC, a patch is predicted given another patch. stepsmin defines how many patches between these two should be ignored during prediction.
#' @param stepsmax The maximum distance between the predicted patch and the given patch.
#' @param emb_scale Scales the impact of a patches context.
#' @examples
#' 
#' #create dummy data
#' path_train_1 <- tempfile()
#' path_train_2 <- tempfile()
#' path_val_1 <- tempfile()
#' path_val_2 <- tempfile()
#' 
#' for (current_path in c(path_train_1, path_train_2,
#'                        path_val_1, path_val_2)) {
#'   dir.create(current_path)
#'   deepG::create_dummy_data(file_path = current_path,
#'                            num_files = 3,
#'                            seq_length = 10,
#'                            num_seq = 5,
#'                            vocabulary = c("a", "c", "g", "t"))
#' }
#' 
#' # create model
#' encoder <- function(maxlen = NULL,
#'                     patchlen = NULL,
#'                     nopatches = NULL,
#'                     eval = FALSE) {
#'   if (is.null(nopatches)) {
#'     nopatches <- nopatchescalc(patchlen, maxlen, patchlen * 0.4)
#'   }
#'   inp <- keras::layer_input(shape = c(maxlen, 4))
#'   stridelen <- as.integer(0.4 * patchlen)
#'   createpatches <- inp %>%
#'     keras::layer_reshape(list(maxlen, 4L, 1L), name = "prep_reshape1", dtype = "float32") %>%
#'     tensorflow::tf$image$extract_patches(
#'       sizes = list(1L, patchlen, 4L, 1L),
#'       strides = list(1L, stridelen, 4L, 1L),
#'       rates = list(1L, 1L, 1L, 1L),
#'       padding = "VALID",
#'       name = "prep_patches"
#'     ) %>%
#'     keras::layer_reshape(list(nopatches, patchlen, 4L),
#'                          name = "prep_reshape2") %>%
#'     tensorflow::tf$reshape(list(-1L, patchlen, 4L),
#'                            name = "prep_reshape3")
#' 
#'   danQ <- createpatches %>%
#'     keras::layer_conv_1d(
#'       input_shape = c(maxlen, 4L),
#'       filters = 320L,
#'       kernel_size = 26L,
#'       activation = "relu"
#'     ) %>%
#'     keras::layer_max_pooling_1d(pool_size = 13L, strides = 13L) %>%
#'     keras::layer_dropout(0.2) %>%
#'     keras::layer_lstm(units = 320, return_sequences = TRUE) %>%
#'     keras::layer_dropout(0.5) %>%
#'     keras::layer_flatten() %>%
#'     keras::layer_dense(925, activation = "relu")
#'   patchesback <- danQ %>%
#'     tensorflow::tf$reshape(list(-1L, tensorflow::tf$cast(nopatches, tensorflow::tf$int16), 925L))
#'   keras::keras_model(inp, patchesback)
#' }
#' 
#' context <- function(latents) {
#'   cres <- latents
#'   cres_dim = cres$shape
#'   predictions <-
#'     cres %>%
#'     keras::layer_lstm(
#'       return_sequences = TRUE,
#'       units = 256,  # WAS: 2048,
#'       name = paste("context_LSTM_1",
#'                    sep = ""),
#'       activation = "relu"
#'     )
#'   return(predictions)
#' }
#' 
#' # train model
#' temp_dir <- tempdir()
#' dir.create(temp_dir)
#' hist <- train_model_cpc(train_type = "CPC",
#'                         ### cpc functions ###
#'                         encoder = encoder,
#'                         context = context,
#'                         #### Generator settings ####
#'                         path_checkpoint = temp_dir,
#'                         path = c(path_train_1, path_train_2),
#'                         path_val = c(path_val_1, path_val_2),
#'                         run_name = "TEST",
#'                         batch_size = 8,
#'                         epochs = 3,
#'                         steps_per_epoch = 6,
#'                         patchlen = 100,
#'                         nopatches = 8)
#'  
#' @export
train_model_cpc <-
  function(arch_type = "CPC",
           ### cpc functions ###
           encoder = NULL,
           context = NULL,
           #### Generator settings ####
           path,
           path_val = NULL,
           format = "fasta",
           delete_used_files = FALSE,
           path_checkpoint = NULL,
           path_tensorboard = NULL,
           train_val_ratio = 0.2,
           run_name,
           
           batch_size = 32,
           epochs = 100,
           steps_per_epoch = 2000,
           shuffle_file_order = F,
           initial_epoch = 1,
           seed = 1234,
           
           path_file_log = T,
           train_val_split_csv = NULL,
           file_limit = NULL,
           proportion_per_seq = NULL,
           max_samples = NULL,
           maxlen = NULL,
           
           patchlen = NULL,
           nopatches = NULL,
           step = NULL,
           preloadGeneratorpath = NULL,
           file_filter = NULL,
           stride = 0.4,
           pretrained_model = NULL,
           learningrate = 0.001,
           learningrate_schedule = NULL,
           k = 5,
           stepsmin = 2,
           stepsmax = 3,
           emb_scale = 0.1) {
    
    # Stride is default 0.4 x patchlen FOR NOW
    stride <- 0.4
    
    patchlen <- as.integer(patchlen)
    
    ########################################################################################################
    ############################### Warning messages if wrong initialization ###############################
    ########################################################################################################
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Model specification ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ## Three options:
    ## 1. Define Maxlen and Patchlen
    ## 2. Define Number of patches and Patchlen
    ## ---> in both cases the respectively missing value will be calculated
    ## 3. Pretrained model is giving specs
    ## error if none of those is fulfilled
    
    if (is.null(pretrained_model)) {
      ## If no pretrained model, patchlen has to be defined
      if (is.null(patchlen)) {
        stop("Please define patchlen")
      }
      ## Either maxlen or number of patches is needed
      if (is.null(maxlen) & is.null(nopatches)) {
        stop("Please define either maxlen or nopatches")
        ## the respectively missing value will be calculated
      } else if (is.null(maxlen) & !is.null(nopatches)) {
        maxlen <- (nopatches - 1) * (stride * patchlen) + patchlen
      } else if (!is.null(maxlen) & is.null(nopatches)) {
        nopatches <-
          as.integer((maxlen - patchlen) / (stride * patchlen) + 1)
      }
      ## if step is not defined, we do not use overlapping sequences
      if (is.null(step)) {
        step = maxlen
      }
    } else if (!is.null(pretrained_model)) {
      specs <-
        readRDS(paste(
          sub("/[^/]+$", "", pretrained_model),
          "modelspecs.rds",
          sep = "/"
        ))
      patchlen          <- specs$patchlen
      maxlen            <- specs$maxlen
      nopatches         <- specs$nopatches
      stride            <- specs$stride
      step              <- specs$step
      k                 <- specs$k
      emb_scale         <- specs$emb_scale
    }
    
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Learning rate schedule ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ## If learning_rate schedule is wanted, all necessary parameters must be given
    LRstop(learningrate_schedule)
    ########################################################################################################
    #################################### Preparation: Data, paths metrics ##################################
    ########################################################################################################
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Path definition ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    # Check if 'a' matches the desired format, already has a datetime
    if (grepl("_(\\d{6}_\\d{6})$", run_name)) {
      runname <- run_name
    } else {
      runname <-
        paste0(run_name , format(Sys.time(), "_%y%m%d_%H%M%S"))
    }
    
    ## Create folder for model
    dir.create(paste(path_checkpoint, runname, sep = "/"))
    dir <- paste(path_checkpoint, runname, sep = "/")
    
    ## Create folder for filelog
    if (!is.na(path_checkpoint)) {
      path_file_log <-
        paste(path_checkpoint, runname, "filelog.csv", sep = "/")
      path_file_logV <-
        paste(path_checkpoint, runname, "filelog_val.csv", sep = "/")
    } else {
      path_file_log <- NULL
    }
    
    GenConfig <-
      GenParams(maxlen, batch_size, step, proportion_per_seq, max_samples)
    GenTConfig <-
      GenTParams(path, shuffle_file_order, path_file_log, seed)
    if (is.null(preloadGeneratorpath)) {
      GenVConfig <- GenVParams(path_val, shuffle_file_order, path_file_logV)
    }
    
    # train train_val_ratio via csv file
    if (!is.null(train_val_split_csv)) {
      if (is.null(path_val)) {
        path_val <- path
      } else {
        if (!all(unlist(path_val) %in% unlist(path))) {
          warning("Train/validation split done via file in train_val_split_csv. Only using files from path argument.")
        }
        path_val <- path
      }
      
      train_val_file <- read.csv2(train_val_split_csv, header = TRUE, stringsAsFactors = FALSE)
      if (dim(train_val_file)[2] == 1) {
        train_val_file <- read.csv(train_val_split_csv, header = TRUE, stringsAsFactors = FALSE)
      }
      train_val_file <- dplyr::distinct(train_val_file)
      
      if (!all(c("file", "type") %in% names(train_val_file))) {
        stop("Column names of train_val_split_csv file must be 'file' and 'type'")
      }
      
      if (length(train_val_file$file) != length(unique(train_val_file$file))) {
        stop("In train_val_split_csv all entires in 'file' column must be unique")
      }
      
      file_filter <- list()
      file_filter[[1]] <- train_val_file %>% dplyr::filter(type == "train")
      file_filter[[1]] <- as.character(file_filter[[1]]$file)
      file_filter[[2]] <- train_val_file %>% dplyr::filter(type == "val" | type == "validation")
      file_filter[[2]] <- as.character(file_filter[[2]]$file)
    }
    
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ File count ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    if (is.null(file_filter) && is.null(train_val_split_csv)) {
      if (is.null(file_limit)) {
        if (is.list(path)) {
          num_files <- 0
          for (i in seq_along(path)) {
            num_files <- num_files + length(list.files(path[[i]]))
          }
        } else {
          num_files <- length(list.files(path))
        }
      } else {
        num_files <- file_limit * length(path)
      }
    } else {
      num_files <- length(file_filter[1])
    }
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Creation of generators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    cat(format(Sys.time(), "%F %R"), ": Preparing the data\n")
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Training Generator ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    fastrain <- get_generator(maxlen = maxlen, batch_size = batch_size, step = step, proportion_per_seq = proportion_per_seq, max_samples = max_samples,
                              path = path, shuffle_file_order = shuffle_file_order, path_file_log = path_file_log, seed = seed,
                              file_filter = file_filter[1], train_type="lm_rds", format=format, delete_used_files = delete_used_files)
      # do.call(get_generator,
      #         c(GenConfig, GenTConfig, ))
    
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Validation Generator ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    if (!is.null(preloadGeneratorpath)) {
      fasval <-
        readPLGpar(preloadGeneratorpath, batch_size, maxlen, 8, seed)
    } else{
      fasval <- get_generator(maxlen = maxlen, batch_size = batch_size, step = step, proportion_per_seq = proportion_per_seq, max_samples = max_samples,
                              path = path_val, shuffle_file_order = shuffle_file_order, path_file_log = path_file_logV, seed = seed,
                              file_filter = file_filter[1], train_type="lm_rds", format=format, delete_used_files = delete_used_files)
        # do.call(
        #   get_generator,
        #   c(
        #     GenConfig,
        #     GenVConfig,
        #     seed = seed,
        #     file_filter = file_filter[2], 
        #     train_type="lm", format=format, 
        #     delete_used_files = delete_used_files
        #   )
        # )
    }
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Creation of metrics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    cat(format(Sys.time(), "%F %R"), ": Preparing the metrics\n")
    train_loss <- tensorflow::tf$keras$metrics$Mean(name = 'train_loss')
    val_loss <- tensorflow::tf$keras$metrics$Mean(name = 'val_loss')
    train_acc <- tensorflow::tf$keras$metrics$Mean(name = 'train_acc')
    val_acc <- tensorflow::tf$keras$metrics$Mean(name = 'val_acc')
    
    ########################################################################################################
    ###################################### History object preparation ######################################
    ########################################################################################################
    
    .GlobalEnv$history <- list(
      params = list(
        batch_size = batch_size,
        epochs = 0,
        steps = steps_per_epoch,
        samples = steps_per_epoch * batch_size,
        verbose = 1,
        do_validation = TRUE,
        metrics = c("loss", "accuracy", "val_loss", "val_accuracy")
      ),
      metrics = list(
        loss = c(),
        accuracy = c(),
        val_loss = c(),
        val_accuracy = c()
      )
    )
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Reformat to S3 object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    class(.GlobalEnv$history) <- "keras_training_history"
    
    ########################################################################################################
    ############################################ Model creation ############################################
    ########################################################################################################

      ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Unsupervised Build from scratch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
      cat(format(Sys.time(), "%F %R"), ": Creating the model\n")
      ## Build encoder
      enc <-
        encoder(maxlen = maxlen,
                patchlen = patchlen,
                nopatches = nopatches)
      
      ## Build model
      model <-
        keras::keras_model(
          enc$input,
          cpcloss(
            enc$output,
            context,
            batch_size = batch_size,
            steps_to_ignore = stepsmin,
            steps_to_predict = stepsmax,
            train_type = arch_type,
            k = k,
            emb_scale = emb_scale
          )
        )
      
      ## Build optimizer
      optimizer <- keras::optimizer_adam(
        learning_rate = learningrate,
        beta_1 = 0.8,
        epsilon = 10 ^ -8,
        decay = 0.999,
        clipnorm = 0.01
      )
      ####~~~~~~~~~~~~~~~~~~~~~~~~~~ Unsupervised Read if pretrained model given ~~~~~~~~~~~~~~~~~~~~~~~~~####
      
      if (!is.null(pretrained_model)) {
      cat(format(Sys.time(), "%F %R"), ": Loading the trained model.\n")
      ## Read model
      #model <- keras::load_model_hdf5(pretrained_model, compile = F)
        # Load the weights into the new model
        model %>% keras::load_model_weights_hdf5(pretrained_model)
        optimizer <- ReadOpt(pretrained_model)
        optimizer$learning_rate$assign(learningrate)
    }
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Saving necessary model objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ## optimizer configuration
    saveRDS(
      optimizer$get_config(),
      paste(dir, "optconfig.rds", sep = "/")
    )
    ## model parameters
    saveRDS(
      list(
        maxlen = maxlen,
        patchlen = patchlen,
        stride = stride,
        nopatches = nopatches,
        step = step,
        batch_size = batch_size,
        epochs = epochs,
        steps_per_epoch = steps_per_epoch,
        train_val_ratio = train_val_ratio,
        max_samples = max_samples,
        k = k,
        emb_scale = emb_scale,
        learningrate = learningrate
      ),
      paste(dir, "modelspecs.rds", sep = "/")
    )
    
    ########################################################################################################
    ######################################## Tensorboard connection ########################################
    ########################################################################################################
    
    if (!is.null(path_tensorboard)) {
      ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Initialize Tensorboard writers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
      logdir <- path_tensorboard
      
      dir.create(file.path(logdir, runname))
      dir.create(file.path(logdir, runname, "train"))
      dir.create(file.path(logdir, runname, "validation"))
      
      writertrain <-
        tensorflow::tf$summary$create_file_writer(file.path(logdir, runname, "train"))
      writerval <-
        tensorflow::tf$summary$create_file_writer(file.path(logdir, runname, "validation"))
      
      ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Write parameters to Tensorboard ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
      tftext <-
        lapply(as.list(match.call())[-1][-c(1, 2)], function(x)
          ifelse(all(nchar(deparse(
            eval(x)
          )) < 20) && !is.null(eval(x)), eval(x), deparse(x)))
      
      with(writertrain$as_default(), {
        tensorflow::tf$summary$text("Specification",
                                    paste(
                                      names(tftext),
                                      tftext,
                                      sep = " = ",
                                      collapse = "  \n"
                                    ),
                                    step = 0L)
      })
    }
    
    ########################################################################################################
    ######################################## Training loop function ########################################
    ########################################################################################################
    
    train_val_loop <-
      function(batches = steps_per_epoch, epoch, train_val_ratio) {
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start of loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
        for (i in c("train", "val")) {
          if (i == "val") {
            ## Calculate steps for validation
            batches <- ceiling(batches * train_val_ratio)
          }
          
          for (b in seq(batches)) {
            if (i == "train") {
              ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Training step ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
              ## If Learning rate schedule specified, calculate learning_rate for current epoch
              if (!is.null(learningrate_schedule)) {
                optimizer$learning_rate$assign(getEpochLR(learningrate_schedule, epoch))
              }
              ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Optimization step ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
              
              with(tensorflow::tf$GradientTape() %as% tape, {
                out <-
                  modelstep(fastrain(),
                            model,
                            arch_type,
                            TRUE)
                l <- out[1]
                acc <- out[2]
              })
              
              gradients <-
                tape$gradient(l, model$trainable_variables)
              optimizer$apply_gradients(purrr::transpose(list(
                gradients, model$trainable_variables
              )))
              train_loss(l)
              train_acc(acc)
              
            } else {
              ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Validation step ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
              if (is.null(preloadGeneratorpath)) {
                out <-
                  modelstep(fasval(),
                            model,
                            arch_type,
                            F)
              } else {
                out <-
                  modelstep(fasval$gen(),
                            model,
                            arch_type,
                            F)
              }
              l <- out[1]
              acc <- out[2]
              val_loss(l)
              val_acc(acc)
              
            }
            
            ## Print status of epoch
            if (b %in% seq(0, batches, by = batches / 10)) {
              cat("-")
              # Specify the path to the folder containing the H5 file
              folder_path <- file.path(path_checkpoint, runname)
              
              # List all files in the folder with the .h5 extension
              h5_files <-
                list.files(folder_path, pattern = "\\.h5$", full.names = TRUE)
              
              # Check if there are any H5 files in the folder
              if (length(h5_files) > 0) {
                # Get the last modification time of the most recent H5 file
                last_modification_time <-
                  file.info(h5_files)$mtime
                
                # Calculate the time difference in hours
                time_difference_hours <-
                  as.numeric(difftime(Sys.time(), last_modification_time, units = "hours"))
                
                # Check if the time difference is greater than 2 hours
                if (min(time_difference_hours) > 2) {
                  savechecks("backup_time",
                             runname,
                             model,
                             optimizer,
                             history,
                             path_checkpoint)
                }
              }
            }
          }
          
          ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End of Epoch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
          if (i == "train") {
            ## Training step
            # Write epoch result metrics value to tensorboard
            if (!is.null(path_tensorboard)) {
              TB_loss_acc(writertrain, train_loss, train_acc, epoch)
              with(writertrain$as_default(), {
                tensorflow::tf$summary$scalar('epoch_lr',
                                              optimizer$learning_rate,
                                              step = tensorflow::tf$cast(epoch, "int64"))
                tensorflow::tf$summary$scalar(
                  'training files seen',
                  nrow(
                    readr::read_csv(
                      path_file_log,
                      col_names = FALSE,
                      col_types = readr::cols()
                    )
                  ) / num_files,
                  step = tensorflow::tf$cast(epoch, "int64")
                )
              })
            }
            # Print epoch result metric values to console
            cat(format(Sys.time(), "%F %R"), ": Training Step Done\n")
            tensorflow::tf$print(" Train Loss",
                                 train_loss$result(),
                                 ", Train Acc",
                                 train_acc$result())
            
            # Save epoch result metric values to history object
            .GlobalEnv$history$params$epochs <- epoch
            .GlobalEnv$history$metrics$loss[epoch] <-
              as.double(train_loss$result())
            .GlobalEnv$history$metrics$accuracy[epoch]  <-
              as.double(train_acc$result())
            
            # Reset states
            train_loss$reset_states()
            train_acc$reset_states()
            
          } else {
            ## Validation step
            # Write epoch result metrics value to tensorboard
            if (!is.null(path_tensorboard)) {
              TB_loss_acc(writerval, val_loss, val_acc, epoch)
            }
            
            # Print epoch result metric values to console
            cat(format(Sys.time(), "%F %R"), ": Validation Step Done\n")
            tensorflow::tf$print(" Validation Loss",
                                 val_loss$result(),
                                 ", Validation Acc",
                                 val_acc$result())
            
            # save results globally for best model saving condition
            if (b == max(seq(batches))) {
              .GlobalEnv$eploss[[epoch]] <- as.double(val_loss$result())
              .GlobalEnv$epacc[[epoch]] <-
                as.double(val_acc$result())
            }
            
            # Save epoch result metric values to history object
            .GlobalEnv$history$metrics$val_loss[epoch] <-
              as.double(val_loss$result())
            .GlobalEnv$history$metrics$val_accuracy[epoch]  <-
              as.double(val_acc$result())
            
            # Reset states
            val_loss$reset_states()
            val_acc$reset_states()
            if (!is.null(preloadGeneratorpath)) {
              fasval$reset()
            }
          }
        }
      }
    
    ########################################################################################################
    ############################################# Training run #############################################
    ########################################################################################################
    
    # initialize global list of validation results for best model saving condition
    .GlobalEnv$eploss <- list()
    .GlobalEnv$epacc <- list()
    
    cat(format(Sys.time(), "%F %R"), ": Starting Training\n")
    
    ## Training loop
    for (i in seq(initial_epoch, (epochs + initial_epoch - 1))) {
      cat(format(Sys.time(), "%F %R"), ": EPOCH", i, " \n")
      
      ## Epoch loop
      train_val_loop(epoch = i, train_val_ratio = train_val_ratio)
      
      ## Save checkpoints
      # best model (smallest loss)
      if (eploss[[i]] == min(unlist(eploss))) {
        savechecks("best", runname, model, optimizer, history, path_checkpoint)
      }
      # backup model every 2 (initally 10) epochs
      if (i %% 2 == 0) {
        savechecks("backup", runname, model, optimizer, history, path_checkpoint)
      }
    }
    
    ########################################################################################################
    ############################################# Final saves ##############################################
    ########################################################################################################
    
    savechecks(cp = "FINAL", runname, model, optimizer, history, path_checkpoint)
    if (!is.null(path_tensorboard)) {
      writegraph <-
        tensorflow::tf$keras$callbacks$TensorBoard(file.path(logdir, runname))
      writegraph$set_model(model)
    }
  }
