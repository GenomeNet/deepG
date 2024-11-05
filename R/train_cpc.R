#' @title Train CPC inspired model
#'   
#' @description
#' Train a CPC (Oord et al.) inspired neural network on genomic data.
#' 
#' @inheritParams generator_fasta_lm
#' @inheritParams generator_fasta_label_folder
#' @inheritParams generator_fasta_label_header_csv
#' @inheritParams train_model
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
#' @param stride The overlap between two patches when splitting the input sequence.
#' @param pretrained_model A pretrained keras model, for which training will be continued
#' @param learningrate A Tensor, floating point value. If a schedule is defines, this value gives the initial learning rate. Defaults to 0.001.
#' @param learningrate_schedule A schedule for a non-constant learning rate over the training. Either "cosine_annealing", "step_decay", or "exp_decay".
#' @param k Value of k for sparse top k categorical accuracy. Defaults to 5.
#' @param stepsmin In CPC, a patch is predicted given another patch. stepsmin defines how many patches between these two should be ignored during prediction.
#' @param stepsmax The maximum distance between the predicted patch and the given patch.
#' @param emb_scale Scales the impact of a patches context.
#' @examplesIf reticulate::py_module_available("tensorflow")
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
#'  
#' @returns A list of training metrics.  
#' @export
train_model_cpc <-
  function(train_type = "CPC",
           ### cpc functions ###
           encoder = NULL,
           context = NULL,
           #### Generator settings ####
           path,
           path_val = NULL,
           path_checkpoint = NULL,
           path_tensorboard = NULL,
           train_val_ratio = 0.2,
           run_name,
           
           batch_size = 32,
           epochs = 100,
           steps_per_epoch = 2000,
           shuffle_file_order = FALSE,
           initial_epoch = 1,
           seed = 1234,
           
           path_file_log = TRUE,
           train_val_split_csv = NULL,
           file_limit = NULL,
           proportion_per_seq = NULL,
           max_samples = NULL,
           maxlen = NULL,
           
           patchlen = NULL,
           nopatches = NULL,
           step = NULL,
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
    runname <-
      paste0(run_name , format(Sys.time(), "_%y%m%d_%H%M%S"))
    
    ## Create folder for model
    if (!is.null(path_checkpoint)) {
      dir.create(paste(path_checkpoint, runname, sep = "/"))
      dir <- paste(path_checkpoint, runname, sep = "/")
      ## Create folder for filelog
      path_file_log <-
        paste(path_checkpoint, runname, "filelog.csv", sep = "/")
    } else {
      path_file_log <- NULL
    }
    
    GenConfig <-
      GenParams(maxlen, batch_size, step, proportion_per_seq, max_samples)
    GenTConfig <-
      GenTParams(path, shuffle_file_order, path_file_log, seed)
    GenVConfig <- GenVParams(path_val, shuffle_file_order)
    
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
      
      train_val_file <- utils::read.csv2(train_val_split_csv, header = TRUE, stringsAsFactors = FALSE)
      if (dim(train_val_file)[2] == 1) {
        train_val_file <- utils::read.csv(train_val_split_csv, header = TRUE, stringsAsFactors = FALSE)
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
    message(format(Sys.time(), "%F %R"), ": Preparing the data\n")
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Training Generator ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    fastrain <-
      do.call(generator_fasta_lm,
              c(GenConfig, GenTConfig, file_filter = file_filter[1]))
    
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Validation Generator ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    fasval <-
      do.call(
        generator_fasta_lm,
        c(
          GenConfig,
          GenVConfig,
          seed = seed,
          file_filter = file_filter[2]
        )
      )
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Creation of metrics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    message(format(Sys.time(), "%F %R"), ": Preparing the metrics\n")
    train_loss <- tensorflow::tf$keras$metrics$Mean(name = 'train_loss')
    val_loss <- tensorflow::tf$keras$metrics$Mean(name = 'val_loss')
    train_acc <- tensorflow::tf$keras$metrics$Mean(name = 'train_acc')
    val_acc <- tensorflow::tf$keras$metrics$Mean(name = 'val_acc')
    
    ########################################################################################################
    ###################################### History object preparation ######################################
    ########################################################################################################
    
    history <- list(
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
    
    eploss <- list()
    epacc <- list()
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Reformat to S3 object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    class(history) <- "keras_training_history"
    
    ########################################################################################################
    ############################################ Model creation ############################################
    ########################################################################################################
    if (is.null(pretrained_model)) {
      ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Unsupervised Build from scratch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
      message(format(Sys.time(), "%F %R"), ": Creating the model\n")
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
            train_type = train_type,
            k = k,
            emb_scale = emb_scale
          )
        )
      
      ## Build optimizer
      optimizer <- # keras::optimizer_adam(
        tensorflow::tf$keras$optimizers$legacy$Adam(
          learning_rate = learningrate,
          beta_1 = 0.8,
          epsilon = 10 ^ -8,
          decay = 0.999,
          clipnorm = 0.01
        )
      ####~~~~~~~~~~~~~~~~~~~~~~~~~~ Unsupervised Read if pretrained model given ~~~~~~~~~~~~~~~~~~~~~~~~~####
      
    } else {
      message(format(Sys.time(), "%F %R"), ": Loading the trained model.\n")
      ## Read model
      model <- keras::load_model_hdf5(pretrained_model, compile = FALSE)
      optimizer <- ReadOpt(pretrained_model)
      optimizer$learning_rate$assign(learningrate)
    }
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Saving necessary model objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ## optimizer configuration
    
    if (!is.null(path_checkpoint)) {
      saveRDS(optimizer$get_config(),
              paste(dir, "optconfig.rds", sep = "/"))
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
    }
    ########################################################################################################
    ######################################## Tensorboard connection ########################################
    ########################################################################################################
    
    if (!is.null(path_tensorboard)) {
      ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Initialize Tensorboard writers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
      logdir <- path_tensorboard
      writertrain <-
        tensorflow::tf$summary$create_file_writer(file.path(logdir, runname, "/train"))
      writerval <-
        tensorflow::tf$summary$create_file_writer(file.path(logdir, runname, "/validation"))
      
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
              
              #with(tensorflow::tf$GradientTape() %as% tape, {
              with(reticulate::`%as%`(tensorflow::tf$GradientTape(), tape), {
                
                out <-
                  modelstep(fastrain(),
                            model,
                            train_type,
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
              out <-
                modelstep(fasval(),
                          model,
                          train_type,
                          FALSE)
              
              l <- out[1]
              acc <- out[2]
              val_loss(l)
              val_acc(acc)
              
            }
            
            ## Print status of epoch
            if (b %in% seq(0, batches, by = batches / 10)) {
              message("-")
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
            tensorflow::tf$print(" Train Loss",
                                 train_loss$result(),
                                 ", Train Acc",
                                 train_acc$result())
            
            # Save epoch result metric values to history object
            history$params$epochs <- epoch
            history$metrics$loss[epoch] <-
              as.double(train_loss$result())
            history$metrics$accuracy[epoch]  <-
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
            tensorflow::tf$print(" Validation Loss",
                                 val_loss$result(),
                                 ", Validation Acc",
                                 val_acc$result())
            
            # save results for best model saving condition
            if (b == max(seq(batches))) {
              eploss[[epoch]] <- as.double(val_loss$result())
              epacc[[epoch]] <-
                as.double(val_acc$result())
            }
            
            # Save epoch result metric values to history object
            history$metrics$val_loss[epoch] <-
              as.double(val_loss$result())
            history$metrics$val_accuracy[epoch]  <-
              as.double(val_acc$result())
            
            # Reset states
            val_loss$reset_states()
            val_acc$reset_states()
          }
        }
        return(list(history,eploss,epacc))
      }
    
    ########################################################################################################
    ############################################# Training run #############################################
    ########################################################################################################
    
    
    message(format(Sys.time(), "%F %R"), ": Starting Training\n")
    
    ## Training loop
    for (i in seq(initial_epoch, (epochs + initial_epoch - 1))) {
      message(format(Sys.time(), "%F %R"), ": EPOCH ", i, " \n")
      
      ## Epoch loop
      out <- train_val_loop(epoch = i, train_val_ratio = train_val_ratio)
      history <- out[[1]]
      eploss <- out[[2]]
      epacc <- out[[3]]
      ## Save checkpoints
      # best model (smallest loss)
      if (eploss[[i]] == min(unlist(eploss))) {
        savechecks("best", runname, model, optimizer, history, path_checkpoint)
      }
      # backup model every 10 epochs
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
