# library(magrittr)
# input_tensor <- array(runif(1*5*4), dim = c(1,5,4))
# input_tensor[1,,]
# target_class_idx <- 1
# model <- keras::load_model_hdf5("/home/rmreches/checkpoints/language_model_1_checkpoints/Ep.002-val_loss0.33-val_acc0.750.hdf5",
#                                 compile = FALSE)
# model <- keras::load_model_weights_hdf5(model,
#                                         "/home/rmreches/checkpoints/language_model_1_checkpoints/Ep.002-val_loss0.33-val_acc0.750.hdf5")
# input_idx = 1
# pred_stepwise = FALSE
# 
# m_steps = 50
# baseline_type = "zero"
# input_seq <- input_tensor
#' 
#' 
#' library(reticulate)
#' py_run_string("import tensorflow as tf")
#' input_idx <- NULL
#' 
#' baseline_seq <- interpolate_seq(m_steps = m_steps,
#'                                 baseline_type = baseline_type,
#'                                 input_seq = input_seq)
#' 
#' baseline_seq <- tensorflow::tf$cast(baseline_seq, dtype = "float32")
#' 
#' if (is.list(input_seq)) {
#'   path_gradients <- list()
#'   avg_grads <- list()
#'   ig <- list()
#'   for (i in 1:length(input_seq)) {
#'     path_gradients[[i]] <- compute_gradients(
#'       model = model,
#'       input_tensor = baseline_seq,
#'       target_class_idx = target_class_idx,
#'       input_idx = i,
#'       pred_stepwise = pred_stepwise)
#'     
#'     avg_grads[[i]] <- integral_approximation(gradients = path_gradients[[i]])
#'     ig[[i]] <- ((input_seq[[i]] - baseline_seq[[i]][1, , ]) * avg_grads[[i]])[1, , ]
#'   }
#' } else {
#'   path_gradients <- compute_gradients(
#'     model = model,
#'     input_tensor = baseline_seq,
#'     target_class_idx = target_class_idx,
#'     input_idx = NULL,
#'     pred_stepwise = pred_stepwise)
#'   
#'   avg_grads <- integral_approximation(gradients = path_gradients)
#'   ig <- ((input_seq - baseline_seq[1, , ]) * avg_grads)[1, , ]
#' }
#' 
#' path_gradients_comp <- compute_gradients(
#'   model = model,
#'   input_tensor = baseline_seq,
#'   target_class_idx = target_class_idx,
#'   input_idx = NULL,
#'   pred_stepwise = pred_stepwise)
#' dim(path_gradients_comp)
#' 
#' in_tensor <- as.array(baseline_seq[1, , ])
#' in_tensor <- array(in_tensor, dim = c(1, dim(baseline)[-1]))
#' dim(in_tensor)
#' in_tensor <- tensorflow::tf$cast(in_tensor, "float32")
#' 
#' path_gradients_one <- compute_gradients(
#'   model = model,
#'   input_tensor = in_tensor,
#'   target_class_idx = target_class_idx,
#'   input_idx = NULL,
#'   pred_stepwise = pred_stepwise)
#' 
#' dim(path_gradients_comp)
#' dim(path_gradients_one)
#' i <- 22
#' path_gradients_comp[i,,] 
#' path_gradients[i,,]
#' 
#' ##### two inputs 
#' 
#' m_path <- "/home/rmreches/checkpoints/target_middle_lm_checkpoints/Ep.035-val_loss0.23-val_acc0.948.hdf5"
#' library(magrittr)
#' input_tensor <- array(runif(1*5*4), dim = c(1,5,4))
#' input_tensor[1,,]
#' target_class_idx <- 1
#' model <- keras::load_model_hdf5(m_path, compile = FALSE)
#' model <- keras::load_model_weights_hdf5(model, m_path)
#' model
#' gen <- dummy_gen(model, 1)
#' input_tensor <- gen()[[1]]
#' input_tensor[[1]] <- tensorflow::tf$cast(input_tensor[[1]], dtype = "float32")
#' input_tensor[[2]] <- tensorflow::tf$cast(input_tensor[[2]], dtype = "float32")
#' 
#' input_idx = NULL
#' pred_stepwise = FALSE
#' m_steps = 50
#' baseline_type = "zero"
#' input_seq <- input_tensor
#' target_class_idx <- 1
#' 
#' a <- integrated_gradients(m_steps = 50,
#'                                  baseline_type = "zero",
#'                                  input_seq = input_tensor,
#'                                  target_class_idx = 1,
#'                                  model = model,
#'                                  pred_stepwise = FALSE)
#' 
#' b <- integrated_gradients(m_steps = 50,
#'                           baseline_type = "zero",
#'                           input_seq = input_tensor,
#'                           target_class_idx = 1,
#'                           model = model,
#'                           pred_stepwise = TRUE)
#' 
#' a[[2]]
#' b[[2]]
#' m_steps = 50
#' baseline_type = "zero"
#' input_seq = input_tensor
#' target_class_idx = 1
#' model = model
#' pred_stepwise = TRUE
#' 
#' 
#' 
#' #' library(reticulate)
#' #' library(magrittr)
#' #' seq_1 <- rep("T", 1000)
#' #' seq_2 <- rep("T", 1000)
#' #' index <- (1:100)*10
#' #' index
#' #' #index2 <- (1:100)*10 - 5
#' #' #index2
#' #' seq_1[index] <- "A"
#' #' seq_2[index] <- "C"
#' #' #seq_1[index2] <- "A"
#' #' #seq_2[index2] <- "C"
#' #' 
#' #' seq_1 <- paste(seq_1, collapse = "")
#' #' seq_2 <- paste(seq_2, collapse = "")
#' #' df_1 <- data.frame(Header = "A", Sequence = seq_1)
#' #' df_2 <- data.frame(Header = "C", Sequence = seq_2)
#' #' df_1$Header <- as.character(df_1$Header)
#' #' df_1$Sequence <- as.character(df_1$Sequence)
#' #' df_2$Header <- as.character(df_2$Header)
#' #' df_2$Sequence <- as.character(df_2$Sequence)
#' #' microseq::writeFasta(fdta = dplyr::as_tibble(df_1),
#' #'                      out.file =  "/home/rmreches/testData/A.fasta")
#' #' microseq::writeFasta(fdta = dplyr::as_tibble(df_2),
#' #'                      out.file =  "/home/rmreches/testData/C.fasta")
#' #' model <- create_model_lstm_cnn(maxlen = 15, layer_lstm = 32,
#' #'                                bidirectional = FALSE, learning.rate = 1e-04,
#' #'                                layer_dense = c(16, 2), vocabulary.size = 4)
#' #' trainNetwork(model = model,
#' #'              #path = "/home/rmreches/testData/LM/train",
#' #'              #path.val = "/home/rmreches/testData/LM/validation",
#' #'              path =     c("/home/rmreches/testData/A.fasta",
#' #'                           "/home/rmreches/testData/C.fasta"),
#' #'              path.val = c("/home/rmreches/testData/A.fasta",
#' #'                           "/home/rmreches/testData/C.fasta"),
#' #'              checkpoint_path = "/home/rmreches/checkpoints",
#' #'              train_type = "label_folder",
#' #'              run.name = "A_C",
#' #'              validation.split = 0.5,
#' #'              batch.size = 32,
#' #'              vocabulary = c("A", "C", "G", "T"),
#' #'              epochs = 8,
#' #'              steps.per.epoch = 50,
#' #'              step = 7,
#' #'              tensorboard.log = "/home/rmreches/tensorboard",
#' #'              output = list(
#' #'                none = TRUE,
#' #'                checkpoints = TRUE,
#' #'                tensorboard = TRUE,
#' #'                log = FALSE,
#' #'                serialize_model = FALSE,
#' #'                full_model = FALSE
#' #'              ),
#' #'              labelVocabulary = c("A", "C")
#' #' )
#' #' # optimizer <-  keras::optimizer_adam(lr = 0.001)
#' #' # model %>% keras::compile(loss = "categorical_crossentropy",
#' #' #                          optimizer = optimizer, metrics = c("acc"))
#' #' m <- matrix(0, nrow = 15, ncol = 4)
#' #' m[1:15, 4] <- 1
#' #' m[10, 1] <- 1
#' #' m
#' #' input_tensor <- tensorflow::tf$cast(array(m, dim = c(1, 15, 4)), dtype = "float32")
#' #' input_seq <- input_tensor
#' #' #input_tensor <- tensorflow::tf$cast(array(rnorm(50*4), dim = c(1, 50, 4)), dtype = "float32")
#' #' 
#' #' # TODO: more than one input layer
#' #' #' @param baseline_type Baseline sequence, either "zero" for all zeros or "shuffle" for random permutation of input_seq.
#' #' interpolate_seq <- function(m_steps = 50,
#' #'                             baseline_type = "shuffle",
#' #'                             input_seq) {
#' #'   
#' #'   stopifnot(baseline_type %in% c("zero", "shuffle"))
#' #'   if (is.list(input_seq)) {
#' #'     baseline <- list()
#' #'     for (i in 1:length(input_seq)) {
#' #'       input_dim <- dim(input_seq[[i]])
#' #'       if (baseline_type == "zero") {
#' #'         baseline[[i]] <- array(rep(0, prod(input_dim)), dim = input_dim)
#' #'       } else {
#' #'         input_dim <- dim(input_seq[[i]])
#' #'         baseline[[i]] <- array(input_seq[[i]][ , sample(input_dim[2]), ], dim = input_dim)
#' #'       }
#' #'     }
#' #'   } else {
#' #'     if (baseline_type == "zero") {
#' #'       baseline <- array(rep(0, prod(dim(input_seq))), dim = dim(input_seq))
#' #'     } else {
#' #'       baseline <- array(input_seq[ , sample(dim(input_seq)[2]), ], dim = dim(input_seq))
#' #'     }
#' #'   }
#' #'   
#' #'   m_steps <- as.integer(m_steps)
#' #'   alphas <- tensorflow::tf$linspace(start = 0.0, stop = 1.0, num = m_steps + 1L) # Generate m_steps intervals for integral_approximation() below.
#' #'   alphas_x <- alphas[ , tensorflow::tf$newaxis, tensorflow::tf$newaxis]
#' #'   if (is.list(baseline)) {
#' #'     delta <- list()
#' #'     sequences <- list()
#' #'     for (i in 1:length(baseline)) {
#' #'       delta[[i]] <- input_seq[[i]] - baseline[[i]]
#' #'       sequences[[i]] <- baseline[[i]] +  alphas_x * delta[[i]]
#' #'     }
#' #'   } else {
#' #'     delta <- input_seq - baseline
#' #'     sequences <- baseline +  alphas_x * delta
#' #'   }
#' #'   return(sequences)
#' #' }
#' #' 
#' #' baseline_seq <- interpolate_seq(m_steps = 50,
#' #'                                 baseline_type = "zero",
#' #'                                 input_seq = input_tensor)
#' #' 
#' #' 
#' #' # target_class_idx 
#' #' # input_idx, input layer to monitor for > 1 input 
#' #' compute_gradients <- function(input_tensor, target_class_idx, model, input_idx = NULL) {
#' #'   py_run_string("import tensorflow as tf")
#' #'   py$input_tensor <- input_tensor
#' #'   py$input_idx <- as.integer(input_idx - 1)
#' #'   py$target_class_idx <- as.integer(target_class_idx - 1)
#' #'   py$model <- model
#' #'   if (!is.null(input_idx)) {
#' #'     py_run_string(
#' #'       "with tf.GradientTape() as tape:
#' #'          tape.watch(input_tensor[input_idx])
#' #'          probs = model(input_tensor)[:, target_class_idx]
#' #'   ")
#' #'   } else {
#' #'     py_run_string(
#' #'       "with tf.GradientTape() as tape:
#' #'          tape.watch(input_tensor)
#' #'          probs = model(input_tensor)[:, target_class_idx]
#' #'   ")
#' #'   }  
#' #'   grad <- py$tape$gradient(py$probs, py$input_tensor)
#' #'   if (!is.null(input_idx)) {
#' #'     return(grad[[input_idx]])
#' #'   } else {
#' #'     return(grad)
#' #'   }  
#' #' }
#' #' 
#' #' w <- compute_gradients(input_tensor = input_tensor, target_class_idx = 1, model = model)
#' #' w
#' #' 
#' #' path_gradients <- compute_gradients(
#' #'   model = model,
#' #'   input_tensor = baseline_seq,
#' #'   target_class_idx = 1,
#' #'   input_idx = 1)
#' #' 
#' #' if (is.list(input_seq)) {
#' #'   path_gradients <- list()
#' #'   for (i in 1:length(input_seq)) {
#' #'     path_gradients[[i]] <- compute_gradients(
#' #'       model = model,
#' #'       input_tensor = baseline_seq,
#' #'       target_class_idx = 1,
#' #'       input_idx = i)
#' #'   }
#' #' } else {
#' #'   path_gradients <- compute_gradients(
#' #'     model = model,
#' #'     input_tensor = baseline_seq,
#' #'     target_class_idx = 1,
#' #'     input_idx = NULL)
#' #' }
#' #' 
#' #' 
#' #' integral_approximation <- function(gradients) {
#' #'   py_run_string("import tensorflow as tf")
#' #'   py$gradients <- gradients
#' #'   # riemann_trapezoidal
#' #'   py_run_string("grads = (gradients[:-1] + gradients[1:]) / tf.constant(2.0)")
#' #'   py_run_string("integrated_gradients = tf.math.reduce_mean(grads, axis=0)")
#' #'   return(py$integrated_gradients)
#' #' }
#' #' 
#' #' avg_grads <- integral_approximation(gradients = path_gradients)
#' #' 
#' #' baseline_type <- "zero"
#' #' if (baseline_type == "zero") {
#' #'   baseline <- array(rep(0, prod(dim(input_seq))), dim = dim(input_seq))
#' #' } else {
#' #'   baseline <- as.array(input_seq)[ , sample(dim(input_seq)[2]), ]
#' #' }
#' #' 
#' #' integrated_grads <- (input_seq - baseline) * avg_grads
#' #' 
#' #' integrated_gradients <- function(m_steps = 50,
#' #'                                  baseline_type = "zero",
#' #'                                  input_seq,
#' #'                                  target_class_idx,
#' #'                                  model, 
#' #'                                  input_idx = NULL
#' #' ) {
#' #'   
#' #'   baseline_seq <- interpolate_seq(m_steps = m_steps,
#' #'                                   baseline_type = baseline_type,
#' #'                                   input_seq = input_seq)
#' #'   
#' #'   if (is.list(input_seq)) {
#' #'     path_gradients <- list()
#' #'     avg_grads <- list()
#' #'     ig <- list()
#' #'     for (i in 1:length(input_seq)) {
#' #'       path_gradients[[i]] <- compute_gradients(
#' #'         model = model,
#' #'         input_tensor = baseline_seq,
#' #'         target_class_idx = target_class_idx,
#' #'         input_idx = i)
#' #'       
#' #'       avg_grads[[i]] <- integral_approximation(gradients = path_gradients[[i]])
#' #'       ig[[i]] <- ((input_seq[[i]] - baseline_seq[[i]][1, , ]) * avg_grads[[i]])[1, , ]
#' #'     }
#' #'   } else {
#' #'     path_gradients <- compute_gradients(
#' #'       model = model,
#' #'       input_tensor = baseline_seq,
#' #'       target_class_idx = target_class_idx,
#' #'       input_idx = NULL)
#' #'     
#' #'     avg_grads <- integral_approximation(gradients = path_gradients)
#' #'     ig <- ((input_seq - baseline_seq[1, , ]) * avg_grads)[1, , ]
#' #'   }
#' #'   
#' #'   return(ig)
#' #' }  
#' #' 
#' #' integrated_grads <- integrated_gradients(m_steps = 50,
#' #'                            baseline_type = "zero",
#' #'                            input_seq = input_seq,
#' #'                            target_class_idx = 1,
#' #'                            model = model, 
#' #'                            input_idx = NULL
#' #' ) 
#' #' integrated_grads
#' #' 
#' #' #### heatmap ######
#' #' 
#' #' heatmaps_integrated_grad <- function(integrated_grads,
#' #'                                      input_seq) {
#' #'   
#' #'   if (is.list(input_seq)) {
#' #'     num_input <- length(input_seq)
#' #'     attribution_mask <- list()
#' #'     nuc_matrix <- list()
#' #'     nuc_seq <- list()
#' #'     sum_nuc <- list()
#' #'     for (i in 1:length(input_seq)) {
#' #'       py$integrated_grads <- integrated_grads[[i]]
#' #'       py_run_string("attribution_mask = tf.reduce_sum(tf.math.abs(integrated_grads), axis=-1)")
#' #'       py_run_string("sum_nuc = tf.reduce_sum(integrated_grads, axis=-1)")
#' #'       attribution_mask[[i]] <- py$attribution_mask
#' #'       attribution_mask[[i]] <- as.matrix(attribution_mask[[i]], nrow = 1) %>% as.data.frame()
#' #'       colnames(attribution_mask[[i]]) <- "abs_sum"
#' #'       sum_nuc[[i]] <- py$sum_nuc 
#' #'       sum_nuc[[i]] <- as.matrix(sum_nuc[[i]], nrow = 1) %>% as.data.frame()
#' #'       colnames(sum_nuc[[i]]) <- "sum"
#' #'       
#' #'       nuc_matrix[[i]] <- as.matrix(integrated_grads[[i]])
#' #'       nuc_seq[[i]] <- apply(input_seq[[i]][1, ,], 1, which.max) %>% as.character()
#' #'       nuc_seq[[i]] <- nuc_seq[[i]] %>% stringr::str_replace_all("1", "A") %>%
#' #'         stringr::str_replace_all("2", "C") %>%
#' #'         stringr::str_replace_all("3", "G") %>%
#' #'         stringr::str_replace_all("4", "T")
#' #'       rownames(nuc_matrix[[i]]) <- nuc_seq[[i]]
#' #'       colnames(nuc_matrix[[i]]) <- c("A", "C", "G", "T")
#' #'     }
#' #'     
#' #'   } else {
#' #'     num_input <- 1
#' #'     py$integrated_grads <- integrated_grads
#' #'     py_run_string("attribution_mask = tf.reduce_sum(tf.math.abs(integrated_grads), axis=-1)")
#' #'     py_run_string("sum_nuc = tf.reduce_sum(integrated_grads, axis=-1)")
#' #'     #py_run_string("mean_nuc = tf.reduce_mean(integrated_grads, axis=-1)")
#' #'     
#' #'     attribution_mask <- py$attribution_mask
#' #'     attribution_mask <- as.matrix(attribution_mask, nrow = 1) %>% as.data.frame()
#' #'     colnames(attribution_mask) <- "abs_sum"
#' #'     
#' #'     sum_nuc <- py$sum_nuc 
#' #'     sum_nuc <- as.matrix(sum_nuc, nrow = 1) %>% as.data.frame()
#' #'     colnames(sum_nuc) <- "sum"
#' #'     
#' #'     nuc_matrix <- as.matrix(integrated_grads)
#' #'     nuc_seq <- apply(input_seq[1, ,], 1, which.max) %>% as.character()
#' #'     nuc_seq <- nuc_seq %>% stringr::str_replace_all("1", "A") %>%
#' #'       stringr::str_replace_all("2", "C") %>%
#' #'       stringr::str_replace_all("3", "G") %>%
#' #'       stringr::str_replace_all("4", "T")
#' #'     rownames(nuc_matrix) <- nuc_seq
#' #'     colnames(nuc_matrix) <- c("A", "C", "G", "T")
#' #'     
#' #'   }
#' #'   
#' #'   ig_min <- Inf 
#' #'   ig_max <- -Inf
#' #'   
#' #'   if (num_input == 1) {
#' #'     ig_min <- min(c(ig_min, keras::k_min(integrated_grads)$numpy())) 
#' #'     ig_max <- max(c(ig_max, keras::k_max(integrated_grads)$numpy()))
#' #'   } else {
#' #'     for (i in 1:num_input) {
#' #'       ig_min <- min(c(ig_min, keras::k_min(integrated_grads[[i]])$numpy())) 
#' #'       ig_max <- max(c(ig_max, keras::k_max(integrated_grads[[i]])$numpy()))
#' #'     }
#' #'   }
#' #'   col_fun <- circlize::colorRamp2(c(ig_min, mean(c(ig_max, ig_min)) , ig_max), c("blue", "white", "red"))
#' #'   
#' #'   hm_list <- list()
#' #'   if (num_input == 1) {
#' #'     row_ha = columnAnnotation(abs_sum = attribution_mask[,1], sum = sum_nuc[,1], mean = mean_nuc[,1])
#' #'     
#' #'     hm_list[[1]] <- Heatmap(matrix = t(nuc_matrix),
#' #'                        name = "hm",
#' #'                        top_annotation = row_ha,
#' #'                        col = col_fun,
#' #'                        cluster_rows = FALSE,
#' #'                        cluster_columns = FALSE,
#' #'                        column_names_rot = 0
#' #'     )
#' #'   } else {
#' #'     for (i in 1:num_input) {
#' #'       row_ha = columnAnnotation(abs_sum = attribution_mask[[i]][,1], sum = sum_nuc[[i]][,1])
#' #'       hm_list[[i]] <- Heatmap(matrix = t(nuc_matrix[[i]]),
#' #'                          name = paste0("hm_", i),
#' #'                          top_annotation = row_ha,
#' #'                          col = col_fun,
#' #'                          cluster_rows = FALSE,
#' #'                          cluster_columns = FALSE,
#' #'                          column_names_rot = 0
#' #'       )
#' #'     }
#' #'   }
#' #'   hm_list
#' #' }
#' #' 
#' #' py$integrated_grads <- integrated_grads
#' #' py_run_string("attribution_mask = tf.reduce_sum(tf.math.abs(integrated_grads), axis=-1)")
#' #' attribution_mask <- py$attribution_mask
#' #' 
#' #' m <- as.matrix(integrated_grads[1,,])
#' #' nuc_seq <- apply(input_seq[1, ,], 1, which.max) %>% as.character()
#' #' nuc_seq <- nuc_seq %>% stringr::str_replace_all("1", "A") %>%
#' #'   stringr::str_replace_all("2", "C") %>%
#' #'   stringr::str_replace_all("3", "G") %>%
#' #'   stringr::str_replace_all("4", "T")
#' #' nuc_seq
#' #' rownames(m) <- nuc_seq
#' #' colnames(m) <- c("A", "C", "G", "T")
#' #' 
#' #' col_fun <- circlize::colorRamp2(c(min(m), max(m)), c("blue", "white", "red"))
#' #' hm <- Heatmap(matrix = m,
#' #'               name = "hm",
#' #'               col = col_fun,
#' #'               cluster_rows = FALSE,
#' #'               cluster_columns = FALSE
#' #' )
#' #' hm
#' #' 
#' #' m <- as.matrix(attribution_mask)
#' #' colnames(m) <- nuc_seq
#' #' 
#' #' col_fun <- circlize::colorRamp2(c(min(m), max(m)), c("white", "red"))
#' #' hm <- Heatmap(matrix = m,
#' #'               name = "hm",
#' #'               col = col_fun,
#' #'               cluster_rows = FALSE,
#' #'               cluster_columns = FALSE,
#' #'               column_names_rot = 90
#' #' )
#' #' hm
