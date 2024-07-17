#' F1 metric
#' 
#' Compute F1 metric. If loss is `"categorical_crossentropy"`, number of targets must be 2. If
#' loss is `"binary_crossentropy"` and number of targets > 1, will flatten `y_true` and `y_pred` matrices 
#' to a single vector (rather than computing separate F1 scores for each class).
#'
#' @param num_targets Size of model output.
#' @param loss Loss function of model.
#' @examples 
#' y_true <- c(1,0,0,1,1,0,1,0,0)  
#' y_pred <-  c(0.9,0.05,0.05,0.9,0.05,0.05,0.9,0.05,0.05)
#'
#' f1_metric <- f1_wrapper(3L, "binary_crossentropy")
#' f1_metric$update_state(y_true, y_pred)
#' f1_metric$result()  
#' 
#' # add metric to a model
#' 
#' num_targets <- 1
#' model <- create_model_lstm_cnn(maxlen = 20,
#'                                layer_lstm = 8,
#'                                bal_acc = FALSE,
#'                                last_layer_activation = "sigmoid",
#'                                loss_fn = "binary_crossentropy",
#'                                layer_dense = c(8, num_targets))
#' 
#' f1_metric <- f1_wrapper(num_targets, loss = model$loss)
#' model %>% keras::compile(loss = model$loss, 
#'                          optimizer = model$optimizer,
#'                          metrics = c(model$metrics, f1_metric))
#' @returns A keras metric.                          
#' @export
f1_wrapper <- function(num_targets = 2, loss = "binary_crossentropy") {
  
  stopifnot(loss %in% c("binary_crossentropy", "categorical_crossentropy"))
  
  if (loss == "binary_crossentropy" & num_targets > 1) {
    message("Will flatten y_true and y_pred matrices to a single vector for evaluation
            rather than computing separate F1 scores for each class and taking the mean.")
  }
  
  if (loss == "categorical_crossentropy" & num_targets != 2) {
    stop("Output size must be two, when loss is categorical_crossentropy")
  }
  
  f1_stateful <- reticulate::PyClass("f1",
                                     inherit = tensorflow::tf$keras$metrics$Metric,
                                     list(
                                       
                                       `__init__` = function(self, num_targets, loss) {
                                         super()$`__init__`(name = "f1")
                                         self$num_targets <- num_targets
                                         self$f1_score <- 0
                                         self$loss <- loss
                                         self$rc <- tensorflow::tf$keras$metrics$Recall()
                                         self$pr <- tensorflow::tf$keras$metrics$Precision()
                                         NULL
                                       },
                                       
                                       update_state = function(self, y_true, y_pred, sample_weight = NULL) {
                                         if (self$loss == "binary_crossentropy") {
                                           self$rc$update_state(y_true, y_pred)
                                           self$pr$update_state(y_true, y_pred)
                                         } else {
                                           self$rc$update_state(y_true[ , 1], y_pred[ , 1])
                                           self$pr$update_state(y_true[ , 1], y_pred[ , 1])
                                         }
                                         NULL
                                       },
                                       
                                       result = function(self) {
                                         self$f1_score <- self$compute_f1()
                                         return(self$f1_score)
                                       },
                                       
                                       compute_f1 = function(self) {
                                         f1 <- (2 * self$pr$result() * self$rc$result())/(self$pr$result() + self$rc$result() + tensorflow::tf$constant(1e-15))
                                         return(f1)
                                       },
                                       
                                       reset_state = function(self) {
                                         self$rc$reset_state()
                                         self$pr$reset_state()
                                         #self$f1_score$assign(0)
                                         NULL
                                       }
                                       
                                     ))
  
  return(f1_stateful(num_targets = num_targets, loss = loss))
}


#' Balanced accuracy metric
#'
#' Compute balanced accuracy as additional score. Useful for imbalanced data. Only implemented for 
#' model with mutually exclusive targets.
#'
#' @param num_targets Number of targets.
#' @param cm_dir Directory of confusion matrix used to compute balanced accuracy.
#' @examples 
#' 
#' y_true <- c(1,0,0,1,
#'             0,1,0,0,
#'             0,0,1,0) %>% matrix(ncol = 3)
#' y_pred <- c(0.9,0.1,0.2,0.1,
#'             0.05,0.7,0.2,0.0,
#'             0.05,0.2,0.6,0.9) %>% matrix(ncol = 3)
#' 
#' cm_dir <- tempfile() 
#' dir.create(cm_dir)
#' bal_acc_metric <- balanced_acc_wrapper(num_targets = 3L, cm_dir = cm_dir)
#' bal_acc_metric$update_state(y_true, y_pred)
#' bal_acc_metric$result()
#' as.array(bal_acc_metric$cm)
#'
#' @returns A keras metric.                          
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
                                                   #self$store_cm()
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
                                                   balanced_acc <- tensorflow::tf$math$divide(acc_sum, tensorflow::tf$math$count_nonzero(col_sums, dtype= acc_sum$dtype))
                                                   return(balanced_acc)
                                                 },
                                                 
                                                 reset_state = function(self) {
                                                   self$store_cm()
                                                   self$count <- self$count + 1
                                                   self$cm$assign_sub(self$cm)
                                                   NULL
                                                 },
                                                 
                                                 store_cm = function(self) {
                                                   #if (self$count > 0) {
                                                   if (self$count %% 2 != 0) {
                                                     file_name <- file.path(self$cm_dir, paste0("cm_val_", floor(self$count/2), ".rds"))
                                                   } else {
                                                     file_name <- file.path(self$cm_dir, paste0("cm_train_", floor(self$count/2), ".rds"))
                                                   }
                                                   saveRDS(keras::k_eval(self$cm), file_name)
                                                   NULL
                                                   #}
                                                 }
                                                 
                                               ))
  return(balanced_acc_stateful(num_targets = num_targets, cm_dir = cm_dir))
}


#' Mean AUC score
#'
#' Compute AUC score as additional metric. If model has several output neurons with binary crossentropy loss, will use the average score.
#'
#' @param model_output_size Number of neurons in model output layer.
#' @param loss Loss function of model, for which metric will be applied to; must be `"binary_crossentropy"`
#' or `"categorical_crossentropy"`.
#' @examples 
#' y_true <- c(1,0,0,1,1,0,1,0,0) %>% matrix(ncol = 3)
#' y_pred <- c(0.9,0.05,0.05,0.9,0.05,0.05,0.9,0.05,0.05) %>% matrix(ncol = 3)
#' 
#' auc_metric <- auc_wrapper(3L, "binary_crossentropy")
#' 
#' auc_metric$update_state(y_true, y_pred)
#' auc_metric$result()  
#' 
#' # add metric to a model
#' num_targets <- 4
#' model <- create_model_lstm_cnn(maxlen = 20,
#'                                layer_lstm = 8,
#'                                bal_acc = FALSE,
#'                                last_layer_activation = "sigmoid",
#'                                loss_fn = "binary_crossentropy",
#'                                layer_dense = c(8, num_targets))
#' 
#' auc_metric <- auc_wrapper(num_targets, loss = model$loss)
#' model %>% keras::compile(loss = model$loss, 
#'                          optimizer = model$optimizer,
#'                          metrics = c(model$metrics, auc_metric))
#'
#' @returns A keras metric.                          
#' @export
auc_wrapper <- function(model_output_size,
                        loss = "binary_crossentropy") {
  
  multi_label <- FALSE
  stopifnot(loss %in% c("binary_crossentropy", "categorical_crossentropy"))
  
  if (loss == "categorical_crossentropy" & model_output_size != 2) {
    stop("Output size must be two, when loss is categorical_crossentropy")
  }
  
  if (loss == "categorical_crossentropy") {
    label_weights <- c(1L, 0L)
  } else {
    label_weights <- NULL
  }
  
  if (loss == "binary_crossentropy" & model_output_size > 1) {
    multi_label <- TRUE
  }
  
  metric_name <- ifelse(loss == "binary_crossentropy" & model_output_size > 1,
                        "mean_AUC", "AUC")
  
  auc_metric <- tensorflow::tf$keras$metrics$AUC(label_weights = label_weights,
                                                 multi_label = multi_label)
  
  return(auc_metric)
}

#' Loss function for label noise
#' 
#' Implements approach from this [paper](https://arxiv.org/abs/1609.03683) and code from 
#' [here](https://github.com/giorgiop/loss-correction/blob/15a79de3c67c31907733392085c333547c2f2b16/loss.py#L16-L21). 
#' Can be used if labeled data contains noise, i.e. some of the data is labeled wrong.
#'
#' @param noise_matrix Matrix of noise distribution.
#' @importFrom magrittr "%>%"
#' @examples 
#' # If first label contains 5% wrong labels and second label no noise
#' noise_matrix <- matrix(c(0.95, 0.05, 0, 1), nrow = 2, byrow = TRUE)
#' noisy_loss <- noisy_loss_wrapper(noise_matrix)
#' 
#' @returns A function implementing noisy loss.                          
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

cpcloss <- function(latents,
                    context,
                    target_dim = 64,
                    emb_scale = 0.1 ,
                    steps_to_ignore = 2,
                    steps_to_predict = 3,
                    steps_skip = 1,
                    batch_size = 32,
                    k = 5,
                    train_type = "cpc") {
  # define empty lists for metrics
  loss <- list()
  acc <- list()
  # create context tensor
  ctx <- context(latents)
  c_dim <- latents$shape[[2]]
  
  # loop for different distances of predicted patches
  for (i in seq(steps_to_ignore, (steps_to_predict - 1), steps_skip)) {
    # define patches to be deleted
    c_dim_i <- c_dim - i - 1
    if (train_type == "Self-GenomeNet") {
      steps_to_ignore <- 1
      steps_to_predict <- 2
      steps_skip <- 1
      target_dim <- ctx$shape[[3]]
      ctx_conv <-
        ctx %>% keras::layer_conv_1d(kernel_size = 1, filters = target_dim)
      logits <- tf$zeros(list(0L,as.integer(batch_size*2)))
      for (j in seq_len(c_dim - (i + 1))) {
        basepos <- ctx_conv[, j,] %>% keras::k_reshape(c(-1, target_dim))
        targetpos <-
          ctx[, (c_dim - j - i), ] %>% keras::k_reshape(c(-1, target_dim))
        logits_j <- tensorflow::tf$matmul(basepos, tensorflow::tf$transpose(targetpos))
        logits <- tensorflow::tf$concat(list(logits, logits_j), axis = 0L)
      }
      # define labels
      labels <-
        rep(c(seq(batch_size, (
          2 * batch_size - 1
        )), (seq(
          0, (batch_size - 1)
        ))), (c_dim - (i + 1))) %>% as.integer()
    } else {
      c_dim_i <- c_dim - i - 1
      # define total number of elements in context tensor
      total_elements <- batch_size * c_dim_i
      # add conv layer and reshape tensor for matrix multiplication
      targets <-
        latents %>% keras::layer_conv_1d(kernel_size = 1, filters = target_dim) %>% keras::k_reshape(c(-1, target_dim))
      # add conv layer and reshape for matrix multiplication
      preds_i <-
        ctx %>% keras::layer_conv_1d(kernel_size = 1, filters = target_dim)
      preds_i <- preds_i[, (1:(c_dim - i - 1)),]
      preds_i <- keras::k_reshape(preds_i, c(-1, target_dim)) * emb_scale
      
      # define logits normally
      logits <- tensorflow::tf$matmul(preds_i, tensorflow::tf$transpose(targets))
      
      # get position of labels
      b <- floor(seq(0, total_elements - 1) / c_dim_i)
      col <- seq(0, total_elements - 1) %% c_dim_i
      
      # define labels
      labels <- b * c_dim + col + (i + 1)
      labels <- as.integer(labels)
      
    }
    # calculate loss and accuracy for each step
    loss[[length(loss) + 1]] <-
      tensorflow::tf$nn$sparse_softmax_cross_entropy_with_logits(labels, logits) %>%
      tensorflow::tf$stack(axis = 0) %>% tensorflow::tf$reduce_mean()
    acc[[length(acc) + 1]] <-
      tensorflow::tf$keras$metrics$sparse_top_k_categorical_accuracy(tensorflow::tf$cast(labels, dtype = "int64"), logits, as.integer(k)) %>%
      tensorflow::tf$stack(axis = 0) %>% tensorflow::tf$reduce_mean()
  }
  # convert to tensor for output
  loss <- loss %>% tensorflow::tf$stack(axis = 0) %>% tensorflow::tf$reduce_mean()
  acc <- acc %>% tensorflow::tf$stack(axis = 0) %>% tensorflow::tf$reduce_mean()
  return(tensorflow::tf$stack(list(loss, acc)))
}

#' Stochastic Gradient Descent with Warm Restarts
#' 
#' Compute the learning Rate for a given epoch using Stochastic Gradient Descent with Warm Restarts. Implements approach from this [paper](https://arxiv.org/abs/1608.03983).
#'
#' @param lrmin Lower limit of the range for the learning rate.
#' @param lrmax Upper limit of the range for the learning rate.
#' @param restart Number of epochs until a restart is conducted.
#' @param mult Factor, by which the number of epochs until a restart is increased at every restart.
#' @param epoch Epoch, for which the learning rate shall be calculated.
#' @examples 
#' sgdr(lrmin = 5e-10, lrmax = 5e-2, restart = 50,
#' mult = 1, epoch = 5)
#' 
#' @returns A numeric value.
#' @export
sgdr <- function(lrmin = 5e-10,
                 lrmax = 5e-2,
                 restart = 50,
                 mult = 1,
                 epoch = NULL) {
  iter <- c()
  position <- c()
  i <- 0
  while (length(iter) < epoch) {
    iter <- c(iter, rep(i, restart * mult ^ i))
    position <- c(position, c(1:(restart * mult ^ i)))
    i <- i + 1
  }
  restart2 <- (restart * mult ^ iter[epoch])
  epoch <- position[epoch]
  return(lrmin + 1 / 2 * (lrmax - lrmin) * (1 + cos((epoch / restart2) * pi)))
}

#' Step Decay
#' 
#' Compute the learning Rate for a given epoch using Step Decay.
#'
#' @param lrmax Upper limit of the range for the learning rate.
#' @param newstep Number of epochs until the learning rate is reduced.
#' @param mult Factor, by which the number of epochs until a restart is decreased after a new step.
#' @param epoch Epoch, for which the learning rate shall be calculated.
#' @examples
#' stepdecay(lrmax = 0.005, newstep = 50,
#' mult = 0.7, epoch = 3)
#' 
#' @returns A numeric value.
#' @export
stepdecay <- function(lrmax = 0.005,
                      newstep = 50,
                      mult = 0.7,
                      epoch = NULL) {
  return(lrmax * (mult ^ (floor((
    epoch
  ) / newstep))))
}

#' Exponential Decay
#' 
#' Compute the learning Rate for a given epoch using Exponential Decay.
#'
#' @param lrmax Upper limit of the range for the learning rate.
#' @param mult Factor, by which the number of epochs until a restart is decreased after a new step.
#' @param epoch Epoch, for which the learning rate shall be calculated.
#' @examples
#' exp_decay(lrmax = 0.005, mult = 0.1, epoch = 8) 
#' 
#' @returns A numeric value.
#' @export
exp_decay <- function(lrmax = 0.005,
                      mult = 0.1,
                      epoch = NULL) {
  return(lrmax * exp(-mult * epoch))
}


euclidean_distance <- function(vects) {
  x <- vects[[1]]
  y <- vects[[2]]
  sum_square <- tensorflow::tf$math$reduce_sum(tensorflow::tf$math$square(x - y), axis=1L, keepdims=TRUE)
  return(tensorflow::tf$math$sqrt(tensorflow::tf$math$maximum(sum_square, tensorflow::tf$keras$backend$epsilon())))
}

cosine_similarity <- function(vects) {
  x <- vects[[1]]
  y <- vects[[2]]
  xy_dot <- tensorflow::tf$math$reduce_sum(x*y, axis=1L, keepdims=TRUE)
  x_norm <- tensorflow::tf$math$sqrt(tensorflow::tf$math$reduce_sum(tensorflow::tf$math$square(x), axis=1L, keepdims=TRUE))
  y_norm <- tensorflow::tf$math$sqrt(tensorflow::tf$math$reduce_sum(tensorflow::tf$math$square(y), axis=1L, keepdims=TRUE))
  return(xy_dot/(x_norm*y_norm))
}

#' Contrastive loss 
#' 
#' Contrastive loss as used here: https://keras.io/examples/vision/siamese_contrastive/.
#'
#' @param margin Integer, baseline for distance for which pairs should be classified as dissimilar.
#' @examples
#' cl <- loss_cl(margin=1)
#' 
#' @returns A function implementing contrastive loss.
#' @export
loss_cl <- function(margin=1) {
  
  contrastive_loss <- function(y_true, y_pred) {
    
    square_pred <- tensorflow::tf$math$square(y_pred)
    margin_square <- tensorflow::tf$math$square(tensorflow::tf$math$maximum(margin - (y_pred), 0))
    l <- tensorflow::tf$math$reduce_mean(
      (1 - y_true) * square_pred + (y_true) * margin_square
    )
    return(l)
  }
  
  return(contrastive_loss)
  
}

#' Focal loss for two or more labels
#' 
#' @param y_true Vector of true values.
#' @param y_pred Vector of predicted values.
#' @param gamma Focusing parameter.
#' @param alpha Vector of weighting factors.
#' @examples
#' y_true <- matrix(c(0, 1, 0, 0, 0, 1), nrow = 2, byrow = TRUE)
#' y_pred <- matrix(c(0.15, 0.8, 0.05,
#'                    0.08, 0.02, 0.9), nrow = 2, byrow = TRUE) 
#' fl <- focal_loss_multiclass(y_true, y_pred)
#' fl$numpy()
#' 
#' @returns A function implementing focal loss.
#' @export
focal_loss_multiclass <- function(y_true, y_pred, gamma = 2.5, alpha = c(1)) {
  y_pred <- keras::k_clip(y_pred, keras::k_epsilon(), 1.0 - keras::k_epsilon())
  cd_loss <- -y_true * keras::k_log(y_pred) # categorical cross entropy
  fl_loss <- alpha * keras::k_pow(1. - y_pred, gamma) * cd_loss
  return(keras::k_sum(fl_loss, axis = -1))
}
