#' F1 metric
#' 
#' Compute F1 metric. If loss is "categorical_crossentropy", number of targets must be 2. If
#' loss is "binary_crossentropy" and number of targets > 1, will flatten y_true and y_pred matrices 
#' to a single vector (rather than computing separate F1 scores for each class).
#'
#'@param num_targets Size of model output.
#'@param loss Loss function of model.
#'@examples 
#' y_true <- c(1,0,0,1,1,0,1,0,0)  
#' y_pred <- c(0.9,0.1,0.1,0.9,0.1,0.1,0.9,0.1,0.1)
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
#'@export
f1_wrapper <- function(num_targets = 2, loss = "binary_crossentropy") {
  
  stopifnot(loss %in% c("binary_crossentropy", "categorical_crossentropy"))
  
  if (loss == "binary_crossentropy" & num_targets > 1) {
    warning("Will flatten y_true and y_pred matrices to a single vector for evaluation
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
#' Compute balanced accuracy as additional score. Useful for imbalanced data.
#'
#'@param num_targets Number of targets.
#'@param cm_dir Directory of confusion matrix used to compute balanced accuracy.
#'
#'@export
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
                                                   if (self$count > 0) {
                                                     if (self$count %% 2 == 0) {
                                                       file_name <- file.path(self$cm_dir, paste0("cm_val_", floor(self$count/2), ".rds"))
                                                     } else {
                                                       file_name <- file.path(self$cm_dir, paste0("cm_train_", floor(self$count/2), ".rds"))
                                                     }
                                                     saveRDS(keras::k_eval(self$cm), file_name)
                                                     NULL
                                                   }
                                                 }
                                                 
                                               ))
  return(balanced_acc_stateful(num_targets = num_targets, cm_dir = cm_dir))
}


#' Mean AUC score
#'
#' Compute AUC score as additional metric. If model has several output neurons with binary crossentropy loss, will use the average score.
#'
#' @param model_output_size Number of neurons in model output layer.
#' @param loss Loss function of model, for which metric will be applied to; must be "binary_crossentropy"
#' or "categorical_crossentropy".
#' @examples 
#' y_true <- c(1,0,0,1,1,0,1,0,0) %>% matrix(ncol = 3)
#' y_pred <- c(0.9,0.1,0.1,0.9,0.1,0.1,0.9,0.1,0.1) %>% matrix(ncol = 3)
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
#'
#' @param noise_matrix Matrix of noise distribution.
#' @importFrom magrittr "%>%"
#' @examples 
#' # If first label contains 5\% wrong labels and second label no noise
#' noise_matrix <- matrix(c(0.95, 0.05, 0, 1), nrow = 2, byrow = TRUE)
#' noisy_loss <- noisy_loss_wrapper(noise_matrix)
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
