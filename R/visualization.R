#' Interpolation between baseline and prediction
#'
#' @param baseline_type Baseline sequence, either "zero" for all zeros or "shuffle" for random permutation of input_seq.
#' @param m_steps Number of steps between baseline and original input.
#' @param input_seq Input tensor.
#' @keywords internal
interpolate_seq <- function(m_steps = 50,
                            baseline_type = "shuffle",
                            input_seq) {

  stopifnot(baseline_type %in% c("zero", "shuffle"))
  if (is.list(input_seq)) {
    baseline <- list()
    for (i in 1:length(input_seq)) {
      input_dim <- dim(input_seq[[i]])
      if (baseline_type == "zero") {
        baseline[[i]] <- array(rep(0, prod(input_dim)), dim = input_dim)
      } else {
        input_dim <- dim(input_seq[[i]])
        baseline[[i]] <- array(input_seq[[i]][ , sample(input_dim[2]), ], dim = input_dim)
      }
    }
  } else {
    if (baseline_type == "zero") {
      baseline <- array(rep(0, prod(dim(input_seq))), dim = dim(input_seq))
    } else {
      baseline <- array(input_seq[ , sample(dim(input_seq)[2]), ], dim = dim(input_seq))
    }
  }

  m_steps <- as.integer(m_steps)
  alphas <- tensorflow::tf$linspace(start = 0.0, stop = 1.0, num = m_steps + 1L) # Generate m_steps intervals for integral_approximation() below.
  alphas_x <- alphas[ , tensorflow::tf$newaxis, tensorflow::tf$newaxis]
  if (is.list(baseline)) {
    delta <- list()
    sequences <- list()
    for (i in 1:length(baseline)) {
      delta[[i]] <- input_seq[[i]] - baseline[[i]]
      sequences[[i]] <- baseline[[i]] +  alphas_x * delta[[i]]
    }
  } else {
    delta <- input_seq - baseline
    sequences <- baseline +  alphas_x * delta
  }
  return(sequences)
}

#' Compute gradients
#'
#' @param input_idx  Input layer to monitor for > 1 input.
#' @param target_class_idx Index of class to compute gradient for.
#' @param model Model to compute gradient for.
#' @param pred_stepwise Whether to do predictions with batch_size 1 rather than all at once. Can be used if
#' input is too big to handle at once.
#' @keywords internal
compute_gradients <- function(input_tensor, target_class_idx, model, input_idx = NULL, pred_stepwise = FALSE) {

  # if (is.list(input_tensor)) {
  #   stop("Stepwise predictions only supported for single input layer yet")
  # }

  py_run_string("import tensorflow as tf")
  py$input_tensor <- input_tensor
  py$input_idx <- as.integer(input_idx - 1)
  py$target_class_idx <- as.integer(target_class_idx - 1)
  py$model <- model

  if (!is.null(input_idx)) {
    py_run_string(
      "with tf.GradientTape() as tape:
             tape.watch(input_tensor[input_idx])
             probs = model(input_tensor)[:, target_class_idx]
    ")
  } else {
    py_run_string(
      "with tf.GradientTape() as tape:
             tape.watch(input_tensor)
             probs = model(input_tensor)[:, target_class_idx]
    ")
  }

  grad <- py$tape$gradient(py$probs, py$input_tensor)
  if (!is.null(input_idx)) {
    return(grad[[input_idx]])
  } else {
    return(grad)
  }
}

integral_approximation <- function(gradients) {
  py_run_string("import tensorflow as tf")
  py$gradients <- gradients
  # riemann_trapezoidal
  py_run_string("grads = (gradients[:-1] + gradients[1:]) / tf.constant(2.0)")
  py_run_string("integrated_gradients = tf.math.reduce_mean(grads, axis=0)")
  return(py$integrated_gradients)
}

#' Compute integrated gradients 
#' 
#' Computes integrated gradients scores for model and an input sequence.
#' This can be used to visualize what part of the input is import for the models decision.
#' Code is R implementation of python code from [here](https://www.tensorflow.org/tutorials/interpretability/integrated_gradients).
#' Tensorflow implementation is based on this [paper](https://arxiv.org/abs/1703.01365).
#' 
#' @param baseline_type Baseline sequence, either `"zero"` for all zeros or `"shuffle"` for random permutation of `input_seq`.
#' @param m_steps Number of steps between baseline and original input.
#' @param input_seq Input tensor.
#' @param target_class_idx Index of class to compute gradient for
#' @param model Model to compute gradient for.
#' @param pred_stepwise Whether to do predictions with batch size 1 rather than all at once. Can be used if
#' input is too big to handle at once. Only supported for single input layer.
#' @param num_baseline_repeats Number of different baseline estimations if baseline_type is `"shuffle"` (estimate integrated
#' gradient repeatedly for different shuffles). Final result is average of \code{num_baseline} single calculations.
#' @examples 
#' model <- create_model_lstm_cnn(layer_lstm = 8, layer_dense = 3, maxlen = 20, verbose = FALSE)
#' random_seq <- sample(0:3, 20, replace = TRUE)
#' input_seq <- array(keras::to_categorical(random_seq), dim = c(1, 20, 4))
#' integrated_gradients(
#'   input_seq = input_seq,
#'   target_class_idx = 3,
#'   model = model)
#' @export
integrated_gradients <- function(m_steps = 50,
                                 baseline_type = "zero",
                                 input_seq,
                                 target_class_idx,
                                 model,
                                 pred_stepwise = FALSE,
                                 num_baseline_repeats = 1) {

  library(reticulate)
  py_run_string("import tensorflow as tf")
  input_idx <- NULL
  if (num_baseline_repeats > 1 & baseline_type == "zero") {
    warning('Ignoring num_baseline_repeats if baseline is of type "zero". Did you mean to use baseline_type = "shuffle"?')
  }

  if (num_baseline_repeats == 1 | baseline_type == "zero") {

    baseline_seq <- interpolate_seq(m_steps = m_steps,
                                    baseline_type = baseline_type,
                                    input_seq = input_seq)

    if (is.list(baseline_seq)) {
      for (i in 1:length(baseline_seq)) {
        baseline_seq[[i]] <- tensorflow::tf$cast(baseline_seq[[i]], dtype = "float32")
      }
    } else {
      baseline_seq <- tensorflow::tf$cast(baseline_seq, dtype = "float32")
    }

    if (is.list(input_seq)) {
      path_gradients <- list()
      avg_grads <- list()
      ig <- list()

      if (pred_stepwise) {
        path_gradients <- gradients_stepwise(
          model = model,
          baseline_seq = baseline_seq,
          target_class_idx = target_class_idx)
      } else {

        path_gradients <- compute_gradients(
          model = model,
          input_tensor = baseline_seq,
          target_class_idx = target_class_idx,
          input_idx = NULL,
          pred_stepwise = pred_stepwise)
      }

      for (i in 1:length(input_seq)) {
        avg_grads[[i]] <- integral_approximation(gradients = path_gradients[[i]])
        ig[[i]] <- ((input_seq[[i]] - baseline_seq[[i]][1, , ]) * avg_grads[[i]])[1, , ]
      }
    } else {

      if (pred_stepwise) {
        path_gradients <- gradients_stepwise(model = model,
                                             baseline_seq = baseline_seq,
                                             target_class_idx = target_class_idx,
                                             input_idx = NULL)
      } else {
        path_gradients <- compute_gradients(
          model = model,
          input_tensor = baseline_seq,
          target_class_idx = target_class_idx,
          input_idx = NULL,
          pred_stepwise = pred_stepwise)
      }

      avg_grads <- integral_approximation(gradients = path_gradients)
      ig <- ((input_seq - baseline_seq[1, , ]) * avg_grads)[1, , ]
    }
  } else {
    ig_list <- list()
    for (i in 1:num_baseline_repeats) {
      ig_list[[i]] <- integrated_gradients(m_steps = m_steps,
                                           baseline_type = "shuffle",
                                           input_seq = input_seq,
                                           target_class_idx = target_class_idx,
                                           model = model,
                                           pred_stepwise = pred_stepwise,
                                           num_baseline_repeats = 1)
    }
    ig_stacked <- tensorflow::tf$stack(ig_list, axis = 0L)
    ig <- tensorflow::tf$reduce_mean(ig_stacked, axis = 0L)
  }

  return(ig)
}

#' Compute gradients stepwise (one batch at a time)
#'
#' @keywords internal
gradients_stepwise <- function(model = model, baseline_seq, target_class_idx,
                               input_idx = NULL) {

  if (is.list(baseline_seq)) {
    first_dim <- dim(baseline_seq[[1]])[1]
    num_input_layers <- length(baseline_seq)

    l <- list()
    for (j in 1:first_dim) {
      input_list <- list()
      for (k in 1:length(baseline_seq)) {
        input <- as.array(baseline_seq[[k]][j, , ])
        input <- array(input, dim = c(1, dim(baseline_seq[[k]])[-1]))
        input <- tensorflow::tf$cast(input, baseline_seq[[k]]$dtype)
        input_list[[k]] <- input
      }
      output <- compute_gradients(
        model = model,
        input_tensor = input_list,
        target_class_idx = target_class_idx,
        input_idx = NULL)
      for (m in 1:length(output)) {
        output[[m]] <- tensorflow::tf$squeeze(output[[m]])
      }
      l[[j]] <- output
    }

    path_gradients <- vector("list", num_input_layers)
    for (n in 1:num_input_layers) {
      temp_list <- vector("list", first_dim)
      for (p in 1:first_dim){
        temp_list[[p]] <- l[[p]][[n]]
      }
      path_gradients[[n]] <- tensorflow::tf$stack(temp_list)
    }

  } else {
    l <- list()
    for (j in 1:dim(baseline_seq)[1]) {
      input <- as.array(baseline_seq[j, , ])
      input <- array(input, dim = c(1, dim(baseline_seq)[-1]))
      input <- tensorflow::tf$cast(input, baseline_seq$dtype)
      output <- compute_gradients(
        model = model,
        input_tensor = input,
        target_class_idx = target_class_idx,
        input_idx = NULL)
      output <- tensorflow::tf$squeeze(output)
      l[[j]] <- output
    }
    path_gradients <- tensorflow::tf$stack(l)
  }
  return(path_gradients)
}


#' Heatmap of integrated gradient scores
#' 
#' Creates a heatmap from output of \code{\link{integrated_gradients}} function. The first row contains 
#' the column-wise absolute sums of IG scores and the second row the sums. Rows 3 to 6 contain the IG scores for each 
#' position and each nucleotide. The last row contains nucleotide information.
#'
#' @param integrated_grads Matrix of integrated gradient scores (output of \code{\link{integrated_gradients}} function).
#' @param input_seq Input sequence for model. Should be the same as \code{input_seq} input for corresponding
#' \code{\link{integrated_gradients}} call that computed input for \code{integrated_grads} argument.
#' @examples 
#' model <- create_model_lstm_cnn(layer_lstm = 8, layer_dense = 3, maxlen = 20, verbose = FALSE)
#' random_seq <- sample(0:3, 20, replace = TRUE)
#' input_seq <- array(keras::to_categorical(random_seq), dim = c(1, 20, 4))
#' ig <- integrated_gradients(
#'   input_seq = input_seq,
#'   target_class_idx = 3,
#'   model = model)
#' heatmaps_integrated_grad(integrated_grads = ig,
#'                          input_seq = input_seq)
#' @export
heatmaps_integrated_grad <- function(integrated_grads,
                                     input_seq) {

  if (is.list(input_seq)) {
    for (i in 1:length(input_seq)) {
      input_seq[[i]] <- tensorflow::tf$cast(input_seq[[i]], dtype = "float32")
    }

    for (i in 1:length(integrated_grads)) {
      integrated_grads[[i]] <- tensorflow::tf$cast(integrated_grads[[i]], dtype = "float32")
    }
  } else {
    input_seq <- tensorflow::tf$cast(input_seq, dtype = "float32")
    integrated_grads <- tensorflow::tf$cast(integrated_grads, dtype = "float32")
  }


  if (is.list(input_seq)) {
    num_input <- length(input_seq)
    attribution_mask <- list()
    nuc_matrix <- list()
    nuc_seq <- list()
    sum_nuc <- list()
    for (i in 1:length(integrated_grads)) {
      py$integrated_grads <- integrated_grads[[i]]
      py_run_string("attribution_mask = tf.reduce_sum(tf.math.abs(integrated_grads), axis=-1)")
      py_run_string("sum_nuc = tf.reduce_sum(integrated_grads, axis=-1)")
      attribution_mask[[i]] <- py$attribution_mask
      attribution_mask[[i]] <- as.matrix(attribution_mask[[i]], nrow = 1) %>% as.data.frame()
      colnames(attribution_mask[[i]]) <- "abs_sum"
      sum_nuc[[i]] <- py$sum_nuc
      sum_nuc[[i]] <- as.matrix(sum_nuc[[i]], nrow = 1) %>% as.data.frame()
      colnames(sum_nuc[[i]]) <- "sum"

      if (length(dim(integrated_grads[[i]])) == 3) {
        nuc_matrix[[i]] <- as.matrix(integrated_grads[[i]][1, , ])
      }
      if (length(dim(integrated_grads[[i]])) == 2) {
        nuc_matrix[[i]] <- as.matrix(integrated_grads[[i]])
      }
      amb_nuc <- (apply(input_seq[[i]][1, ,], 1, max) %>% as.character()) != "1"
      nuc_seq[[i]] <- apply(input_seq[[i]][1, ,], 1, which.max) %>% as.character()
      nuc_seq[[i]] <- nuc_seq[[i]] %>% stringr::str_replace_all("1", "A") %>%
        stringr::str_replace_all("2", "C") %>%
        stringr::str_replace_all("3", "G") %>%
        stringr::str_replace_all("4", "T")
      nuc_seq[[i]][amb_nuc] <- "0"
      rownames(nuc_matrix[[i]]) <- nuc_seq[[i]]
      colnames(nuc_matrix[[i]]) <- c("A", "C", "G", "T")
    }

  } else {
    num_input <- 1
    py$integrated_grads <- integrated_grads
    py_run_string("attribution_mask = tf.reduce_sum(tf.math.abs(integrated_grads), axis=-1)")
    py_run_string("sum_nuc = tf.reduce_sum(integrated_grads, axis=-1)")
    #py_run_string("mean_nuc = tf.reduce_mean(integrated_grads, axis=-1)")

    attribution_mask <- py$attribution_mask
    attribution_mask <- as.matrix(attribution_mask, nrow = 1) %>% as.data.frame()
    colnames(attribution_mask) <- "abs_sum"

    sum_nuc <- py$sum_nuc
    sum_nuc <- as.matrix(sum_nuc, nrow = 1) %>% as.data.frame()
    colnames(sum_nuc) <- "sum"

    if (length(dim(integrated_grads)) == 3) {
      nuc_matrix <- as.matrix(integrated_grads[1, , ])
    }
    if (length(dim(integrated_grads)) == 2) {
      nuc_matrix <- as.matrix(integrated_grads)
    }
    amb_nuc <- (apply(input_seq[1, ,], 1, max) %>% as.character()) != "1"
    nuc_seq <- apply(input_seq[1, ,], 1, which.max) %>% as.character()
    nuc_seq <- nuc_seq %>% stringr::str_replace_all("1", "A") %>%
      stringr::str_replace_all("2", "C") %>%
      stringr::str_replace_all("3", "G") %>%
      stringr::str_replace_all("4", "T")
    nuc_seq[amb_nuc] <- "0"
    rownames(nuc_matrix) <- nuc_seq
    colnames(nuc_matrix) <- c("A", "C", "G", "T")
  }

  if (num_input == 1) {
    ig_min <- keras::k_min(integrated_grads)$numpy()
    ig_max <- keras::k_max(integrated_grads)$numpy()
    col_fun <- circlize::colorRamp2(c(ig_min, mean(c(ig_max, ig_min)) , ig_max), c("blue", "white", "red"))
  } else {
    col_fun <- list()
    for (i in 1:num_input) {
      ig_min <- keras::k_min(integrated_grads[[i]])$numpy()
      ig_max <- keras::k_max(integrated_grads[[i]])$numpy()
      col_fun[[i]] <- circlize::colorRamp2(c(ig_min, mean(c(ig_max, ig_min)) , ig_max), c("blue", "white", "red"))
    }
  }

  hm_list <- list()
  if (num_input == 1) {
    row_ha = ComplexHeatmap::columnAnnotation(abs_sum = attribution_mask[,1], sum = sum_nuc[,1]) # mean = mean_nuc[,1]
    if (length(unique(row.names(nuc_matrix))) == 4) {
      nuc_col <- c("A" = "red", "C" = "green", "G" = "blue", "T" = "yellow")
    }
    if (length(unique(row.names(nuc_matrix))) == 5) {
      nuc_col <- c("A" = "red", "C" = "green", "G" = "blue", "T" = "yellow", "0" = "white")
    }
    ha <- ComplexHeatmap::HeatmapAnnotation(nuc = row.names(nuc_matrix), col = list(nuc = nuc_col))
    hm_list[[1]] <- ComplexHeatmap::Heatmap(matrix = t(nuc_matrix),
                                            name = "hm",
                                            top_annotation = row_ha,
                                            bottom_annotation = ha,
                                            col = col_fun,
                                            cluster_rows = FALSE,
                                            cluster_columns = FALSE,
                                            column_names_rot = 0
    )
  } else {
    for (i in 1:num_input) {
      row_ha <- ComplexHeatmap::columnAnnotation(abs_sum = attribution_mask[[i]][,1], sum = sum_nuc[[i]][,1])
      if (length(unique(row.names(nuc_matrix[[i]]))) == 4) {
        nuc_col <- c("A" = "red", "C" = "green", "G" = "blue", "T" = "yellow")
      }
      if (length(unique(row.names(nuc_matrix[[i]]))) == 5) {
        nuc_col <- c("A" = "red", "C" = "green", "G" = "blue", "T" = "yellow", "0" = "white")
      }
      ha <- ComplexHeatmap::HeatmapAnnotation(nuc = row.names(nuc_matrix[[i]]), col = list(nuc = nuc_col))
      hm_list[[i]] <- ComplexHeatmap::Heatmap(matrix = t(nuc_matrix[[i]]),
                                              name = paste0("hm_", i),
                                              top_annotation = row_ha,
                                              bottom_annotation = ha,
                                              col = col_fun[[i]],
                                              cluster_rows = FALSE,
                                              cluster_columns = FALSE,
                                              column_names_rot = 0
      )
    }
  }
  hm_list
}
