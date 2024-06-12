#' Rds data generator
#' 
#' Creates training batches from rds files. Rds files must contain a  
#' list of length 2 (input/target) or of length 1 (for language model).
#' If \code{target_len} is not NULL will take the last \code{target_len} entries of 
#' the first list element as targets and the rest as input.    
#' 
#' @inheritParams generator_fasta_label_header_csv
#' @param rds_folder Path to input files.
#' @param target_len Number of target nucleotides for language model.
#' @param delete_used_files Whether to delete file once used. Only applies for rds files. 
#' @examples 
#' # create 3 rds files
#' rds_folder <- tempfile()
#' dir.create(rds_folder)
#' batch_size <- 7
#' maxlen <- 11
#' voc_len <- 4
#' for (i in 1:3) {
#'   x <- sample(0:(voc_len-1), maxlen*batch_size, replace = TRUE)
#'   x <- keras::to_categorical(x, num_classes = voc_len)
#'   x <- array(x, dim = c(batch_size, maxlen, voc_len))
#'   y <- sample(0:2, batch_size ,replace = TRUE)
#'   y <- keras::to_categorical(y, num_classes = 3)
#'   xy_list <- list(x, y)
#'   file_name <- paste0(rds_folder, "/file_", i, ".rds")
#'   saveRDS(xy_list, file_name) 
#' }
#' 
#' # create generator
#' gen <- generator_rds(rds_folder, batch_size = 2)
#' z <- gen()
#' x <- z[[1]]
#' y <- z[[2]]
#' x[1, , ]
#' y[1, ]
#' 
#' @export
generator_rds <- function(rds_folder, batch_size, path_file_log = NULL,
                          max_samples = NULL,
                          proportion_per_seq = NULL,
                          target_len = NULL,
                          seed = NULL, delete_used_files = FALSE,
                          reverse_complement = FALSE,
                          sample_by_file_size = FALSE,
                          n_gram = NULL, n_gram_stride = 1,
                          reverse_complement_encoding = FALSE,
                          add_noise = NULL,
                          reshape_xy = NULL) {
  
  if (!is.null(reshape_xy)) {
    reshape_xy_bool <- TRUE
    reshape_x_bool <- ifelse(is.null(reshape_xy$x), FALSE, TRUE)
    reshape_y_bool <- ifelse(is.null(reshape_xy$y), FALSE, TRUE)
    reshape_sw_bool <- ifelse(is.null(reshape_xy$sw), FALSE, TRUE)
  } else {
    reshape_xy_bool <- FALSE
  }
  
  if (!is.null(seed)) set.seed(seed)
  is_lm <- !is.null(target_len)
  
  if (!is.null(n_gram) & is_lm && (target_len < n_gram)) {
    stop("target_len needs to be at least as big as n_gram.")
  }
  
  rds_files <- list_fasta_files(rds_folder, format = "rds", file_filter = NULL)
  num_files <- length(rds_files)
  
  read_success <- FALSE
  while (!read_success) {
    tryCatch(
      expr = {
        if (sample_by_file_size) {
          file_prob <- file.info(rds_files)$size/sum(file.info(rds_files)$size)
          file_index <- sample(1:num_files, size = 1, prob = file_prob)
        } else {
          file_prob <- NULL
          rds_files <- sample(rds_files)
          file_index <- 1
        }
        rds_file <- readRDS(rds_files[file_index])
        read_success <- TRUE
      },
      error = function(e){ 
        
      }
    )
  }
  
  if (is.list(rds_file)) {
    x_complete <- rds_file[[1]]
    include_sw <- ifelse(length(rds_file) == 3, TRUE, FALSE) 
  } else {
    x_complete <- rds_file
    include_sw <- FALSE
  }
  
  if (!is_lm) y_complete <- rds_file[[2]]
  # TODO: adjust for different input format (input mix of 3D and 1D etc.)
  multi_input <- ifelse(is.list(x_complete), TRUE, FALSE)
  multi_output <- ifelse(length(rds_file) > 1 && is.list(rds_file[[2]]), TRUE, FALSE)
  
  if (multi_input) {
    x_dim_list <- list()
    size_splits_in <- list()
    for (i in 1:length(x_complete)) {
      x_dim_list[[i]] <- dim(x_complete[[i]])
      size_splits_in[[i]] <- x_dim_list[[i]][length(x_dim_list[[i]])] 
      if (i > 1) {
        if (length(x_dim_list[[i]]) != length(x_dim_list[[i-1]])) {
          stop("rds generator only works if separate inputs have same dimension size")
        }
      }
    }
    x_complete <- tensorflow::tf$concat(x_complete, 
                                        axis = as.integer(length(x_dim_list[[1]]) - 1)) %>% as.array()
  } 
  x_dim_start <- dim(x_complete)
  
  if (!is_lm) {
    if (multi_output) {
      y_dim_list <- list()
      size_splits_out <- list()
      for (i in 1:length(y_complete)) {
        y_dim_list[[i]] <- dim(y_complete[[i]])
        size_splits_out[[i]] <- y_dim_list[[i]][length(y_dim_list[[i]])] 
        if (i > 1) {
          if (length(y_dim_list[[i]]) != length(y_dim_list[[i-1]])) {
            stop("rds generator only works if separate outputs have same dimension size")
          }
        }
      }
      y_complete <- tensorflow::tf$concat(y_complete,
                                          axis = as.integer(length(y_dim_list[[1]]) - 1)) %>% as.array()
    } 
    y_dim_start <- dim(y_complete)
    if (is.null(y_dim_start)) y_dim_start <- length(y_complete)
    if (x_dim_start[1] != y_dim_start[1]) {
      stop("Different number of samples for input and target")
    }
  }
  
  if (include_sw) {
    sw_complete <- rds_file[[3]]
    multi_sw <- is.list(sw_complete)
    if (multi_sw) {
      # TODO:  have temporal/non-temporal in same batch 
      sw_temporal <- ifelse(dim(sw_complete[[1]])[2] == 1, FALSE, TRUE)
      sw_dim <- lapply(sw_complete,dim)
      
      sw_dim_list <- list()
      size_splits_sw <- list()
      for (i in 1:length(sw_complete)) {
        sw_dim_list[[i]] <- dim(sw_complete[[i]])
        size_splits_sw[[i]] <- sw_dim_list[[i]][length(sw_dim_list[[i]])] 
      }
      sw_complete <- tensorflow::tf$concat(sw_complete, 
                                           axis = as.integer(length(sw_dim_list[[1]]) - 1)) %>% as.array()
      
    } else {
      sw_temporal <- ifelse(dim(sw_complete)[2] == 1, FALSE, TRUE)
      sw_dim <- dim(sw_complete) 
    }
    sw_dim_start <- dim(sw_complete)
  }
  
  sample_index <- 1:x_dim_start[1]
  
  if (!is.null(proportion_per_seq)) {
    sample_index <- sample(sample_index, min(length(sample_index), length(sample_index) * proportion_per_seq))
  }
  
  if (!is.null(max_samples)) {
    sample_index <- sample(sample_index, min(length(sample_index), max_samples))
  }
  
  if (!is.null(path_file_log)) {
    utils::write.table(x = rds_files[1], file = path_file_log, row.names = FALSE, col.names = FALSE)
  }
  
  x_dim <- x_dim_start
  if (!is_lm) y_dim <- y_dim_start
  if (include_sw) sw_dim <- sw_dim_start
  
  function() {
    
    # TODO: adjust for multi input/output
    x_index <- 1
    x <- array(0, c(batch_size, x_dim[-1]))
    if (is_lm) {
      y <- vector("list", target_len)
    } else {
      y <- array(0, c(batch_size, y_dim[2]))
    }
    
    if (include_sw) {
      sw <- array(0, c(batch_size, sw_dim[2]))
    }
    
    while (x_index <= batch_size) {
      
      # go to next file if sample index empty
      if (length(sample_index) == 0) {
        if (num_files == 1) {
          sample_index <<- 1:x_dim[1]
        } else {
          
          if (delete_used_files) file.remove(rds_files[file_index])
          
          read_success <- FALSE
          while (!read_success) {
            tryCatch(
              expr = {
                if (sample_by_file_size) {
                  file_index <<- sample(1:num_files, size = 1, prob = file_prob)
                } else {
                  file_index <<- file_index + 1
                  if (file_index > num_files) {
                    file_index <<- 1
                    rds_files <<- sample(rds_files)
                  }
                }
                rds_file <<- readRDS(rds_files[file_index])
                read_success <- TRUE
              },
              error = function(e){ 
                
              }
            )
          }
          
          if (multi_input) {
            # combine inputs in one tensor
            x_complete <<- tensorflow::tf$concat(rds_file[[1]], 
                                                 axis = as.integer(length(x_dim_list[[1]]) - 1)) %>% as.array()
          } else {
            x_complete <<- rds_file[[1]]
          } 
          x_dim <<- dim(x_complete)
          
          if (include_sw) {
            if (multi_sw) {
              sw_complete <- tensorflow::tf$concat(rds_file[[3]], 
                                                   axis = as.integer(length(sw_dim_list[[1]]) - 1)) %>% as.array()
            } else {
              sw_complete <<- rds_file[[3]]
            }
          }
          
          if (!is_lm) {
            if (multi_output) {
              y_complete <- tensorflow::tf$concat(rds_file[[2]], 
                                                  axis = as.integer(length(y_dim_list[[1]]) - 1)) %>% as.array()
            } else {
              y_complete <<- rds_file[[2]]
            }
            y_dim <<- dim(y_complete)
            #if (is.null(y_dim)) y_dim <<- length(y_complete)
          }
          
          if (!is_lm && (x_dim[1] != y_dim[1])) {
            print(x_dim)
            print(y_dim)
            stop("Different number of samples for input and target")
          }
          
          sample_index <<- 1:x_dim[1]
          if (!is.null(proportion_per_seq)) {
            if (length(sample_index) > 1) {
              sample_index <<- sample(sample_index, min(length(sample_index), max(1, floor(length(sample_index) * proportion_per_seq))))
            }
          }
          
          if (!is.null(max_samples)) {
            if (length(sample_index) > 1) {
              sample_index <<- sample(sample_index, min(length(sample_index), max_samples))
            }
          }
        }
        
        # log file
        if (!is.null(path_file_log)) {
          utils::write.table(x = rds_files[file_index], file = path_file_log, append = TRUE, col.names = FALSE, row.names = FALSE)
        }
      }
      
      if (length(sample_index) == 0) next
      if (length(sample_index) == 1) {
        index <- sample_index
      } else {
        index <- sample(sample_index, min(batch_size - x_index + 1, length(sample_index)))
      }
      
      #subsetting
      subset_index <- x_index:(x_index + length(index) - 1)
      
      if (length(x_dim) == 4) {
        x[subset_index, , , ] <- x_complete[index, , , ]
      }
      
      if (length(x_dim) == 3) {
        x[subset_index, , ] <- x_complete[index, , ]
      }
      
      if (length(x_dim) == 2) {
        x[subset_index, ] <- x_complete[index, ]
      }
      
      if (!is_lm) {
        y[subset_index, ] <- y_complete[index, ]
      }
      
      if (include_sw) {
        sw[subset_index, ] <- sw_complete[index, ]
      }
      
      x_index <- x_index + length(index)
      sample_index <<- setdiff(sample_index, index)
      
    }
    
    if (is_lm) {
      for (m in 1:target_len) {
        if (batch_size == 1) {
          y[[m]] <- matrix(x[ , x_dim[2] - target_len + m, ], nrow = 1)
        } else {
          y[[m]] <- x[ , x_dim[2] - target_len + m, ]
        }
      }
      
      x <- x[ , 1:(x_dim[2] - target_len), ]
      if (batch_size == 1) {
        dim(x) <- c(1, dim(x))
      }
    }
    
    if (!is.null(n_gram) & is_lm) {
      y <- do.call(rbind, y)
      y_list <- list()
      for (i in 1:batch_size) {
        index <- (i-1)  + (1 + (0:(target_len-1)) * batch_size)
        n_gram_matrix <- n_gram_of_matrix(input_matrix = y[index, ], n = n_gram)
        y_list[[i]] <- n_gram_matrix
      }
      y_tensor <- keras::k_stack(y_list, axis = 1L) %>% keras::k_eval()
      y <- vector("list", dim(y_tensor)[2])
      for (i in 1:dim(y_tensor)[2]) {
        y_subset <- y_tensor[ , i, ]
        if (batch_size == 1) y_subset <- matrix(y_subset, nrow = 1)
        y[[i]] <- y_subset
      }
      
      if (is.list(y) & length(y) == 1) {
        y <- y[[1]]
      }
      
      if (n_gram_stride > 1 & is.list(y)) {
        stride_index <- 0:(length(y)-1) %% n_gram_stride == 0
        y <- y[stride_index]
        
      }
    }
    
    if (!is.null(add_noise)) {
      noise_args <- c(add_noise, list(x = x))
      x <- do.call(add_noise_tensor, noise_args)
    }
    
    if (reverse_complement_encoding){
      x_1 <- x
      x_2 <- array(x_1[ , (dim(x)[2]):1, 4:1], dim = dim(x))
      x <- list(x_1, x_2)
    }
    
    if (multi_input) {
      x <- tensorflow::tf$split(x, num_or_size_splits = size_splits_in, axis = as.integer(length(x_dim)-1))
    }
    
    if (multi_output) {
      y <- tensorflow::tf$split(y, num_or_size_splits = size_splits_out, axis = as.integer(length(y_dim)-1))
    }
    
    if (include_sw && multi_sw) {
      sw <- tensorflow::tf$split(sw, num_or_size_splits = size_splits_sw, axis = as.integer(length(sw_dim)-1))
    }
    
    if (reshape_xy_bool) {
      if (reshape_x_bool) x <- reshape_xy$x(x)
      if (reshape_y_bool) y <- reshape_xy$y(y)
      if (reshape_sw_bool) sw <- reshape_xy$sw(sw)
    }
    
    if (include_sw) {
      return(list(x, y, sw))
    } else {
      return(list(x, y))
    }
    
  }
}
