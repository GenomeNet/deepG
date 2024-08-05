context("training")

test_that("Sucessful training from a dummy model", {
   
   testthat::skip_if_not_installed("tensorflow")
   
   # language model
   maxlen <- 30
   batch_size <- 10

   model <- create_model_lstm_cnn(
      maxlen = maxlen,
      layer_dense = 4,
      layer_lstm = 8,
      solver = "adam",
      vocabulary_size = 4,
      compile = TRUE)

   trainedNetwork <- train_model(reduce_lr_on_plateau = FALSE,
                                 train_type = "lm",
                                 path = "fasta",
                                 path_val = "fasta",
                                 model = model,
                                 train_val_ratio = 0.2,
                                 steps_per_epoch = 3,
                                 batch_size = batch_size,
                                 epochs = 1)

   expect_type(trainedNetwork, "list")
   expect_type(trainedNetwork[["metrics"]][["loss"]], "double")
   expect_gte(trainedNetwork[["metrics"]][["loss"]], 0)
   expect_type(trainedNetwork[["metrics"]][["val_loss"]], "double")
   expect_gte(trainedNetwork[["metrics"]][["val_loss"]], 0)

   # label folder
   model <- create_model_lstm_cnn(
      maxlen = maxlen,
      kernel_size = c(4,4),
      pool_size = c(2,2),
      filters = c(2,4),
      layer_dense = 2,
      layer_lstm = 8,
      solver = "adam",
      vocabulary_size = 4,
      compile = TRUE)

   trainedNetwork <- train_model(reduce_lr_on_plateau = FALSE,
                                 train_type = "label_folder",
                                 path = rep("fasta", 2),
                                 path_val = rep("fasta", 2),
                                 model = model,
                                 vocabulary_label = c("A", "B"),
                                 train_val_ratio = 0.2,
                                 steps_per_epoch = 3,
                                 batch_size = batch_size,
                                 epochs = 1)

   expect_type(trainedNetwork, "list")
   expect_type(trainedNetwork[["metrics"]][["loss"]], "double")
   expect_gte(trainedNetwork[["metrics"]][["loss"]], 0)
   expect_type(trainedNetwork[["metrics"]][["val_loss"]], "double")
   expect_gte(trainedNetwork[["metrics"]][["val_loss"]], 0)
   
   # train/val split with csv
   maxlen <- 3
   model <- create_model_lstm_cnn(
      maxlen = maxlen,
      layer_dense = 2,
      layer_lstm = 6,
      solver = "adam",
      vocabulary_size = 4,
      compile = TRUE)
   
   files <- as.list(list.files(c("fasta_2", "fasta"), full.names = TRUE))
   
   # target csv
   A <- sample(c(0,1), length(files), replace = TRUE)
   B <- ifelse(A == 0, 1, 0)
   target_df <- data.frame(file = basename(unlist(files)), A = A, B = B)
   target_from_csv <- tempfile(fileext = ".csv")
   write.csv(target_df, target_from_csv, row.names = FALSE)
   
   # train/val csv
   train_val_split_csv <- tempfile(fileext = ".csv")
   ttv_df <- data.frame(file = basename(unlist(files)), 
                        type = rep(c("train", "validation"), each = 2))
   train_files <- ttv_df$file[ttv_df$type == "train"]
   val_files <- ttv_df$file[ttv_df$type == "validation"]
   write.csv(ttv_df, train_val_split_csv, row.names = FALSE)
   
   path_file_log <- tempfile(fileext = ".csv")
   
   trainedNetwork <- train_model(reduce_lr_on_plateau = FALSE,
                                 train_type = "label_csv",
                                 path = files,
                                 path_val = files,
                                 model = model,
                                 vocabulary_label = c("A", "B"),
                                 train_val_ratio = 0.2,
                                 steps_per_epoch = 3,
                                 batch_size = 4,
                                 max_samples = 2,
                                 target_from_csv = target_from_csv,
                                 train_val_split_csv = train_val_split_csv,
                                 path_file_log = path_file_log,
                                 epochs = 2)
   
   file_log <- read.csv(path_file_log)
   train_file_log <- unique(basename(file_log[,1]))
   expect_true(all(sort(train_file_log) == sort(train_files)))
   expect_true(all(sort(train_file_log) != sort(val_files)))
   

})
