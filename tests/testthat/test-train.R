context("training")

test_that("Sucessful training from a dummy model", {
   
   skip_if_no_keras()
   skip_on_travis()
   
   maxlen <- 30
   batch_size <- 10
   
   model <- create_model_lstm_cnn(
      maxlen = maxlen,
      layer_dense = 4,
      layer_lstm = 8,
      solver = "adam",
      vocabulary_size = 4,
      compile = TRUE)
   
   trainedNetwork <- train(reduce_lr_on_plateau = FALSE,
                           train_type = "lm",
                           path = "/home/rmreches/deepG/tests/testthat/fasta",
                           path_val = "/home/rmreches/deepG/tests/testthat/fasta",
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
   
   trainedNetwork <- train(reduce_lr_on_plateau = FALSE,
                           train_type = "label_folder",
                           path = rep("/home/rmreches/deepG/tests/testthat/fasta", 2),
                           path_val = rep("/home/rmreches/deepG/tests/testthat/fasta", 2),
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

})
