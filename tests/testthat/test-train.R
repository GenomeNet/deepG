context("training")

test_that("Sucessful training from a dummy model", {
  
 skip_if_no_keras()

 data("parenthesis")
 maxlen <- 30
 batch.size <- 10
 preprocessed <- preprocessSemiRedundant(substr(parenthesis, 1, 100),
                                                maxlen = maxlen)
 
 expect_error(trainNetwork(""))
 expect_error(trainNetwork(dataset = preprocessed, maxlen = 0))
 expect_error(trainNetwork(dataset = preprocessed, maxlen = ""))
 expect_error(trainNetwork(dataset = preprocessed, dropout.rate = ""))
 expect_error(trainNetwork(dataset = preprocessed, dropout.rate = 0))
 expect_error(trainNetwork(dataset = preprocessed, dropout.rate = 1))
 expect_error(trainNetwork(dataset = preprocessed, layer.size = ""))
 expect_error(trainNetwork(dataset = preprocessed, layer.size = 1))
 expect_error(trainNetwork(dataset = preprocessed, layer_lstm = ""))
 expect_error(trainNetwork(dataset = preprocessed, layer_lstm  = 1))
 expect_error(trainNetwork(dataset = preprocessed, batch.size = ""))
 expect_error(trainNetwork(dataset = preprocessed, batch.size  = 1))
 expect_error(trainNetwork(dataset = "", path = ""))

 skip_on_travis()
 
 model <- create_model_lstm_cnn(
   maxlen = maxlen,
   layer.size = 2,
   layers.lstm = 2,
   solver = "adam",
   use.cnn = FALSE,
   num_targets = 7,
   vocabulary.size = 7,
   compile = TRUE)
 
 trainedNetwork <- trainNetwork(dataset = preprocessed,
                                path = "",
                                model = model,
                                batch.size = batch.size,
                                epochs = 1,
                                tensorboard.log = "",
                                output = list(none = TRUE, # no output 
                                              checkpoints = TRUE, 
                                              tensorboard = TRUE,
                                              log = TRUE,
                                              serialize_model = TRUE,
                                              full_model = TRUE
                                ))

 expect_type(trainedNetwork, "list")
 expect_equal(length(trainedNetwork),2)
 expect_type(trainedNetwork[1], "list")
 expect_equal(length(trainedNetwork[[1]]),7)
 expect_type(trainedNetwork[2], "list")
 expect_equal(length(trainedNetwork[[2]]),4)

 expect_type(trainedNetwork[[1]][["batch_size"]],"integer")
 expect_equal(trainedNetwork[[1]][["batch_size"]],10)
 expect_type(trainedNetwork[[1]][["epochs"]],"integer")
 expect_equal(trainedNetwork[[1]][["epochs"]],1)
 expect_type(trainedNetwork[[1]][["metrics"]],"character")
 expect_equal(trainedNetwork[[1]][["metrics"]],c("loss","acc", "val_loss", "val_acc"))

 expect_type(trainedNetwork[[2]][["loss"]],"double")
 expect_type(trainedNetwork[[2]][["val_loss"]],"double")
})
