context("predict")

test_that("Prediction of next character", {
  
  skip_if_no_keras()
  
  example.model <-
    keras::load_model_hdf5("example_model_cpu_new_full_model.hdf5")
  sequence <- strrep("A", 100)
  vocabulary <- c("l", "p", "a", "c", "g", "t")
 
})

test_that("Prediction of replacement of n characters", {
  
  skip_if_no_keras()
  
  example.model <-
    keras::load_model_hdf5("example_model_cpu_new_full_model.hdf5")
  sequence <- strrep("A", 100)
  vocabulary <- c("l", "p", "a", "c", "g", "t")
})

test_that("Evaluation of a model on .fasta files", {
  
  skip_if_no_keras()
})

