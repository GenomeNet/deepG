context("predict")

test_that("Prediction of next character", {
  
  skip_if_no_keras()
  
  example.model <-
    keras::load_model_hdf5("example_model_cpu_new_full_model.hdf5")
  sequence <- strrep("A", 100)
  vocabulary <- c("l", "p", "a", "c", "g", "t")
  
  expect_error(predictNextNucleotide())
  expect_error(predictNextNucleotide(sequence = ""))
  expect_error(predictNextNucleotide(model = ""))
  
  predicted_NextNucleotide <- predictNextNucleotide(sequence = sequence, model = example.model, 
                                                    vocabulary = vocabulary)
  
  expect_message(predictNextNucleotide(sequence = sequence, model = example.model, verbose = T,  vocabulary = vocabulary))
  expect_silent(predictNextNucleotide(sequence = sequence, model = example.model,  vocabulary = vocabulary))
  expect_s4_class(predicted_NextNucleotide, "prediction")
  expect_type(predicted_NextNucleotide@next_char, "character")
  expect_type(predicted_NextNucleotide@probability, "double")
  expect_type(predicted_NextNucleotide@index, "integer")
  expect_type(predicted_NextNucleotide@alternative_probability, "double")
  expect_type(predicted_NextNucleotide@solution, "character")
})

test_that("Prediction of replacement of n characters", {
  
  skip_if_no_keras()
  
  example.model <-
    keras::load_model_hdf5("example_model_cpu_new_full_model.hdf5")
  sequence <- strrep("A", 100)
  vocabulary <- c("l", "p", "a", "c", "g", "t")
  
  expect_error(replaceChar())
  expect_error(replaceChar(sequence = "", model = ""))
  expect_type(replaceChar(sequence = sequence, model = example.model,  vocabulary = vocabulary), "character")
  expect_equivalent(nchar(replaceChar(sequence = sequence, model = example.model, vocabulary = vocabulary)), 100)
  
  expect_silent(replaceChar(sequence = sequence, model = example.model, vocabulary = vocabulary))
})

test_that("Evaluation of a model on .fasta files", {
  
  skip_if_no_keras()
  
  expect_error(evaluateFasta())
  expect_error(evaluateFasta(""))
})

