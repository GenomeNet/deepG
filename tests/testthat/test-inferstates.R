context("infer-states")

test_that("Check Cell States", {
  
  # skip_if_no_keras()
  # 
  # sequence <- preprocessSemiRedundant(strrep("ATGC", 100), maxlen = 80, vocabulary =  c("l", "p", "a", "c", "g", "t"))
  # states <- getStates(model.path = "example_model_cpu_new_full_model.hdf5", x = sequence$X, maxlen = 80)
  # 
  # expect_error(getStates())
  # expect_error(getStates(model.path = ""))
  # expect_error(getStates(x = ""))
  # expect_error(getStates(model.path = "", x = ""))
  # 
  # expect_message(getStates(model.path = "example_model_cpu_new_full_model.hdf5", x = sequence$X, verbose = T, maxlen = 80))
  # expect_type(states, "double")
  
  #file.remove("output_states.csv")
})

test_that("Check Cell States from FASTA files", {
  
  skip_if_no_keras()
  
  expect_error(getStatesFromFasta())
  expect_error(getStatesFromFasta(""))
})
