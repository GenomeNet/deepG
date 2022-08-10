context("preprocess")

test_that("Check preprocessing", {
  
  z <- seq_encoding_lm(sequence = c(1,0,5,1,3,4,3,1,4,1,2),
                       maxlen = 5,
                       vocabulary = c("a", "c", "g", "t"),
                       start_ind = c(1,3),
                       ambiguous_nuc = "equal",
                       target_len = 1,
                       output_format = "target_right")
  
  x <- z[[1]]
  y <- z[[2]]
  
  expect_equivalent(x[1,1,], c(1,0,0,0))
  expect_equivalent(x[1,2,], c(0,0,0,0))
  expect_equivalent(x[1,3,], rep(0.25, 4))
  expect_equivalent(x[1,4,], c(1,0,0,0))
  expect_equivalent(x[1,5,], c(0,0,1,0))
  expect_equivalent(y[1,], c(0,0,0,1))
  
  expect_equivalent(x[2,1,], rep(0.25, 4))
  expect_equivalent(x[2,2,], c(1,0,0,0))
  expect_equivalent(x[2,3,], c(0,0,1,0))
  expect_equivalent(x[2,4,], c(0,0,0,1))
  expect_equivalent(x[2,5,], c(0,0,1,0))
  expect_equivalent(y[2,], c(1,0,0,0))
  
  # use character string as input
  z <- seq_encoding_lm(sequence = NULL,
                       maxlen = 5,
                       vocabulary = c("a", "c", "g", "t"),
                       start_ind = c(1,3),
                       ambiguous_nuc = "zero",
                       target_len = 1,
                       output_format = "target_right",
                       char_sequence = "ACTaaTNTNaZ")
  
  
  x <- z[[1]]
  y <- z[[2]]
  
  expect_equivalent(apply(x[1,,], 1, which.max), c(1,2,4,1,1))
  expect_equivalent(apply(x[2,,], 1, which.max), c(4,1,1,4,1))
  expect_equivalent(y[1,], c(0,0,0,1))
  expect_equivalent(y[2,], c(0,0,0,1))
  
  x <- seq_encoding_label(sequence = c(1,0,5,1,3,4,3,1,4,1,2),
                          maxlen = 5,
                          vocabulary = c("a", "c", "g", "t"),
                          start_ind = c(1,3),
                          ambiguous_nuc = "equal")
  
  expect_equivalent(x[1,1,], c(1,0,0,0))
  expect_equivalent(x[1,2,], c(0,0,0,0))
  expect_equivalent(x[1,3,], rep(0.25, 4))
  expect_equivalent(x[1,4,], c(1,0,0,0))
  expect_equivalent(x[1,5,], c(0,0,1,0))
  
  expect_equivalent(x[2,1,], rep(0.25, 4))
  expect_equivalent(x[2,2,], c(1,0,0,0))
  expect_equivalent(x[2,3,], c(0,0,1,0))
  expect_equivalent(x[2,4,], c(0,0,0,1))
  expect_equivalent(x[2,5,], c(0,0,1,0))
  
  # use character string as input
  x <- seq_encoding_label(maxlen = 5,
                          vocabulary = c("a", "c", "g", "t"),
                          start_ind = c(1,3),
                          ambiguous_nuc = "equal",
                          char_sequence = "ACTaaTNTNaZ")
  
  expect_equivalent(apply(x[1,,], 1, which.max), c(1,2,4,1,1))
  expect_equivalent(apply(x[2,,], 1, which.max), c(4,1,1,4,1))
  expect_equivalent(x[2,5,], rep(0.25, 4))
  
})
