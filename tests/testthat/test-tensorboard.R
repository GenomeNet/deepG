context("tensorboard")

test_that("Tensorboard", {
  
  expect_error(ecoliEvaluation())
  expect_error(ecoliEvaluation(""))
})
