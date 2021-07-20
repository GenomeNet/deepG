context("utlis")

test_that("Check functions in utlis.R", {
  
  expect_error(writeHdf5())
  expect_error(writeHdf5(""))
  
  expect_error(writeDict())
  expect_error(writeDict(""))
  
  expect_error(messagef())
  expect_error(calculateStepsPerEpoch())
})
