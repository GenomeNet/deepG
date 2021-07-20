context("gpu-management")

test_that("Check the functions in gpu-management", {
  expect_error(startGPUSession())
  expect_error(startGPUSession(""))
})
