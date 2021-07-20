context("tSNE")

test_that("Check tSNE functions", {
  
  expect_error(generateStatesFromFolder())
  expect_error(generateStatesFromFolder(""))

  expect_error(extractCellFromStates())
  expect_error(extractCellFromStates(""))
  
  #expect_error(plotTsne())
  #expect_error(plotTsne(""))
  
  #file.remove("tsne.pdf")
})
