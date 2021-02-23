test_that("mockSVGenes generates a SpatialExperiment", {
  output <- mockSVGenes(20, 2, 30)

  expect_s4_class(output, "SpatialExperiment")

  rd <- SummarizedExperiment::rowData(output)
  expect_equal(nrow(rd), 20)

  cd <- SummarizedExperiment::colData(output)
  expect_lt(abs(nrow(cd) - 30), 3)
})
