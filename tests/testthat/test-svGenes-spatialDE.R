test_that("spatialDE identifies variable genes", {
  set.seed(42)
  
  se <- mockSVGenes(100, 2, 100)
  output <- svGenes(se, method = "spatialde")
  
  expect_s4_class(output, "SpatialExperiment")
  
  rd <- rowData(output)
  expect_named(rd, c("gene", "spatialde"), ignore.order = TRUE)
  
  sde <- rd$spatialde
  expect_s3_class(sde, "data.frame")
  expect_true(all(c("g", "qval") %in% colnames(sde)))
  expect_identical(sde$g, rd$gene)

  expect_true(all(sde[1:2, "qval"] < 0.05))
  expect_true(all(sde[-c(1:2), "qval"] > 0.05))
})
