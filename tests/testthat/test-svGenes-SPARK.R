
test_that("SPARK identifies variable genes", {
  set.seed(42)
  
  se <- mockSVGenes(100, 2, 100)
  output <- svGenes(se, method = "spark")
  
  expect_s4_class(output, "SpatialExperiment")
  
  rd <- rowData(output)
  expect_named(rd, c("gene", "spark"), ignore.order = TRUE)
  
  sde <- rd$spark
  expect_s3_class(sde, "data.frame")
  expect_true("adjusted_pvalue" %in% colnames(sde))
  expect_identical(rownames(sde), rd$gene)
  
  expect_true(all(sde[1:2, "qval"] < 0.05))
  expect_true(all(sde[3:20, "qval"] > 0.05))
})
