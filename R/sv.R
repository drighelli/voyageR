#' Identifies spatially variable genes.
#'
#' You can select one among the following statistical tests:
#' * spatialde, from [Svensson et al, 2018](https://www.nature.com/articles/nmeth.4636).
#' * spark, from [Sun et al, 2020](https://www.nature.com/articles/s41592-019-0701-7).
#'
#' @param x Spatial gene expression object.
#' @param method Name of the method for testing spatial variabily.
#' @param ... Additional arguments.
#'
#' @return Extended object with DE information.
#'
#' @importFrom methods setGeneric
#' @export
#' @docType methods
#' @rdname svGenes-methods
setGeneric("svGenes", function(x, method, ...) standardGeneric("svGenes"))

#' Identifies spatially variable gene from the given SpatialExperiment.
#'
#' @param x SpatialExpression object.
#' @param method Name of the method for testing spatial variability.
#' @param assay Extract gene expression from this assay.
#'
#' @return Extended SpatialExpression object.
#'
#' @importFrom methods setMethod
#' @importFrom SummarizedExperiment assay
#' @importFrom SpatialExperiment spatialCoords
#' @export
#'
#' @rdname svGenes-methods
#' @aliases svGenes,SpatialExperiment,ANY-method
#' 
#' @examples 
#' svGenes(mockSVGenes(100, 4, 40), method = "spatialde")
#' 
setMethod("svGenes",
          "SpatialExperiment",
          function(x,
                   method = c("spatialde", "spark"),
                   assay = "counts") {
            method <- match.arg(method)
            counts <- SummarizedExperiment::assay(x, assay)
            switch(method,
                   spatialde = spatialde_svg(x, counts),
                   spark = spark_svg(x, counts))
          })

spatialde_svg <- function(x, counts) {
  pacman::p_load_gh("sales-lab/spatialDE@wrap_functions")
  
  coordinates <- as.data.frame(spatialCoords(x))
  
  sample_info <- coordinates
  sample_info$total_counts <- colSums(counts)
  
  stabilized <- spatialDE::stabilize(counts)
  regressed <- spatialDE::regress_out(sample_info, stabilized)
  spatialDE::run(coordinates, regressed)
}

spark_svg <- function(x, counts) {
  
}
