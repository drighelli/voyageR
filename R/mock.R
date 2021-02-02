#' Generate SpatialExperiment instance including spatially variable genes.
#'
#' @param tot_genes Total number of gene profiles to be generated.
#' @param de_genes Number of spatially variable genes.
#' @param lambda Expected number of points over which genes are sampled.
#'
#' @return A [SpatialExperiment] object.
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom S4Vectors SimpleList DataFrame
#' @export
#' 
#' @examples 
#' mockSVGenes(20, 2, 20)
#' 
mockSVGenes <- function(tot_genes, de_genes, lambda) {
  pp <- trendsceek::sim_pois(lambda)
  
  low_expr <- rep(10, de_genes)
  high_expr <- rep(50, de_genes)
  pp <- trendsceek::add_markdist_step(pp, low_expr, high_expr)
  
  coordinates <- data.frame(x = pp$x, y = pp$y)

  counts <- rbind(
    t(pp$marks),
    matrix(stats::rpois((tot_genes - de_genes) * nrow(coordinates), 10),
           ncol = nrow(coordinates))
  )
  
  genes <- paste0("gene", seq.int(tot_genes))
  rownames(counts) <- genes
  colnames(counts) <- paste0("spot", seq.int(nrow(coordinates)))
  
  rowData <- DataFrame(gene = genes)
  
  SpatialExperiment(assays = SimpleList(counts = counts),
                    spatialCoords = coordinates,
                    rowData = rowData)
}
