#' Generate SpatialExperiment instance including spatially variable genes.
#'
#' @param tot_genes Total number of gene profiles to be generated.
#' @param de_genes Number of spatially variable genes.
#' @param grid_size Genes will be spatially arranged on a size x size grid.
#'
#' @return A [SpatialExperiment] object.
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom stats rnbinom runif
#' @export
#' 
#' @examples 
#' mockSVGenes(1000, 10, 100)
#' 
mockSVGenes <- function(tot_genes, de_genes, grid_size) {
  coordinates <- data.frame(
    x = rep(seq.int(grid_size), grid_size),
    y = rep(seq.int(grid_size), each=grid_size)
  )
  
  mu <- 2^runif(tot_genes, -1, 5)
  counts <- matrix(rnbinom(tot_genes*grid_size*grid_size, mu=mu, size=10),
                   nrow=tot_genes)
  m <- (grid_size/2)+1
  mask <- coordinates$x < m & coordinates$y < m
  counts[seq.int(de_genes), mask] <- counts[seq.int(de_genes), mask] + 20
  
  rownames(counts) <- paste0("gene", seq.int(tot_genes))
  colnames(counts) <- paste("spot", coordinates$x, coordinates$y, sep = "_")
  
  SpatialExperiment(assays = SimpleList(counts = counts),
                    spatialCoords = coordinates)
}
