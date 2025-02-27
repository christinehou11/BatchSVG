#' Example SpatialExperiment Dataset for Batch Effect Analysis
#'
#' `spe_sub4` is an example `SpatialExperiment` object used to demonstrate 
#' batch effect correction and spatially variable gene (SVG) selection in the 
#' `BatchSVG` package. It takes the four samples ("V11L05-333_B1",
#' "V11L05-333_D1","V11L05-335_D1","V11L05-336_A1") from the full spatial 
#' transcriptomics data stored in the `humanHippocampus2024` data package.
#'
#' @format A `SpatialExperiment` object with:
#' \describe{
#'    \item{Dimensions}{2082 genes (rows) × 18,945 spots (columns).}
#'    \item{Assays}{\code{counts} (raw gene expression) and \code{logcounts} 
#'        (log-transformed normalized counts).}
#'    \item{Row metadata (rowData)}{Includes gene annotations such as 
#'        \code{source}, \code{type}, and \code{gene_search}.}
#'    \item{Column metadata (colData)}{Includes 150 variables such as 
#'        \code{sample_id}, \code{in_tissue}, and \code{nmf100}.}
#'    \item{Reduced dimensions}{Contains PCA, t-SNE, and UMAP embeddings 
#'        labeled as \code{10x_pca}, \code{10x_tsne}, and \code{10x_umap}.}
#'    \item{Spatial coordinates}{Available in \code{spatialCoords} with pixel 
#'        row and column positions.}
#'    \item{Image metadata}{Includes fields such as \code{sample_id}, 
#'        \code{image_id}, and \code{scaleFactor}.}}
#'        
#' @usage data(svgs_sub4)
#'
#' @source Processed spatial transcriptomics data for demonstration purposes.
#'
#' @examples
#' data(svgs_sub4)
"spe_sub4"

#' Filtered Spatially Variable Genes (SVGs) Subset
#'
#' `svgs_sub4` is a subset of spatially variable genes (SVGs) derived from the 
#' `res_ranks` dataset. It includes genes that were ranked within the top 2,000 
#' in at least two of the selected samples from the original dataset.
#'
#' @format A data frame with:
#' \describe{
#'   \item{gene_id}{Character. Ensembl Gene ID representing each gene.}
#'   \item{n}{Integer. The number of selected samples in which the gene ranked 
#'   within the top 2,000 spatially variable genes.}
#' }
#'
#' @source Filtered subset from `res_ranks` (ranked `nnSVG()` results).
#' 
#' @usage data(svgs_sub4)
#' 
#' @examples
#' data(svgs_sub4)
"svgs_sub4"