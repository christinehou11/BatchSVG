#' @rdname svg_nSD
#'
#' @name svg_nSD
#'
#' @title SVGs Plots with Relative Change in Deviance and Rank Difference
#'
#' @description Function to the spatially variable genes in relative change in
#'    deviance and rank difference colored by the number of standard deviation.
#'
#' @importFrom dplyr left_join filter
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom stats sd
#' @importFrom ggplot2
#'
#' @param batch_df \code{data.frame} : Input data frame generated from
#'    `featureSelection()` function using \code{SpatialExperiment} 
#'    object containing the raw data.
#'    
#' @param nSD_dev \code{integer}: Number of standard deviation (nSD) on relative
#'    change in deviance. The default value is 5.
#'
#' @param nSD_rank \code{integer}: Number of standard deviation (nSD) on rank
#'     difference. The default value is 5.
#'
#' @return The output values are returned as
#'     a list of data frames containing the deviance and rank with and without 
#'     the bacth effect, the corresponding difference, the corresponding nSD, 
#'     gene, gene name, and whether the gene is outlier defined by the chosen 
#'     deviance and rank cutoff. The length of the list of data frames is equal
#'     to the number of batch effects included.
#'
#' @export
#'
#' @examples
#' sample_id_df <- utils::read.csv(
#'     system.file("extdata","sample_id_df.csv",package = "BatchSVG"),
#'     row.names = 1, check.names = FALSE)
#' 
#' plots <- svg_nSD(sample_id_df, nSD_dev = 5, nSD_rank)
#' 

biasDetect <- function(batch_df, nSD_dev = 5, nSD_rank = 5) {
    
    stopifnot(
        
    )
    
    plots
    }