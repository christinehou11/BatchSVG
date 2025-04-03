#' @keywords internal
#' @title Identify Batch-Biased Features in Spatially Variable Genes
#' @description `BatchSVG` is a feature-based Quality Control (QC) to identify 
#' SVGs on spatial transcriptomics data with specific types of batch effect. 
#' Regarding to the spatial transcriptomics data experiments, the batch can be 
#' defined as "sample", "sex", and etc.The `BatchSVG` method is based on 
#' binomial deviance model (Townes et al, 2019) and applies cutoffs based on 
#' the number of standard deviation (nSD) of relative change in deviance and 
#' rank difference as the data-driven thresholding approach to detect the 
#' batch-biased outliers.
#' @docType package
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot scale_color_manual aes labs
#' @name BatchSVG-package
#' @aliases BatchSVG-package
#' @author Christine Hou
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL