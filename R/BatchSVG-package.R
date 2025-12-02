#' @keywords internal
#' @title Identify Batch-Biased Features in Spatially Variable Genes
#' @description `BatchSVG` is a method to identify batch-biased spatially variable genes (SVGs) in spatial transcriptomics data. The batch variable can be defined as sample, donor sex, or other batch effects of interest. The BatchSVG method is based on the binomial deviance model (Townes et al, 2019).
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