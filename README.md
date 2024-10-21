## BiasDetect

`BiasDetect` is a feature-based Quality Control (QC) for bias genes identification on spatial transcriptomics and snRNA-seq data with some types of batch effect. The `BiasDetect` method is based on binomial deviance model ([Townes et al, 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1861-6)) and applies cutoffs based on the standard deviation (SD) of deviance and rank difference metrics as the data-driven thresholding approach to find the batch-biased features.

#### Installation

`BiasDetect` is a R package available in *Bioconductor* version 3.19 and later. You can install `BiasDetect` by using the following commands in R session from *Bioconductor*:

Install `BiasDetect`

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("BiasDetect")

## Check that you have a valid Bioconductor installation
BiocManager::valid()
```

Load required packages

``` r
library("BiasDetect")
library("SummarizedExperiment")
library("SpatialExperiment")
library("ggplot2")
library("PRECAST")
library("scater")
library("ggspavis")
```

Additionally, you can install development version from [GitHub](https://christinehou11.github.io/BiasDetect):

``` r
BiocManager::install("christinehou11/BiasDetect")
```

#### Tutorial

A detailed tutorial is available in [Find Biased Features](https://jac-thom.github.io/findBiasedFeatures/) written by Jacqui Thompson.

### Code of Conduct
  
Please note that the BiasDetect project is released with a [Contributor Code of Conduct](https://christinehou11.github.io/BiasDetect/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
