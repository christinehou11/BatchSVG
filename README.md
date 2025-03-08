## BatchSVG

`BatchSVG` is a feature-based Quality Control (QC) to identify 
SVGs on spatial transcriptomics data with specific types of batch effect. 
Regarding to the spatial transcriptomics data experients, the batch can be 
defined as "sample", "sex", and etc.The `BatchSVG` method is based on 
binomial deviance model ([Townes et al, 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1861-6))) 
and applies cutoffs based on 
the number of standard deviation (nSD) of relative change in deviance and 
rank difference as the data-driven thresholding approach to detect the 
batch-biased outliers.

#### Installation

`BatchSVG` is a R package available in *Bioconductor* version 3.19 and later. 
You can install `BatchSVG` by using the following commands in R session from 
*Bioconductor*.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("BatchSVG")

## Check that you have a valid Bioconductor installation
BiocManager::valid()
```

The development version can be obtained by

``` r
BiocManager::install("christinehou11/BatchSVG")
```

### Code of Conduct
  
Please note that the BatchSVG project is released with 
a [Contributor Code of Conduct](https://christinehou11.github.io/BatchSVG/CODE_OF_CONDUCT.html). 
By contributing to this project, you agree to abide by its terms.
