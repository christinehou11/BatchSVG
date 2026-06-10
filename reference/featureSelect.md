# Feature selection

The function computes batch-adjusted deviance values, ranks the genes
accordingly, and quantifies batch effects in terms of standard
deviations from the mean difference. The list follows the order of batch
effects provided in `batch_effects`.

## Usage

``` r
featureSelect(input, batch_effects = NULL, VGs = NULL, verbose = TRUE)
```

## Arguments

- input:

  A `SpatialExperiment` object containing spatial transcriptomics data..

- batch_effects:

  A character vector specifying column names in `colData(input)` that
  indicate batch effects. Must match existing column names.

- VGs:

  A character vector specifying the variable genes (VGs) to be analyzed.
  Only genes present in this vector will be retained for feature
  selection.

- verbose:

  Logical (TRUE or FALSE). Default is TRUE. If TRUE, progress messages
  will be printed; If FALSE, messages will be suppressed.

## Value

A named list where each element corresponds to a batch effect. Each
batch contains a data frame with the following columns:

- **"gene_id"**: Gene identifier.

- **"gene_name"**: Gene name.

- **"dev_default"**: Deviance score without batch correction.

- **"dev\_"**: Deviance score with batch correction.

- **"rank_default"**: Rank of the gene based on deviance without batch
  correction.

- **"rank\_"**: Rank of the gene based on deviance with batch
  correction.

- **"d_diff"**: Relative change in deviance between default and
  batch-corrected models.

- **"nSD_dev\_"**: number of standard deviation of relative change in
  deviance for the batch.

- **"r_diff"**: Rank difference between default and batch-corrected
  models.

- **"nSD_rank\_"**: number of standard deviation of rank difference for
  the batch.

## Examples

``` r
library(spatialLIBD)
#> Loading required package: SpatialExperiment
#> Loading required package: SingleCellExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
spatialLIBD_spe <- fetch_data(type = "spe")
#> 
#> adding rname 'https://www.dropbox.com/s/f4wcvtdq428y73p/Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata?dl=1'
#> 
#> 2026-06-10 21:13:05.914904 loading file /home/runner/.cache/R/BiocFileCache/37b92a53361a_Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata%3Fdl%3D1
libd_svg <- read.csv(
    system.file("extdata","libd-all_nnSVG_p-05-features-df.csv",
              package = "BatchSVG"),
    row.names = 1, check.names = FALSE)
   
list_batch_df <- featureSelect(input = spatialLIBD_spe, 
   batch_effects = "subject", VGs = libd_svg$gene_id)
#> Running feature selection with batch...
#> Batch Effect: subject
#> Running feature selection without batch...
#> Calculating deviance and rank difference...
```
