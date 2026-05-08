# SVGs Plots

This function generates visualizations to assess the impact of batch
effects on spatially variable genes (SVGs) by analyzing changes in
deviance and rank. The function bins the deviations into normalized
standard deviation (nSD) intervals and creates histograms and scatter
plots to illustrate the distribution of batch effects.

## Usage

``` r
svg_nSD(list_batch_df, sd_interval_dev, sd_interval_rank)
```

## Arguments

- list_batch_df:

  A named list of data frames, where each data frame corresponds to a
  batch effect and contains columns with deviance and rank differences.

- sd_interval_dev:

  A numeric vector specifying the interval widths for standard deviation
  bins for each batch when analyzing the relative change in deviance.
  The order of values must correspond to the order of batches in
  `list_batch_df` . If a single value is provided, it is applied to all
  batches; otherwise, it must have the same length as `list_batch_df`.

- sd_interval_rank:

  `vector`: A numeric vector specifying the interval widths for standard
  deviation bins when analyzing rank differences. The order of values
  must correspond to the order of batches in `list_batch_df`. If a
  single value is provided, it is applied to all batches; otherwise, it
  must have the same length as `list_batch_df`.

## Value

A combined `ggplot` object containing:

- **Deviance Plots**:

  - Histogram of deviance differences across SVGs, colored by nSD
    intervals.

  - Scatter plot comparing deviance values before and after batch
    correction.

- **Rank Plots**:

  - Histogram of rank differences across SVGs, colored by nSD intervals.

  - Scatter plot comparing ranks before and after batch correction.

The function arranges plots for each batch in a grid format for easy
comparison.

## Examples

``` r
# use the result generated from featureSelect()
data(list_batch_df)
plots <- svg_nSD(list_batch_df = list_batch_df, 
   sd_interval_dev = 3, sd_interval_rank = 3)
```
