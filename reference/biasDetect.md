# Biased Genes Identification

Function to identify the bias genes based on user-selected threshold of
number of standard deviation in relative change in deviance and rank
difference.

## Usage

``` r
biasDetect(
  list_batch_df,
  threshold = c("both", "dev", "rank"),
  nSD_dev = NULL,
  nSD_rank = NULL,
  plot_point_size = 3,
  plot_point_shape = 16,
  plot_text_size = 3,
  plot_palette = "YlOrRd"
)
```

## Arguments

- list_batch_df:

  `list` : The list of data frame(s) generated from `featureSelection()`
  function. The length of the data frame list should be at least one.

- threshold:

  A character string specifying the filtering criterion. Must be one of:

  - `"dev"`: Filters genes based on the deviance threshold only.

  - `"rank"`: Filters genes based on the rank threshold only.

  - `"both"`: Filters genes based on either the deviance or rank
    threshold. Default is "both".

- nSD_dev:

  `integer`: A numeric vector specifying the number of standard
  deviation (nSD) for each batch when analyzing the relative change in
  deviance. The order of values must correspond to the order of batches
  in `list_batch_df`. Required if `threshold` is "dev" or "both". If a
  single value is provided, it is applied to all batches; otherwise, it
  must have the same length as `list_batch_df`.

- nSD_rank:

  `vector`: A numeric vector specifying the number of standard deviation
  (nSD) for each batch when analyzing rank differences. The order of
  values must correspond to the order of batches in `list_batch_df`.
  Required if `threshold` is "rank" or "both". If a single value is
  provided, it is applied to all batches; otherwise, it must have the
  same length as `list_batch_df`.

- plot_point_size:

  `vector`: A numeric vector specifying point sizes in plots. If asingle
  value is provided, it is applied to all batches.

- plot_point_shape:

  `vector`: A numeric vector specifying point shapes in plots. If a
  single value is provided, it is applied to all batches.

- plot_text_size:

  `vector`: A numeric vector specifying text label size in plots.
  Default is `3`.

- plot_palette:

  `vector`: A character string vector specifying the color palette for
  plots. Default is `"YlOrRd"`.

## Value

A named list where each element corresponds to a batch and contains:

- `"Plot"`: A diagnostic plot (either deviance, rank, or both).

- `"Table"`: A filtered data frame containing outlier genes based on the
  specified threshold.

## Examples

``` r
# use the result generated from featureSelect()
data(list_batch_df)
biaGenes <- biasDetect(list_batch_df = list_batch_df, threshold = "both", 
   nSD_dev = 3, nSD_rank = 3)
```
