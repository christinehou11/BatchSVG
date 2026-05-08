# List of data frames of nSD(s) on the batch of subject

The `list_batch_df` dataset contains results from
[`featureSelect()`](https://christinehou11.github.io/BatchSVG/reference/featureSelect.md)
applied to the spatial transcriptomics data from the `spatialLIBD`
package.

## Usage

``` r
data(list_batch_df)
```

## Format

A named list of data frames, where each element corresponds to a batch
effect:

- **"gene_id"**: Gene identifier.

- **"gene_name"**: Gene name.

- **"dev_default"**: Deviance score without batch correction.

- **"dev\_(batch name)"**: Deviance score with batch correction.

- **"rank_default"**: Rank of the gene based on deviance without batch
  correction.

- **"rank\_(batch name)"**: Rank of the gene based on deviance with
  batch correction.

- **"d_diff"**: Relative change in deviance between default and
  batch-corrected models.

- **"nSD_dev\_(batch name)"**: number of standard deviation of relative
  change in deviance for the batch.

- **"r_diff"**: Rank difference between default and batch-corrected
  models.

- **"nSD_rank\_(batch name)"**: number of standard deviation of rank
  difference for the batch.

## Source

<https://github.com/christinehou11/BatchSVG/blob/main/inst/scripts/make-list_batch_df.R>
