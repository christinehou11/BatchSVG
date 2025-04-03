#' List of data frames of nSD(s) on the batch of subject
#' 
#' The `list_batch_df` dataset contains results from `featureSelect()` applied 
#' to the spatial transcriptomics data from the `spatialLIBD` package.
#' 
#' @format A named list of data frames, where each element corresponds to a 
#'    batch effect:
#'    \itemize{
#'    \item \strong{"gene_id"}: Gene identifier.
#'    \item \strong{"gene_name"}: Gene name.
#'    \item \strong{"dev_default"}: Deviance score without batch correction.
#'    \item \strong{"dev_(batch name)"}: Deviance score with batch correction.
#'    \item \strong{"rank_default"}: Rank of the gene based on deviance 
#'        without batch correction.
#'    \item \strong{"rank_(batch name)"}: Rank of the gene based on deviance 
#'        with batch correction.
#'    \item \strong{"d_diff"}: Relative change in deviance between default 
#'        and batch-corrected models.
#'    \item \strong{"nSD_dev_(batch name)"}: number of standard deviation of  
#'        relative change in deviance for the batch.
#'    \item \strong{"r_diff"}: Rank difference between default and 
#'        batch-corrected models.
#'    \item \strong{"nSD_rank_(batch name)"}: number of standard deviation of 
#'        rank difference for the batch.
#'    }
#'    
#' @usage data(list_batch_df)
#' @source 
#' \url{https://github.com/christinehou11/BatchSVG/blob/main/inst/scripts/make-list_batch_df.R}
"list_batch_df"