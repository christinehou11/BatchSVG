#' @rdname featureSelect
#'
#' @name featureSelect
#'
#' @title Feature selection
#'
#' @description The function computes batch-adjusted deviance values, ranks the 
#'    genes accordingly, and quantifies batch effects in terms of standard 
#'    deviations from the mean difference. The list follows the order of batch 
#'    effects provided in \code{batch_effects}.
#'
#' @importFrom scry devianceFeatureSelection
#' @importFrom dplyr left_join filter
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom stats sd
#'
#' @param input A \code{SpatialExperiment} object containing spatial 
#'    transcriptomics data..
#' @param batch_effects A character vector specifying column names in 
#'    \code{colData(input)} that indicate batch effects. Must match existing 
#'    column names.
#' @param VGs A character vector specifying the variable genes (VGs) to be 
#'     analyzed. Only genes present in this vector will be retained for 
#'     feature selection.
#' @param verbose Logical (TRUE or FALSE). Default is TRUE.
#'    If TRUE, progress messages will be printed; 
#'    If FALSE, messages will be suppressed.
#'
#' @return A named list where each element corresponds to a batch effect. 
#'    Each batch contains a data frame with the following columns:
#'    \itemize{
#'    \item \strong{"gene_id"}: Gene identifier.
#'    \item \strong{"gene_name"}: Gene name.
#'    \item \strong{"dev_default"}: Deviance score without batch correction.
#'    \item \strong{"dev_<batch>"}: Deviance score with batch correction.
#'    \item \strong{"rank_default"}: Rank of the gene based on deviance 
#'        without batch correction.
#'    \item \strong{"rank_<batch>"}: Rank of the gene based on deviance with 
#'        batch correction.
#'    \item \strong{"d_diff"}: Relative change in deviance between default 
#'        and batch-corrected models.
#'    \item \strong{"nSD_dev_<batch>"}: number of standard deviation of  
#'        relative change in deviance for the batch.
#'    \item \strong{"r_diff"}: Rank difference between default and 
#'        batch-corrected models.
#'    \item \strong{"nSD_rank_<batch>"}: number of standard deviation of 
#'        rank difference for the batch.
#'    }
#'
#' @export
#'
#' @examples
#' library(spatialLIBD)
#' spatialLIBD_spe <- fetch_data(type = "spe")
#' libd_svg <- read.csv(
#'     system.file("extdata","libd-all_nnSVG_p-05-features-df.csv",
#'               package = "BatchSVG"),
#'     row.names = 1, check.names = FALSE)
#'    
#' list_batch_df <- featureSelect(input = spatialLIBD_spe, 
#'    batch_effects = "subject", VGs = libd_svg$gene_id)
#' 
featureSelect <- function(input, batch_effects = NULL, VGs = NULL,
                            verbose = TRUE) {
    
    stopifnot(
        inherits(input, c("SpatialExperiment")), 
        is.character(VGs)
    )
    
    if (!is.null(batch_effects)) {
        stopifnot(is.character(batch_effects))
        for (batch in batch_effects) {
            if (!batch %in% names(colData(input))) {
                stop("The batch_effect is not a valid column")
            }
        }
    } else {
        stop("Please provide a valid batch_effect.")
    }
    
    input <- input[rowData(input)$gene_id %in% VGs, ]
    rownames(input) <- rowData(input)$gene_id
    
    if (verbose) message("Running feature selection with batch...")
    bd <- devianceFeatureSelection(input, assay = "counts", fam = "binomial")
    bd_df <- cbind.data.frame("gene_id"=rownames(bd),
        "gene_name"=rowData(bd)$gene_name, "dev"= rowData(bd)$binomial_deviance,
        "rank"=(nrow(bd) + 1) - rank(rowData(bd)$binomial_deviance))
    rownames(bd_df) <- bd_df$gene
    
    results_list <- vector("list", length(batch_effects))
    names(results_list) <- batch_effects  
    
    for (i in seq_along(batch_effects)) {
        batch <- batch_effects[i]
        
        if (verbose) message("Batch Effect: ", batch)
        if (verbose) message("Running feature selection without batch...")
        
        batch_data <- colData(input)[[batch]]
        bd_batch <- devianceFeatureSelection(input, assay = "counts",
            fam = "binomial",batch = as.factor(batch_data))
        
        bd_batch_df <- cbind.data.frame(
            "gene_id"=rownames(bd_batch),
            "gene_name"=rowData(bd_batch)$gene_name,
            "dev"= rowData(bd_batch)$binomial_deviance,
            "rank"=(nrow(bd_batch)+1)-rank(rowData(bd_batch)$binomial_deviance))
        rownames(bd_batch_df) <- bd_batch_df$gene
        
        if (verbose) message("Calculating deviance and rank difference...")
        
        batch_df <- left_join(bd_df, bd_batch_df, by=c("gene_id", "gene_name"),
            suffix=c("_default",paste0("_", batch)))
        
        batch_df$d_diff <- 
            (batch_df$dev_default- batch_df[[paste0("dev_", batch)]]) /
            batch_df[[paste0("dev_", batch)]]
        batch_df[[paste0("nSD_dev_", batch)]] <- 
            (batch_df$d_diff - mean(batch_df$d_diff)) / sd(batch_df$d_diff)
        
        batch_df$r_diff <- 
            batch_df[[paste0("rank_", batch)]] - batch_df$rank_default
        batch_df[[paste0("nSD_rank_", batch)]] <- 
            (batch_df$r_diff - mean(batch_df$r_diff)) / sd(batch_df$r_diff)
        
        results_list[[i]] <- batch_df
    }
    results_list
    }