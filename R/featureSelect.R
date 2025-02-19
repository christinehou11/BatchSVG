#' @rdname featureSelect
#'
#' @name featureSelect
#'
#' @title Feature Selection
#'
#' @description Function to conduct feature selection on spatial
#'     transcriptomics data and calculate the difference of deviance and rank
#'     with and without the selected batch effect.
#'
#' @importFrom scry devianceFeatureSelection
#' @importFrom dplyr left_join filter
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom stats sd
#'
#' @param input  \code{SpatialExperiment}: \code{SpatialExperiment} input data.
#'     It is assumed to have an \code{assays} slot
#'     containing \code{counts} assay for \code{devianceFeatureSelection()} to
#'     successfully operate and calculate the deviance and rank. The
#'     \code{logcounts} is strongly recommended to be included in \code{assays}.
#'     Also, the input data is assumed to set up \code{rownames(input)} as
#'     \code{genes} and \code{gene_name} in \code{rowData(input)}. The input can
#'     be either the raw and complete data object or the filtered data object
#'     containing only spatially varaible genes (SVGs).
#'
#' @param batch_effects \code{list}: Any batch effect (\code{slide} or
#'     \code{subject}) based on what metadata is available within input data.
#'     The name of the batch effect within each input object
#'     should be specified based on different scenarios.
#'
#' @param VGs \code{character}: Spatially Variable Genes (SVGs)
#'     for \code{SpatialExperiment} object. If it is a data frame, it is assumed
#'     to contain one column of identified variable genes with column name as
#'     gene name/ID, e.g. "ENSG00000002330".
#'
#' @return The output values are returned as
#'     a list of data frames containing the deviance and rank with and without 
#'     the bacth effect, the corresponding difference, the corresponding nSD, 
#'     gene, gene name, and whether the gene is outlier defined by the chosen 
#'     deviance and rank cutoff. The length of the list of data frames is equal
#'     to the number of batch effects included.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(SummarizedExperiment)
#' load(system.file("extdata","spe_sub4.rda",package = "BatchSVG"))
#' res_ranks <- read.csv(
#'     system.file("extdata","res_ranks.csv",package = "BatchSVG"),
#'     row.names = 1, check.names = FALSE)
#' res_df_sub4 <- tidyr::pivot_longer(
#'     tibble::rownames_to_column(as.data.frame(res_ranks), var<-"gene_id"), 
#'     colnames(res_ranks), 
#'     names_to = "sample_id", 
#'     values_to = "rank", 
#'     values_drop_na = TRUE)
#' # subset to 4 samples (rank <= 2000, n > 1)
#' res_df2_sub4 <- filter(res_df_sub4, 
#'     sample_id %in%
#'     c("V11L05-333_B1","V11L05-333_D1","V11L05-335_D1","V11L05-336_A1"),
#'     rank <= 2000)
#' svgs_sub4 <- group_by(res_df2_sub4, gene_id) %>% 
#'     tally() %>% 
#'     filter(n > 1)
#' SVGs <- svgs_sub4$gene_id
#' 
#' batch_df <- featureSelect(spe_sub4, 
#'     batch_effects = c("sample_id", "sex"), VGs = SVGs)
#' 
featureSelect <- function(input, batch_effects = NULL, VGs = NULL) {

    stopifnot(
        inherits(input, c("SpatialExperiment")),
        is.character(VGs)
    )

    if (!is.null(batch_effects)) {
        stopifnot(is.character(batch_effects))
        for (batch in batch_effects) {
            if (!batch %in% names(colData(input))) {
                stop(paste("The batch_effect", batch, "is not a valid column"))
            }
        }
    } else { 
        stop("Please provide a valid batch_effect.")
    }
    
    # Re-organize raw input
    input <- input[rowData(input)$gene_id %in% VGs, ]
    rownames(input) <- rowData(input)$gene_id
    
    message("Step 1: Running feature selection without batch...")
    bd <- devianceFeatureSelection(input, assay = "counts", fam = "binomial")
    bd_df <- cbind.data.frame("gene"=rownames(bd),
                            "gene_name"=rowData(bd)$gene_name,
                            "dev"= rowData(bd)$binomial_deviance,
                            "rank"=(nrow(bd)+1)
                                -rank(rowData(bd)$binomial_deviance))
    rownames(bd_df) <- bd_df$gene
    
    results_list <- list()
    
    for (batch in batch_effects) {
        message(paste("Batch Effect:", batch))
        message("Running feature selection...")
        
        batch_data <- colData(input)[[batch]]
        
        bd_batch <- devianceFeatureSelection(input, assay = "counts",
                    fam = "binomial",batch = as.factor(batch_data))
        bd_batch_df <- cbind.data.frame("gene"=rownames(bd_batch),
                                  "gene_name"=rowData(bd_batch)$gene_name,
                                  "dev"= rowData(bd_batch)$binomial_deviance,
                                  "rank"=(nrow(bd_batch)+1)
                                  -rank(rowData(bd_batch)$binomial_deviance))
        rownames(bd_batch_df) <- bd_batch_df$gene
    
        message("Calculating deviance and rank difference...")
        batch_df <- left_join(bd_df, bd_batch_df, by=c("gene", "gene_name"),
                                suffix=c("_default",paste0("_", batch)))
    
        batch_df$d_diff <- 
            (batch_df$dev_default- batch_df[[paste0("dev_", batch)]]) /
            batch_df[[paste0("dev_", batch)]]
        mean_dev <- mean(batch_df$d_diff)
        sd_dev <- sd(batch_df$d_diff)
        batch_df[[paste0("nSD_dev_", batch)]] <- 
            (batch_df$d_diff - mean_dev) / sd_dev
    
        batch_df$r_diff <- 
            batch_df[[paste0("rank_", batch)]] - batch_df$rank_default
        mean_rank <- mean(batch_df$r_diff)
        sd_rank <- sd(batch_df$r_diff)
        batch_df[[paste0("nSD_rank_", batch)]] <- 
            (batch_df$r_diff - mean_rank) / sd_rank
        
        results_list[[batch]] <- batch_df
    }

    results_list
    }


