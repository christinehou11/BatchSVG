#' @rdname biasDetect
#'
#' @name biasDetect
#'
#' @title Biased Genes Identification
#'
#' @description Function to identify the bias genes based on user-selected
#'    threshold of number of standard deviation in relative change in deviance 
#'    and rank difference.
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot scale_color_manual aes labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom RColorBrewer brewer.pal
#'
#' @param list_batch_df \code{list} : The list of data frame(s) generated from
#'    `featureSelection()` function. The length of the data frame list
#'    should be at least one.
#' @param threshold A character string specifying the filtering criterion. 
#'    Must be one of:
#'    - `"dev"`: Filters genes based on the deviance threshold only.
#'    - `"rank"`: Filters genes based on the rank threshold only.
#'    - `"both"`: Filters genes based on either the deviance or rank threshold.
#'        Default is "both".
#'
#' @param nSD_dev \code{integer}: A numeric vector specifying the 
#'    number of standard deviation (nSD) for each batch when analyzing the 
#'    relative change in deviance. The order of values must correspond to 
#'    the order of batches in `list_batch_df`.
#'    Required if `threshold` is "dev" or "both".
#'    If a single value is provided, it is applied to all batches; otherwise, 
#'    it must have the same length as `list_batch_df`.
#'
#' @param nSD_rank \code{vector}: A numeric vector specifying the 
#'    number of standard deviation (nSD) for each batch when analyzing rank 
#'    differences. The order of values must correspond to the order of batches 
#'    in `list_batch_df`.
#'    Required if `threshold` is "rank" or "both".
#'    If a single value is provided, it is applied to all batches; otherwise, 
#'    it must have the same length as `list_batch_df`.
#'
#' @return A named list where each element corresponds to a batch and contains:
#'   - `"Plot"`: A diagnostic plot (either deviance, rank, or both).
#'   - `"Table"`: A filtered data frame containing outlier genes based on the 
#'   specified threshold.
#'
#' @export
#'
#' @examples
#' data(spe_sub4)
#' spe_sub4
#' 
#' data(svgs_sub4)
#' SVGs <- svgs_sub4$gene_id
#' 
#' list_batch_df <- featureSelect(input = spe_sub4, 
#'     batch_effects = c("sample_id", "sex"), VGs = SVGs)
#' 
#' biaGenes <- biasDetect(list_batch_df = list_batch_df, threshold = "both", 
#'     nSD_dev = c(5,4), nSD_rank = c(4,6))
biasDetect <- function(list_batch_df, threshold = "both", 
                        nSD_dev = NULL, nSD_rank = NULL) {
    
    filter_condition <- match.arg(threshold, choices = c("dev", "rank", "both"))
    if (threshold %in% c("dev", "both") && is.null(nSD_dev)) {
        stop("When threshold = 'dev' or 'both', you must specify 'nSD_dev'.")}
    if (threshold %in% c("rank", "both") && is.null(nSD_rank)) {
        stop("When threshold = 'rank' or 'both', you must specify 'nSD_rank'.")}
    
    num_batches <- length(list_batch_df)
    if (!is.null(nSD_dev)) {
        if (!is.numeric(nSD_dev) || any(nSD_dev != floor(nSD_dev))) {
            stop("'nSD_dev' must be an integer.")}
        if (length(nSD_dev) == 1) { nSD_dev <- rep(nSD_dev, num_batches)}}
    
    if (!is.null(nSD_rank)) {
        if (!is.numeric(nSD_rank) || any(nSD_rank != floor(nSD_rank))) {
            stop("'nSD_rank' must be an integer.")}
        if (length(nSD_rank) == 1) { nSD_rank <- rep(nSD_rank, num_batches)}}
    
    stopifnot(is.list(list_batch_df), length(list_batch_df) > 0)
    biased_list <- list()
    
    for (i in seq_along(list_batch_df)) {
        batch <- names(list_batch_df)[i]
        batch_df <- list_batch_df[[batch]]
        stopifnot(is.data.frame(batch_df))
        
        if (!is.null(nSD_dev)) {
        sd_dev <- nSD_dev[i]
        dev_colname <- paste0("nSD_dev_",batch)
        batch_df$nSD_bin_dev <- cut(abs(batch_df[[dev_colname]]), right = FALSE,
            breaks=seq(0,max(batch_df[[dev_colname]]) + sd_dev, 
            by=sd_dev), include.lowest=TRUE)
        col_pal_dev <- brewer.pal(length(unique(batch_df[["nSD_bin_dev"]])), 
            "YlOrRd")
        col_pal_dev[1] <- "grey"
        dev_sd_plot <- ggplot(batch_df,
            aes(x = .data[["dev_default"]], y = .data[[paste0("dev_", batch)]],
                color = .data[["nSD_bin_dev"]]))
        dev_sd_plot <- .theme_dev_point_plot(dev_sd_plot) +
            scale_color_manual(values=col_pal_dev) +
            labs(subtitle = paste0("Batch: ", batch, "; nSD width: ", sd_dev))+
            geom_text_repel(
                aes(label = ifelse(.data[[dev_colname]] > sd_dev,
                .data[["gene_name"]], "")), size = 3, max.overlaps = Inf)
        batch_df$dev_outlier <- batch_df$nSD_dev >= sd_dev
        }
        
        if (!is.null(nSD_rank)) {
        sd_rank <- nSD_rank[i]
        rank_colname <- paste0("nSD_rank_", batch)
        batch_df$nSD_bin_rank <- cut(abs(batch_df[[rank_colname]]), right=FALSE,
            breaks=seq(0,max(batch_df[[rank_colname]]) + sd_rank, 
            by=sd_rank),include.lowest=TRUE)
        col_pal_rank <- brewer.pal(length(unique(batch_df$nSD_bin_rank)), 
            "YlOrRd")
        col_pal_rank[1] <- "grey"
        rank_sd_plot <- ggplot(batch_df, 
            aes(x = .data[["rank_default"]],y = .data[[paste0("rank_", batch)]],
                color = .data[["nSD_bin_rank"]]))
        rank_sd_plot <- .theme_rank_point_plot(rank_sd_plot) +
            scale_color_manual(values = col_pal_rank) +
            labs(subtitle = paste0("Batch: ", batch, "; nSD width: ", sd_rank))+
            geom_text_repel(
                aes(label = ifelse(.data[[rank_colname]] > sd_rank,
                .data[["gene_name"]], "")), size = 3, max.overlaps = Inf)
        batch_df$rank_outlier <- batch_df$nSD_rank >= sd_rank
        }
        
        if (filter_condition == "dev") {
            biased_list[[batch]][["Plot"]] <- dev_sd_plot
            biased_list[[batch]][["Table"]] <- filter(batch_df,
                .data$dev_outlier==TRUE)}
        else if (filter_condition == "rank") {
            biased_list[[batch]][["Plot"]] <- rank_sd_plot
            biased_list[[batch]][["Table"]] <- filter(batch_df,
                .data$rank_outlier==TRUE)}
        else {
            biased_list[[batch]][["Plot"]] <- plot_grid(
                dev_sd_plot, rank_sd_plot, ncol = 2, align = "hv")
            biased_list[[batch]][["Table"]] <- filter(batch_df,
                .data$dev_outlier==TRUE | .data$rank_outlier==TRUE)}
    }
    biased_list
    }
