#' @rdname svg_nSD
#'
#' @name svg_nSD
#'
#' @title SVGs Plots
#'
#' @description This function generates visualizations to assess the impact 
#'    of batch effects on spatially variable genes (SVGs) by analyzing changes 
#'    in deviance and rank. The function bins the deviations into normalized 
#'    standard deviation (nSD) intervals and creates histograms and scatter 
#'    plots to illustrate the distribution of batch effects.
#'
#' @importFrom ggplot2 ggplot scale_fill_manual labs theme_bw aes
#' @importFrom rlang .data
#' @importFrom cowplot plot_grid
#' @importFrom RColorBrewer brewer.pal
#'
#' @param list_batch_df A named list of data frames, where each data frame 
#'    corresponds to a batch effect and contains columns with deviance and 
#'    rank differences.
#'    
#' @param sd_interval_dev A numeric vector specifying the 
#'    interval widths for standard deviation bins for each batch when analyzing 
#'    the relative change in deviance. The order of values must correspond to 
#'    the order of batches in \code{list_batch_df} . 
#'    If a single value is provided, it is applied to all batches; otherwise, 
#'    it must have the same length as \code{list_batch_df}.
#'
#' @param sd_interval_rank \code{vector}: A numeric vector specifying the 
#'    interval widths for standard deviation bins when analyzing rank 
#'    differences. The order of values must correspond to the order of batches 
#'    in \code{list_batch_df}.
#'    If a single value is provided, it is applied to all batches; otherwise, 
#'    it must have the same length as \code{list_batch_df}.
#'
#' @return A combined \code{ggplot} object containing:
#'    \itemize{
#'        \item \strong{Deviance Plots}:
#'        \itemize{
#'          \item Histogram of deviance differences across SVGs, colored by 
#'          nSD intervals.
#'          \item Scatter plot comparing deviance values before and after 
#'          batch correction.
#'        }
#'      \item \strong{Rank Plots}:
#'        \itemize{
#'          \item Histogram of rank differences across SVGs, colored by 
#'          nSD intervals.
#'          \item Scatter plot comparing ranks before and after batch 
#'          correction.
#'        }
#'    }
#'    The function arranges plots for each batch in a grid format for 
#'    easy comparison.
#'
#' @export
#'
#' @examples
#' suppressPackageStartupMessages({
#' library(ExperimentHub)
#' library(SummarizedExperiment)
#' library(tibble)
#' })
#' 
#' ehub <- ExperimentHub()
#' spe <- ehub[["EH9605"]]
#' fix_order <- dplyr::distinct(
#'     as.data.frame(colData(spe)), slide, array, brnum, sample_id, 
#'     position, sex) %>% 
#'     dplyr::arrange(slide, array)
#' sub4 <- fix_order$sample_id[c(14,16, 20,21)]
#' spe_sub4 <- spe[,spe$sample_id %in% sub4]
#' 
#' svgs_sub4 <- utils::read.csv(
#'     system.file("extdata","svgs_sub4.csv",package = "BatchSVG"),
#'     row.names = 1, check.names = FALSE)
#' 
#' spe_sub4 <- spe_sub4[rowData(spe_sub4)$gene_id %in% svgs_sub4$gene_id,]
#' rownames(spe_sub4) <- rowData(spe_sub4)$gene_id
#' 
#' SVGs <- svgs_sub4$gene_id
#' list_batch_df <- featureSelect(input = spe_sub4, 
#'    batch_effects = c("sample_id", "sex"), VGs = SVGs)
#' 
#' plots <- svg_nSD(list_batch_df = list_batch_df, 
#'    sd_interval_dev = c(5,4), sd_interval_rank = c(4,6))
svg_nSD <- function(list_batch_df, sd_interval_dev, sd_interval_rank) {
    # input check
    num_batches <- length(list_batch_df)
    if (length(sd_interval_dev) == 1) {
        sd_interval_dev <- rep(sd_interval_dev, num_batches)}
    if (length(sd_interval_rank) == 1) {
        sd_interval_rank <- rep(sd_interval_rank, num_batches)}
    stopifnot(
        is.list(list_batch_df), length(list_batch_df) > 0,
        all(sd_interval_dev == floor(sd_interval_dev), 
            sd_interval_rank == floor(sd_interval_rank)),
        length(sd_interval_dev) == num_batches, 
        length(sd_interval_rank) == num_batches)
    
    plot_list <- list()
    
    for (i in seq_along(list_batch_df)) {
        batch <- names(list_batch_df)[i]
        batch_df <- list_batch_df[[batch]]
        stopifnot(is.data.frame(batch_df))
        
        # deviance plot
        sd_dev <- sd_interval_dev[i]
        dev_colname <- paste0("nSD_dev_",batch)
        batch_df$nSD_bin_dev <- cut(abs(batch_df[[dev_colname]]), right = FALSE,
            breaks=seq(0,max(batch_df[[dev_colname]]) + sd_dev, 
            by=sd_dev), include.lowest=TRUE)
        col_pal_dev <- brewer.pal(length(unique(batch_df[["nSD_bin_dev"]])), 
                                "YlOrRd")
        col_pal_dev[1] <- "grey"
        
        dev_sd_plot1 <- ggplot(batch_df, 
            aes(x = .data[["d_diff"]], fill=  .data[["nSD_bin_dev"]]))
        dev_sd_plot1 <- .theme_dev_bar_plot(dev_sd_plot1) + 
            scale_fill_manual(values = col_pal_dev) +
            labs(subtitle = paste0("Batch: ", batch, "; nSD width = ", sd_dev))
        
        dev_sd_plot2 <- ggplot(batch_df,
            aes(x = .data[["dev_default"]],y = .data[[paste0("dev_", batch)]],
                color = .data[["nSD_bin_dev"]]))
        dev_sd_plot2 <- .theme_dev_point_plot(dev_sd_plot2,
            point_size = 3, point_shape = 16) + 
            scale_color_manual(values=col_pal_dev) +
            labs(subtitle = paste0("Batch: ",batch, "; nSD width = ", sd_dev))
        
        
        # rank plot
        sd_rank <- sd_interval_rank[i]
        rank_colname <- paste0("nSD_rank_", batch)
        batch_df$nSD_bin_rank <- cut(abs(batch_df[[rank_colname]]), right=FALSE,
            breaks=seq(0,max(batch_df[[rank_colname]]) + sd_rank, 
            by=sd_rank),include.lowest=TRUE)
        col_pal_rank <- brewer.pal(length(unique(batch_df$nSD_bin_rank)), 
            "YlOrRd")
        col_pal_rank[1] <- "grey"
        
        rank_sd_plot1 <- ggplot(batch_df, 
            aes(x = .data[["r_diff"]], fill = .data[["nSD_bin_rank"]]))
        rank_sd_plot1  <- .theme_rank_bar_plot(rank_sd_plot1) + 
            scale_fill_manual(values = col_pal_rank) +
            labs(subtitle = paste0("Batch: ", batch, "; nSD width = ", sd_rank))
        
        rank_sd_plot2 <- ggplot(batch_df, 
            aes(x = .data[["rank_default"]],y = .data[[paste0("rank_", batch)]],
                color = .data[["nSD_bin_rank"]]))
        rank_sd_plot2 <- .theme_rank_point_plot(rank_sd_plot2,
            point_size = 3, point_shape = 16) + 
            scale_color_manual(values = col_pal_rank) +
            labs(subtitle = paste0("Batch: ", batch,"; nSD width = ", sd_rank))
        
        plot_list[[batch]] <- plot_grid(dev_sd_plot1, dev_sd_plot2, 
            rank_sd_plot1, rank_sd_plot2, ncol = 2, align = "hv")
    }
    plot_list
    }
