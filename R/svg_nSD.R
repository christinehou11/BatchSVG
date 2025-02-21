#' @rdname svg_nSD
#'
#' @name svg_nSD
#'
#' @title SVGs Plots with Relative Change in Deviance and Rank Difference
#'
#' @description Function to the spatially variable genes in relative change in
#'    deviance and rank difference colored by the number of standard deviation.
#'
#' @importFrom ggplot2 ggplot geom_histogram scale_fill_manual labs
#'                scale_y_continuous theme_bw theme margin geom_point
#'                scale_x_log10 scale_y_log10 scale_color_manual
#'                geom_abline scale_y_reverse
#' @importFrom gridExtra grid.arrange
#' @importFrom ggpubr ggarrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales pseudo_log_trans
#'
#' @param list_batch_df \code{list} : The list of data frame(s) generated from
#'    `featureSelection()` function. The length of the data frame list
#'    should be at least one.
#'    
#' @param sd_interval_dev \code{vector}: A numeric vector specifying the 
#'    interval widths for standard deviation bins when analyzing the 
#'    relative change in deviance. 
#'    If a single value is provided, it is applied to all batches; otherwise, 
#'    it must have the same length as `list_batch_df`.
#'
#' @param sd_interval_rank \code{vector}: A numeric vector specifying the 
#'    interval widths for standard deviation bins when analyzing rank 
#'    differences.
#'    If a single value is provided, it is applied to all batches; otherwise, 
#'    it must have the same length as `list_batch_df`.
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
#' 
#' load(system.file("extdata","list_batch_df.rda",package = "BatchSVG"))
#' 
#' plots <- svg_nSD(list_batch_df, 
#'     sd_interval_dev = c(5,4), sd_interval_rank = c(4,6))
#' 

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
        
        # select corresponding values
        sd_dev <- sd_interval_dev[i]
        sd_rank <- sd_interval_rank[i]
        
        # deviance plot
        batch_df$nSD_bin_dev <- cut(abs(batch_df$nSD_dev), right = FALSE,
            breaks=seq(0,max(batch_df$nSD_dev) + sd_dev, 
                    by=sd_dev), include.lowest=TRUE)
        
        col_pal_dev <- brewer.pal(length(unique(batch_df[["nSD_bin_dev"]])), 
                                    "YlOrRd")
        col_pal_dev[1] <- "grey"
        
        dev_sd_plot1 <- ggplot(batch_df, 
            aes(x = .data[["d_diff"]], fill=  .data[["nSD_bin_dev"]])) +
            geom_histogram(color = "grey20", bins=50) +
            scale_fill_manual(values = col_pal_dev) +
            labs(x = "\u0394 deviance", y = "# SVGs", 
                fill = "nSD Interval   ") +
            scale_y_continuous(trans = pseudo_log_trans(sigma = 1),
                breaks = 10^(0:4), labels = format(10^(0:4), scientific = F)) +
            theme_bw()
        
        dev_sd_plot2 <- ggplot(batch_df,
            aes(x = .data[["dev_default"]], 
                y = .data[[paste0("dev_", batch)]],
                color = .data[["nSD_bin_dev"]])) +
            geom_point() + scale_x_log10() + scale_y_log10() +
            scale_color_manual(values=col_pal_dev) +
            labs(x= "dev (no batch)", y="dev (batch)") +
            theme_bw() +
            geom_abline(aes(slope = 1, intercept = 0), lty = 2)
        
        dev_sd_plots <- ggarrange(dev_sd_plot1, dev_sd_plot2,
                        ncol = 2,common.legend = TRUE, legend="right")
        tg_dev <- textGrob(paste0("SVGs with relative change in deviance"), 
                        gp = gpar(fontsize = 13, fontface = 'bold'))
        sg_dev <- textGrob(paste0("Batch: ", batch,
                                    "; nSD width = ", sd_dev), 
                        gp = gpar(fontsize = 10))
        margin <- unit(0.5, "line")
        dev_sd_plots <- grid.arrange(tg_dev, sg_dev, dev_sd_plots,
            heights = unit.c(grobHeight(tg) + 1.2*margin, 
                        grobHeight(sg) + margin, unit(1,"null")))
        
        # rank plot
        batch_df$nSD_bin_rank <- cut(abs(batch_df$nSD_rank), right=FALSE,
            breaks=seq(0,max(batch_df$nSD_rank) + sd_rank, 
                    by=sd_rank),include.lowest=TRUE)
        
        col_pal_rank <- brewer.pal(length(unique(batch_df$nSD_bin_rank)), 
                                    "YlOrRd")
        col_pal_rank[1] <- "grey"
        
        rank_sd_plot1 <- ggplot(batch_df, 
            aes(x = .data[["r_diff"]], fill = .data[["nSD_bin_rank"]])) +
            geom_histogram(color = "grey20", bins = 50) +
            scale_fill_manual(values = col_pal_rank) +
            labs(x = "rank difference", y = "# SVGs",
                fill = "nSD Interval   ") +
            scale_y_continuous(trans = pseudo_log_trans(sigma = 1),
                breaks = 10^(0:4), labels = format(10^(0:4), scientific = F)) +
            theme_bw()
        
        rank_sd_plot2 <- ggplot(batch_df, 
            aes(x = .data[["rank_default"]],
                y = .data[[paste0("rank_", batch)]],
                color = .data[["nSD_bin_rank"]])) +
            geom_point() + scale_y_reverse() + 
            scale_color_manual(values = col_pal_rank) +
            geom_abline(aes(slope = -1, intercept = 0), lty = 2) +
            labs(x= "rank (no batch)", y= "rank (batch)") +
            theme_bw()
        
        rank_sd_plots <- ggarrange(rank_sd_plot1, rank_sd_plot2,
                        ncol = 2,common.legend = TRUE, legend="right")
        tg_rank <- textGrob(paste0("SVGs with rank difference"), 
                        gp = gpar(fontsize = 13, fontface = 'bold'))
        sg_rank <- textGrob(paste0("Batch: ", batch,
                                  "; nSD width = ", sd_rank), 
                           gp = gpar(fontsize = 10))
        rank_sd_plots <- grid.arrange(tg_rank, sg_rank, rank_sd_plots,
            heights = unit.c(grobHeight(tg) + 1.2*margin, 
                        grobHeight(sg) + margin, unit(1,"null")))
        
        # combine deviance and rank plots
        batch_plots <- grid.arrange(dev_sd_plots, rank_sd_plots, nrow = 2)
        plot_list[[batch]] <- batch_plots
    }

    final_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))
    
    final_plot
    }
