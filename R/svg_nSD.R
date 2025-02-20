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
#'                scale_y_continuous theme_bw theme margin
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales pseudo_log_trans
#'
#' @param list_batch_df \code{list} : The list of data frame(s) generated from
#'    `featureSelection()` function. The length of the data frame list
#'    should be at least one.
#'    
#' @param sd_interval_dev \code{integer}: Interval of standard deviation on 
#'    relative change in deviance. The default value is 5.
#'
#' @param sd_interval_rank \code{integer}: Number of standard deviation on rank
#'     difference. The default value is 5.
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
#' sample_id_df <- utils::read.csv(
#'     system.file("extdata","sample_id_df.csv",package = "BatchSVG"),
#'     row.names = 1, check.names = FALSE)
#' 
#' plots <- svg_nSD(sample_id_df, nSD_dev = 5, nSD_rank)
#' 

svg_nSD <- function(list_batch_df, 
                        sd_interval_dev = 5, sd_interval_rank = 5) {
  
    stopifnot(
        is.list(list_batch_df), length(list_batch_df) > 0,
        all(sd_interval_dev == floor(sd_interval_dev), 
            sd_interval_rank == floor(sd_interval_rank))
    )
    
    plot_list <- list()
    
    for (batch in names(list_batch_df)) {
        
        batch_df <- list_batch_df[[batch]]
        stopifnot(is.data.frame(batch_df))
        
        # deviance plot
        batch_df$nSD_bin_dev <- cut(abs(batch_df$nSD_dev), right = FALSE,
            breaks=seq(0,max(batch_df$nSD_dev) + sd_interval_dev, 
                    by=sd_interval_dev), include.lowest=TRUE)
        
        col_pal_dev <- brewer.pal(length(unique(batch_df$nSD_bin_dev)), 
                                    "YlOrRd")
        col_pal_dev[1] <- "grey"
        
        dev_sd_plot <- ggplot(batch_df, 
                            aes(x = d_diff, fill=  nSD_bin_dev)) +
            geom_histogram(color = "grey20", bins=50) +
            scale_fill_manual(values = col_pal_dev) +
            labs(x = "\u0394 deviance", y = "# SVGs", 
                fill = paste("nSD Interval: ", batch),
                title = paste("nSD bin width = ", sd_interval_dev), 
                subtitle = paste("Batch: ", batch)) +
            scale_y_continuous(trans = pseudo_log_trans(sigma = 1),
                breaks = 10^(0:4), labels = format(10^(0:4), scientific = F)) +
            theme_bw() + 
            theme(legend.position = "right", plot.margin = margin(5, 10, 5, 5))
        
        # rank plot
        batch_df$nSD_bin_rank <- cut(abs(batch_df$nSD_rank), right=FALSE,
            breaks=seq(0,max(batch_df$nSD_rank) + sd_interval_rank, 
                    by=sd_interval_rank),include.lowest=TRUE)
        
        col_pal_rank <- brewer.pal(length(unique(batch_df$nSD_bin_rank)), 
                                    "YlOrRd")
        col_pal_rank[1] <- "grey"
        
        rank_sd_plot <- ggplot(batch_df, 
                            aes(x = r_diff, fill = nSD_bin_rank)) +
            geom_histogram(color = "grey20", bins = 50) +
            scale_fill_manual(values = col_pal_rank) +
            labs(x = "rank difference", y = NULL, 
                title = paste("nSD bin width = ", sd_interval_rank),
                subtitle = paste("Batch: ", batch)) +
            scale_y_continuous(trans = pseudo_log_trans(sigma = 1),
                breaks = 10^(0:4), labels = format(10^(0:4), scientific = F)) +
            theme_bw() + 
            theme(legend.position = "none")
        
        # combine deviance and rank plots
        plot_list[[batch]] <- arrangeGrob(dev_sd_plot, rank_sd_plot, 
                                ncol = 2, widths = c(1, 1))
    }

    final_plot <- do.call(grid.arrange, c(plot_list, nrow = 2))
    
    final_plot
    }
