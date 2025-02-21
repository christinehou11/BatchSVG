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
#' @importFrom ggplot2 ggplot geom_point scale_color_manual scale_y_log10
#'            geom_abline theme_bw aes labs scale_x_log10 scale_y_reverse
#'            coord_fixed labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr text_grob annotate_figure ggarrange
#'
#' @param list_batch_df \code{list} : The list of data frame(s) generated from
#'    `featureSelection()` function. The length of the data frame list
#'    should be at least one.
#'
#' @param nSD_dev \code{integer}: A numeric vector specifying the 
#'    number of standard deviation (nSD) for each batch when analyzing the 
#'    relative change in deviance. The order of values must correspond to 
#'    the order of batches in `list_batch_df`.
#'    If a single value is provided, it is applied to all batches; otherwise, 
#'    it must have the same length as `list_batch_df`.
#'
#' @param nSD_rank \code{vector}: A numeric vector specifying the 
#'    number of standard deviation (nSD) for each batch when analyzing rank 
#'    differences. The order of values must correspond to the order of batches 
#'    in `list_batch_df`.
#'    If a single value is provided, it is applied to all batches; otherwise, 
#'    it must have the same length as `list_batch_df`.
#'
#' @param visual \code{logical}: Whether to display the detected bias genes by
#'    visualizations of \code{deviance} and \code{rank}. Default is FALSE.
#'    If visual = TRUE, the returned format will be two parallel plots 
#'    presenting bias genes based on both \code{nSD_dev} and \code{nSD_rank}.
#'
#' @return The named list where each element corresponds to a batch. 
#'    Each batch contains:
#'    \itemize{
#'      \item \strong{"Plot"}: A combined ggplot object displaying two scatter 
#'          plots—one for deviance and one for rank—with an enforced 1:1 aspect 
#'          ratio.
#'      \item \strong{"Table"}: A data frame containing identified biased genes,
#'           with columns indicating whether they are outliers based on 
#'           deviance or rank thresholds.
#'    }
#'    The order of batches in the list follows the order of `nSD_dev` and 
#'    `nSD_rank` provided by the user.
#'
#' @export
#'
#' @examples
#' load(system.file("extdata","list_batch_df.rda",package = "BatchSVG"))
#' 
#' biaGenes <- biasDetect(list_batch_df, nSD_dev = c(5,4), nSD_rank = c(4,6))
#' 
#' biaGenes[["sample_id"]][["Table"]] # see biased genes in table format
#' biaGenes[["sex"]][["Plot"]] # see biased genes in plot format

biasDetect <- function(list_batch_df, nSD_dev, nSD_rank) {
    # input check
    num_batches <- length(list_batch_df)
    if (length(nSD_dev) == 1) {
        nSD_dev <- rep(nSD_dev, num_batches)}
    if (length(nSD_rank) == 1) {
        nSD_rank <- rep(nSD_rank, num_batches)}
    stopifnot(
        is.list(list_batch_df), length(list_batch_df) > 0,
        all(nSD_dev == floor(nSD_dev), 
            nSD_rank == floor(nSD_rank)),
        length(nSD_dev) == num_batches, 
        length(nSD_rank) == num_batches)
    
    biased_list <- list()
    
    for (i in seq_along(list_batch_df)) {
        batch <- names(list_batch_df)[i]
        batch_df <- list_batch_df[[batch]]
        stopifnot(is.data.frame(batch_df))
        
        # select corresponding values
        sd_dev <- nSD_dev[i]
        sd_rank <- nSD_rank[i]
        
        # plots
        dev_colname <- paste0("nSD_dev_",batch)
        batch_df$nSD_bin_dev <- cut(abs(batch_df[[dev_colname]]), right = FALSE,
            breaks=seq(0,max(batch_df[[dev_colname]]) + sd_dev, 
            by=sd_dev), include.lowest=TRUE)
        col_pal_dev <- brewer.pal(length(unique(batch_df[["nSD_bin_dev"]])), 
                                "YlOrRd")
        col_pal_dev[1] <- "grey"
        dev_sd_plot <- ggplot(batch_df,
            aes(x = .data[["dev_default"]], 
                y = .data[[paste0("dev_", batch)]],
                color = .data[["nSD_bin_dev"]])) +
            geom_point() + scale_x_log10() + scale_y_log10() +
            scale_color_manual(values=col_pal_dev) +
            geom_text_repel(
                aes(label = ifelse(.data[[dev_colname]] > sd_dev,
                .data[["gene_name"]], "")), size = 3, max.overlaps = Inf) +
            labs(x= "dev (no batch)", y="dev (batch)", color = "nSD Interval",
                title = "Condition: deviance") +
            theme_bw() + coord_fixed(ratio = 1) +
            geom_abline(aes(slope = 1, intercept = 0), lty = 2)

        rank_colname <- paste0("nSD_rank_", batch)
        batch_df$nSD_bin_rank <- cut(abs(batch_df[[rank_colname]]), right=FALSE,
            breaks=seq(0,max(batch_df[[rank_colname]]) + sd_rank, 
            by=sd_rank),include.lowest=TRUE)
        col_pal_rank <- brewer.pal(length(unique(batch_df$nSD_bin_rank)), 
                                  "YlOrRd")
        col_pal_rank[1] <- "grey"
        rank_sd_plot <- ggplot(batch_df, 
            aes(x = .data[["rank_default"]],
                y = .data[[paste0("rank_", batch)]],
                color = .data[["nSD_bin_rank"]])) +
            geom_point() + scale_y_reverse() + 
            scale_color_manual(values = col_pal_rank) +
            geom_text_repel(
                aes(label = ifelse(.data[[rank_colname]] > sd_rank,
                .data[["gene_name"]], "")), size = 3, max.overlaps = Inf) +
            geom_abline(aes(slope = -1, intercept = 0), lty = 2) +
            labs(x= "rank (no batch)", y= "rank (batch)", color = "nSD_bin",
                title = "Condition: rank") +
            theme_bw() + coord_fixed(ratio = 1) 

        batch_plots <- ggarrange(dev_sd_plot, rank_sd_plot, ncol = 2,
            common.legend = TRUE, legend = "bottom")
        biased_list[[batch]][["Plot"]] <- annotate_figure(batch_plots, 
            top = text_grob(paste0("Identified Biased Genes - Batch: ", batch),
            face = "bold", size = 16))
        
        # data frame
        batch_df$dev_outlier <- batch_df$nSD_dev >= nSD_dev
        batch_df$rank_outlier <- batch_df$nSD_rank >= nSD_rank
        biased_list[[batch]][["Table"]] <- filter(batch_df,
            .data$dev_outlier==TRUE | .data$rank_outlier==TRUE)
    }
    biased_list
    }
