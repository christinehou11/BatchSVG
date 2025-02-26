# helper function for plots

#' @importFrom ggplot2 scale_y_continuous element_text theme theme_bw
#'                    geom_histogram labs 
#' @importFrom scales pseudo_log_trans
.theme_dev_bar_plot <- function(plot) {
    plot +
    geom_histogram(color = "grey20", bins = 50) +
    scale_y_continuous(trans = pseudo_log_trans(sigma = 1),
        breaks = 10^(0:4), labels = format(10^(0:4), scientific = F)) +
    theme_bw() + 
    theme(legend.position = "right", aspect.ratio = 1,
        plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
    labs(x = "\u0394 deviance", y = "# SVGs", 
         fill = "nSD Deviance Interval",
         title = "SVGs with relative change in deviance")
    }

#' @importFrom ggplot2 scale_y_continuous element_text theme theme_bw
#'                    geom_histogram labs 
#' @importFrom scales pseudo_log_trans
.theme_rank_bar_plot <- function(plot) {
    plot +
    geom_histogram(color = "grey20", bins = 50) +
    scale_y_continuous(trans = pseudo_log_trans(sigma = 1),
        breaks = 10^(0:4), labels = format(10^(0:4), scientific = F)) +
    theme_bw() + 
    theme(legend.position = "right", aspect.ratio = 1,
        plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
    labs(x = "rank difference", y = "# SVGs",
         fill = "nSD Rank Interval",
         title = "SVGs with rank difference")
    }

#' @importFrom ggplot2 geom_point scale_x_log10 scale_y_log10 geom_abline
#'                theme_bw theme element_text labs aes
.theme_dev_point_plot <- function(plot) {
    plot +
    geom_point() + scale_x_log10() + scale_y_log10() +
    geom_abline(aes(slope = 1, intercept = 0), lty = 2) +
    theme_bw() + 
    theme(legend.position = "right", aspect.ratio = 1,
        plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
    labs(x= "dev (no batch)", y="dev (batch)", 
        color = "nSD Deviance Interval",
        title = "Deviance without vs. with batch")
    }

#' @importFrom ggplot2 geom_point scale_y_reverse geom_abline aes theme_bw
#'                    element_text labs
.theme_rank_point_plot <- function(plot) {
    plot +
    geom_point() + scale_y_reverse() + 
    geom_abline(aes(slope = -1, intercept = 0), lty = 2) +
    theme_bw() + 
    theme(legend.position = "right", aspect.ratio = 1,
        plot.title = element_text(size = 12, face = "bold"), 
        plot.subtitle = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
    labs(x= "dev (no batch)", y="dev (batch)", 
         color = "nSD Rank Interval",
         title = "Rank without vs. with batch")
    }