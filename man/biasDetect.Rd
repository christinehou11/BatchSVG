% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biasDetect.R
\name{biasDetect}
\alias{biasDetect}
\title{Biased Genes Identification}
\usage{
biasDetect(
  list_batch_df,
  threshold = "both",
  nSD_dev = NULL,
  nSD_rank = NULL,
  plot_point_size = 3,
  plot_point_shape = 16,
  plot_text_size = 3,
  plot_pallete = "YlOrRd"
)
}
\arguments{
\item{list_batch_df}{\code{list} : The list of data frame(s) generated from
\code{featureSelection()} function. The length of the data frame list
should be at least one.}

\item{threshold}{A character string specifying the filtering criterion.
Must be one of:
\itemize{
\item \code{"dev"}: Filters genes based on the deviance threshold only.
\item \code{"rank"}: Filters genes based on the rank threshold only.
\item \code{"both"}: Filters genes based on either the deviance or rank threshold.
Default is "both".
}}

\item{nSD_dev}{\code{integer}: A numeric vector specifying the
number of standard deviation (nSD) for each batch when analyzing the
relative change in deviance. The order of values must correspond to
the order of batches in \code{list_batch_df}.
Required if \code{threshold} is "dev" or "both".
If a single value is provided, it is applied to all batches; otherwise,
it must have the same length as \code{list_batch_df}.}

\item{nSD_rank}{\code{vector}: A numeric vector specifying the
number of standard deviation (nSD) for each batch when analyzing rank
differences. The order of values must correspond to the order of batches
in \code{list_batch_df}.
Required if \code{threshold} is "rank" or "both".
If a single value is provided, it is applied to all batches; otherwise,
it must have the same length as \code{list_batch_df}.}

\item{plot_point_size}{\code{vector}: A numeric vector specifying point sizes
in plots. If asingle value is provided, it is applied to all batches.}

\item{plot_point_shape}{\code{vector}: A numeric vector specifying point
shapes in plots. If a single value is provided, it is applied to all
batches.}

\item{plot_text_size}{\code{vector}: A numeric vector specifying text label
size in plots. Default is \code{3}.}

\item{plot_pallete}{\code{vector}: A character string vector specifying the
color palette for plots. Default is \code{"YlOrRd"}.}
}
\value{
A named list where each element corresponds to a batch and contains:
\itemize{
\item \code{"Plot"}: A diagnostic plot (either deviance, rank, or both).
\item \code{"Table"}: A filtered data frame containing outlier genes based on the
specified threshold.
}
}
\description{
Function to identify the bias genes based on user-selected
threshold of number of standard deviation in relative change in deviance
and rank difference.
}
\examples{
library(ExperimentHub)
library(SummarizedExperiment)
library(tibble)

ehub <- ExperimentHub()
spe <- ehub[["EH9605"]]
fix_order <- dplyr::distinct(
    as.data.frame(colData(spe)), slide, array, brnum, sample_id, 
    position, sex) \%>\% 
    dplyr::arrange(slide, array)
sub4 <- fix_order$sample_id[c(14,16, 20,21)]
spe_sub4 <- spe[,spe$sample_id \%in\% sub4]

svgs_sub4 <- utils::read.csv(
    system.file("extdata","svgs_sub4.csv",package = "BatchSVG"),
    row.names = 1, check.names = FALSE)

spe_sub4 <- spe_sub4[rowData(spe_sub4)$gene_id \%in\% svgs_sub4$gene_id,]
rownames(spe_sub4) <- rowData(spe_sub4)$gene_id

SVGs <- svgs_sub4$gene_id
list_batch_df <- featureSelect(input = spe_sub4, 
   batch_effects = c("sample_id", "sex"), VGs = SVGs)

biaGenes <- biasDetect(list_batch_df = list_batch_df, threshold = "both", 
   nSD_dev = c(5,4), nSD_rank = c(4,6))
}
