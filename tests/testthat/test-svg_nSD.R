suppressPackageStartupMessages({
library(ExperimentHub)
library(SummarizedExperiment)
library(tibble)
})

ehub <- ExperimentHub()
spe <- ehub[["EH9605"]]
fix_order <- dplyr::distinct(
    as.data.frame(colData(spe)), slide, array, brnum, sample_id, 
    position, sex) %>% 
    dplyr::arrange(slide, array)
sub4 <- fix_order$sample_id[c(14,16, 20,21)]
spe_sub4 <- spe[,spe$sample_id %in% sub4]

svgs_sub4 <- utils::read.csv(
    system.file("extdata","svgs_sub4.csv",package = "BatchSVG"),
    row.names = 1, check.names = FALSE)

spe_sub4 <- spe_sub4[rowData(spe_sub4)$gene_id %in% svgs_sub4$gene_id,]
rownames(spe_sub4) <- rowData(spe_sub4)$gene_id

SVGs <- svgs_sub4$gene_id
list_batch_df <- featureSelect(input = spe_sub4, 
    batch_effects = c("sample_id", "sex"), VGs = SVGs)

plots <- svg_nSD(list_batch_df = list_batch_df, sd_interval_dev = c(4,7),
    sd_interval_rank = c(6,7))

test_that("The example results have the correct number of batch effects", {
    expect_true(length(plots) == 2)
})