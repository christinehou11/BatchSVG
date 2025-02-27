library(ExperimentHub)
library(SummarizedExperiment)
library(tibble)

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

bias_dev <- biasDetect(list_batch_df = list_batch_df, threshold = "dev",
        nSD_dev = 7)

bias_rank <- biasDetect(list_batch_df = list_batch_df, threshold = "rank",
        nSD_rank = 6)

bias_both <- biasDetect(list_batch_df = list_batch_df, threshold = "both",
        nSD_dev = c(4,8), nSD_rank = c(6,5))

test_that("Both example results have two batch effects" , {
    expect_true(length(bias_dev) == 2)
    expect_true(length(bias_rank) == 2)
    expect_true(length(bias_both) == 2)
})

test_that("All example results have table and plot saved for each batch", {
    expect_true(is.data.frame(bias_dev$sample_id$Table))
    expect_true("Plot" %in% names(bias_dev$sample_id))
    
    expect_true(is.data.frame(bias_rank$sample_id$Table))
    expect_true("Plot" %in% names(bias_rank$sample_id))
    
    expect_true(is.data.frame(bias_both$sample_id$Table))
    expect_true("Plot" %in% names(bias_both$sample_id))
})

test_that("The data frames for all example results have correct column(s)", {
    expect_true("rank_outlier" %in% colnames(bias_both$sample_id$Table) & 
                "dev_outlier" %in% colnames(bias_both$sample_id$Table))
    expect_true(!("rank_outlier" %in% colnames(bias_dev$sample_id$Table)) & 
                "dev_outlier" %in% colnames(bias_dev$sample_id$Table))
    expect_true("rank_outlier" %in% colnames(bias_rank$sample_id$Table) & 
                !("dev_outlier" %in% colnames(bias_rank$sample_id$Table)))
})

