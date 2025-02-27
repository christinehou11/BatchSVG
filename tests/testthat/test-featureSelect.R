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

test_that("The example result has correct class", {
    expect_true(is.list(list_batch_df))
})

test_that("The example result has corrent length with correct name", {
    expect_true(length(list_batch_df) == 2)
    expect_true(names(list_batch_df)[1] == "sample_id")
})

test_that("The function generates data frame including correct column for
            every batch effect", {
    expect_true(is.data.frame(list_batch_df$sex))
    expect_true("nSD_dev_sex" %in% colnames(list_batch_df$sex))
    expect_true("nSD_rank_sample_id" %in% colnames(list_batch_df$sample_id))
})