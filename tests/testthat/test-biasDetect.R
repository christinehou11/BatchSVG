library(spatialLIBD)

spatialLIBD_spe <- fetch_data(type = "spe")
libd_svg <- read.csv(
    system.file("extdata","libd-all_nnSVG_p-05-features-df.csv",
        package = "BatchSVG"),
    row.names = 1, check.names = FALSE)

list_batch_df <- featureSelect(input = spatialLIBD_spe, 
    batch_effects = "subject", VGs = libd_svg$gene_id)

bias_dev <- biasDetect(list_batch_df = list_batch_df, threshold = "dev",
        nSD_dev = 3)

bias_rank <- biasDetect(list_batch_df = list_batch_df, threshold = "rank",
        nSD_rank = 3)

bias_both <- biasDetect(list_batch_df = list_batch_df, threshold = "both",
        nSD_dev = 3, nSD_rank = 3)

test_that("Both example results have one batch effect" , {
    expect_true(length(bias_dev) == 1)
    expect_true(length(bias_rank) == 1)
    expect_true(length(bias_both) == 1)
})

test_that("All example results have table and plot saved for each batch", {
    expect_true(is.data.frame(bias_dev$subject$Table))
    expect_true("Plot" %in% names(bias_dev$subject))
})

test_that("The data frames for all example results have correct column(s)", {
    expect_true("rank_outlier" %in% colnames(bias_both$subject$Table) & 
                "dev_outlier" %in% colnames(bias_both$subject$Table))
})

