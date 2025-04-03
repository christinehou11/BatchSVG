suppressMessages(library(spatialLIBD))

spatialLIBD_spe <- fetch_data(type = "spe")
libd_svg <- read.csv(
    system.file("extdata","libd-all_nnSVG_p-05-features-df.csv",
        package = "BatchSVG"),
    row.names = 1, check.names = FALSE)

list_batch_df <- featureSelect(input = spatialLIBD_spe, 
    batch_effects = "subject", VGs = libd_svg$gene_id)

test_that("The example result has correct class", {
    expect_true(is.list(list_batch_df))
})

test_that("The example result has corrent length with correct name", {
    expect_true(length(list_batch_df) == 1)
    expect_true(names(list_batch_df)[1] == "subject")
})

test_that("The function generates data frame including correct column for
            every batch effect", {
    expect_true(is.data.frame(list_batch_df$subject))
    expect_true("nSD_dev_subject" %in% colnames(list_batch_df$subject))
    expect_true("nSD_rank_subject" %in% colnames(list_batch_df$subject))
})