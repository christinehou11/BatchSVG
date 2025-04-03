data(list_batch_df)

plots <- svg_nSD(list_batch_df = list_batch_df, sd_interval_dev = 3,
    sd_interval_rank = 3)

test_that("The example results have the correct number of batch effects", {
    expect_true(length(plots) == 1)
})