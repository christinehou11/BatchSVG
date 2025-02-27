data("spe_sub4")
data("svgs_sub4")
SVGs <- svgs_sub4$gene_id

list_batch_df <- featureSelect(input = spe_sub4, 
    batch_effects = c("sample_id", "sex"), VGs = SVGs)

plots <- svg_nSD(list_batch_df = list_batch_df, sd_interval_dev = c(4,7),
    sd_interval_rank = c(6,7))

test_that("multiplication works", {
    expect_true(length(plots) == 2)
})