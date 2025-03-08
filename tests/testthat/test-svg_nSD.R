# library(spatialLIBD)
# 
# spatialLIBD_spe <- fetch_data(type = "spe")
# libd_svg <- read.csv(
#     system.file("extdata","libd-all_nnSVG_p-05-features-df.csv",
#         package = "BatchSVG"),
#     row.names = 1, check.names = FALSE)
# 
# list_batch_df <- featureSelect(input = spatialLIBD_spe, 
#     batch_effects = "subject", VGs = libd_svg$gene_id)

load(system.file("extdata","list_batch_df.rda", package = "BatchSVG"))

plots <- svg_nSD(list_batch_df = list_batch_df, sd_interval_dev = 3,
    sd_interval_rank = 3)

test_that("The example results have the correct number of batch effects", {
    expect_true(length(plots) == 1)
})