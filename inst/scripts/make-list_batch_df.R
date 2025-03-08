# run locally
library(BatchSVG)
library(spatialLIBD)

spatialLIBD_spe <- fetch_data(type = "spe")
libd_nnsvgs <- read.csv(
    system.file("extdata","libd-all_nnSVG_p-05-features-df.csv",
        package = "BatchSVG"),
    row.names = 1, check.names = FALSE)

list_batch_df <- featureSelect(spatiaLIBD_spe, batch_effects = "subject",
    VGs = libd_nnsvgs$gene_id)

save(list_batch_df, file = "~/Desktop/BatchSVG/inst/extdata/list_batch_df.rda")