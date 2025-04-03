#' Generate example list_batch_df for BatchSVG
#'
#' This script generates `list_batch_df` using spatialLIBD data
#' and the BatchSVG package. It is stored in `inst/extdata/` for examples
#' and testing purposes.
#'
#' Run locally before package submission.
#' 
#' Keep here for testing and example convenience. Use ?list_batch_df to read
#' information page.
#' 
library("BatchSVG")
library("spatialLIBD")
library("here")

spatialLIBD_spe <- fetch_data(type = "spe")
libd_nnsvgs <- read.csv(
    system.file("extdata","libd-all_nnSVG_p-05-features-df.csv",
        package = "BatchSVG"),
    row.names = 1, check.names = FALSE)

list_batch_df <- featureSelect(spatiaLIBD_spe, batch_effects = "subject",
    VGs = libd_nnsvgs$gene_id, verbose = FALSE)

save(list_batch_df, file = here("inst/extdata/list_batch_df.rda"))