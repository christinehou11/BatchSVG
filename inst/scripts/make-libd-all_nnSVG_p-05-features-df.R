library("spatialLIBD")
library("here")
library("nnSVG")

spatialLIBD_spe <- fetch_data(type = "spe")
nnsvg <- nnSVG(spatialLIBD_spe)

write.csv(nnsvg, "~/Desktop/inst/extdata/libd-all_nnSVG_p-05-features-df.csv")