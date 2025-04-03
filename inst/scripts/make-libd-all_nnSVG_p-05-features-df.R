library("spatialLIBD")
# spatialLIBD source
# Pardo B, Spangler A, Weber LM, Hicks SC, Jaffe AE, Martinowich K, Maynard KR, 
# Collado-Torres L (2022). “spatialLIBD: an R/Bioconductor package to visualize 
# spatially-resolved transcriptomics data.” BMC Genomics. 
# doi:10.1186/s12864-022-08601-w, https://doi.org/10.1186/s12864-022-08601-w.
# 
# Maynard KR, Collado-Torres L, Weber LM, Uytingco C, Barry BK, Williams SR, 
# II JLC, Tran MN, Besich Z, Tippani M, Chew J, Yin Y, Kleinman JE, Hyde TM, 
# Rao N, Hicks SC, Martinowich K, Jaffe AE (2021). “Transcriptome-scale spatial 
# gene expression in the human dorsolateral prefrontal cortex.” Nature 
# Neuroscience. doi:10.1038/s41593-020-00787-0, 
# https://www.nature.com/articles/s41593-020-00787-0.
# 
# Huuki-Myers LA, Spangler A, Eagles NJ, Montgomergy KD, Kwon SH, Guo B, 
# Grant-Peters M, Divecha HR, Tippani M, Sriworarat C, Nguyen AB, 
# Ravichandran P, Tran MN, Seyedian A, Consortium P, Hyde TM, Kleinman JE, 
# Battle A, Page SC, Ryten M, Hicks SC, Martinowich K, Collado-Torres L, 
# Maynard KR (2024). “A data-driven single-cell and spatial transcriptomic 
# map of the human prefrontal cortex.” Science. doi:10.1126/science.adh1938, 
# https://doi.org/10.1126/science.adh1938.
# 
# Kwon SH, Parthiban S, Tippani M, Divecha HR, Eagles NJ, Lobana JS, 
# Williams SR, Mark M, Bharadwaj RA, Kleinman JE, Hyde TM, Page SC, Hicks SC, 
# Martinowich K, Maynard KR, Collado-Torres L (2023). “Influence of 
# Alzheimer’s disease related neuropathology on local microenvironment gene 
# expression in the human inferior temporal cortex.” GEN Biotechnology. 
# doi:10.1089/genbio.2023.0019, https://doi.org/10.1089/genbio.2023.0019.

library("nnSVG")
# nnSVG source
# Weber L.M. et al. (2023), "nnSVG for the scalable identification of spatially 
# variable genes using nearest-neighbor Gaussian processes", Nature 
# Communications, 14, 4059

library("here")

spatialLIBD_spe <- fetch_data(type = "spe")
# clean data processed by LIBD collaborators

set.seed(123)
nnsvg <- nnSVG(spatialLIBD_spe)

write.csv(nnsvg, here("inst/extdata/libd-all_nnSVG_p-05-features-df.csv"))