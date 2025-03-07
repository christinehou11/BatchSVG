# Run this on JHPCE
# /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/nnSVG
#
# Article citation:
# Jacqueline R. Thompson, Erik D. Nelson, Madhavi Tippani, Anthony D. Ramnauth, 
# Heena R. Divecha, Ryan A. Miller, Nicholas J. Eagles, Elizabeth A. Pattie, 
# Sang Ho Kwon, Svitlana V. Bach, Uma M. Kaipa, Jianing Yao, Christine Hou, 
# Joel E. Kleinman, Leonardo Collado-Torres, Shizhong Han, Kristen R. Maynard,
# Thomas M. Hyde, Keri Martinowich, Stephanie C. Page, and Stephanie C. Hicks.
# An integrated single-nucleus and spatial transcriptomics atlas reveals the 
# molecular landscape of the human hippocampus. biorxiv. 2025. 
# DOI: 10.1101/2024.04.26.590643

library('here')
library('tidyverse')

# load nnSVG results
load(here("nnSVG_outs_HE_only.rda"))

res_df_sub <- pivot_longer(
    rownames_to_column(as.data.frame(res_ranks), var<-"gene_id"), 
    colnames(res_ranks), 
    names_to="sample_id", 
    values_to="rank", 
    values_drop_na=TRUE)

res_df_sub <- filter(res_df_sub,
    sample_id %in% 
        c("V11L05-333_B1", "V11L05-333_D1", "V11L05-335_D1", "V11L05-336_A1"), 
    rank <= 2000) # top 2k sig features

svgs_sub4 <- group_by(res_df_sub, gene_id) |>
    tally() |> 
    filter(n>1)

write.csv(svgs_sub4, "~/Desktop/BatchSVG/inst/extdata/BatchSVG/svgs_sub4.csv")