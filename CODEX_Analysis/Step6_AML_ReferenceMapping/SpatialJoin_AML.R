# spatial join of qupath data to mesmer

library(sf)
library(readr)
library(nngeo)

# read in seurat metadata
AML1_183_md <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Spatial_Join_QuPath/seurat_metadata/AML1_183_md.csv")
AML1_382_md <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Spatial_Join_QuPath/seurat_metadata/AML1_382_md.csv")
AML2_191_md <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Spatial_Join_QuPath/seurat_metadata/AML2_191_md.csv")
AML2_380_md <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Spatial_Join_QuPath/seurat_metadata/AML2_380_md.csv")
AML3_1329_md <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Spatial_Join_QuPath/seurat_metadata/AML3_1329_md.csv")
AML3_1443_md <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Spatial_Join_QuPath/seurat_metadata/AML3_1443_md.csv")

# read in qupath export
qp_191 <- read_tsv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/AML_CODEX_Analysis_NBMProject/191_AllData_QuPath_bone_dists.txt")
qp_382 <- read_tsv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Spatial_Join_QuPath/382_AllData_QuPath_bone_dists.txt")
qp_183 <- read_tsv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Spatial_Join_QuPath/183_AllData_QuPath_bone_dists.txt")
qp_380 <- 
qp_1329 <- read_tsv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/AML_CODEX_Analysis_NBMProject/1329_AllData_QuPath_bone_dists.txt")
qp_1443 <-read_tsv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/AML_CODEX_Analysis_NBMProject/1443_AllData_QuPath_bone_dists.txt")

  
# create function for spatial join
  
qp_spatial_join <- function(qupath_md,mesmer_md) {
  qupath_md$`Centroid X µm` <- qupath_md$`Centroid X µm` / 0.50694
  qupath_md$`Centroid Y µm` <- qupath_md$`Centroid Y µm` / 0.50694
  # Make df & df1 sf objects, and keep the coordinates as columns just in case.
  df <-  qupath_md %>% st_as_sf(coords = c("Centroid X µm", "Centroid Y µm"), remove = FALSE) %>%
  st_set_crs(2193)
  df <- df %>% dplyr::filter(ROI  == "Geometry")
  df1 <- mesmer_md %>% st_as_sf(coords = c("x.coord", "y.coord"), remove = FALSE) %>%
  st_set_crs(2193)
  # crop dfs to be the same
  df <- df %>% dplyr::filter(`Centroid X µm` < max(df1$x.coord) & `Centroid Y µm` < max(df1$y.coord) & `Centroid X µm` > min(df1$x.coord) & `Centroid Y µm` > min(df1$y.coord))
  # Join df with df1, based on the nearest feature:
  df_near <- st_join(df1, df, join = st_is_within_distance, dist = 1)
  # write dataframe
  df_classified <- df_near %>% dplyr::select(colnames(mesmer_md), Class,  `Distance to annotation with Bone µm`)
  df_classified %>% ggplot(aes(x = Class, y = NPM1C, fill = Class)) + geom_violin() + theme_classic() + 
    scale_fill_manual(values = c("red", "steelblue")) + 
    ggtitle("AML") + 
    theme(plot.title = element_text(hjust = 0.5))
  p1 <- ggplot() + 
    geom_sf(data = df, color = 'red', size = 3) + 
    geom_sf(data = df1, color = 'black', alpha = .05) + 
    theme_void()
  ggsave(p1, paste0("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Spatial_Join_QuPath/output/",unique(mesmer_md$orig.ident),".pdf"))
  write.csv(df_classified, paste0("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Spatial_Join_QuPath/output/df_classified_bonedists", unique(mesmer_md$orig.ident), ".csv"))
  return(df_classified)
}

qp_spatial_join(qupath_md = qp_183, mesmer_md = AML1_183_md)
qp_spatial_join(qupath_md = qp_382, mesmer_md = AML1_382_md)
qp_spatial_join(qupath_md = qp_1329, mesmer_md = AML3_1329_md)
qp_spatial_join(qupath_md = qp_1443, mesmer_md = AML3_1443_md)
qp_spatial_join(qupath_md = qp_191, mesmer_md = AML2_191_md)


