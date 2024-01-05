library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(scCustomize)
options(Seurat.object.assay.version = "v4")
setwd("Seurat_Analysis/Bone_Marrow")

# read in channelnames
channel_names <- read_csv(file = "FUSION/NBM/NBM27_H10_CITRATE_REIMAGE/H10/Scan1/.temp/MarkerList.txt", col_names = FALSE)
sample_names <- c('H10', 'H14', 'H26', 'H27', 'H32', 'H33', 'H35', 'H36', 'H37', 'H38', 'H39', 'H41')
first_file <- paste0('Seurat_Analysis/Bone_Marrow/Osteo-MSC_Masks/', sample_names[1], '_Quant/combined_markers.csv')
first_data <- read_csv(first_file)

#Specify the directory of your data
data_list1 <- list()
for (sample_name in sample_names) {
  file_path <- paste0('Seurat_Analysis/Bone_Marrow/Osteo-MSC_Masks/', sample_name, '_Quant/combined_markers.csv')
  gc_csd_raw <- read_csv(file_path)
  gc_csd_raw$indexes <- paste0(sample_name,"_Osteo_",gc_csd_raw$indexes)
  colnames(gc_csd_raw) <- colnames(first_data)
  data_list1[[sample_name]] <- gc_csd_raw
}
gc_csd_raw1 <- do.call(rbind, data_list1)
data_list2 <- list()
for (sample_name in sample_names) {
  file_path <- paste0('Seurat_Analysis/Bone_Marrow/Fibro-MSC_Masks/', sample_name, '_Quant/combined_markers.csv')
  gc_csd_raw <- read_csv(file_path)
  gc_csd_raw$indexes <- paste0(sample_name,"_Fibro_",gc_csd_raw$indexes)
  colnames(gc_csd_raw) <- colnames(first_data)
  data_list2[[sample_name]] <- gc_csd_raw
}
gc_csd_raw2 <- do.call(rbind, data_list2)
gc_csd_raw <- rbind(gc_csd_raw1, gc_csd_raw2)

colnames(gc_csd_raw)[5:58] <- channel_names$X1
colnames(gc_csd_raw)[2] <- "Size"
colnames(gc_csd_raw)[1] <- "CellID"
colnames(gc_csd_raw)[3:4] <- c("x.coord", "y.coord")
gc_csd_raw[,3:58] <- gc_csd_raw[,3:58] / gc_csd_raw$Size 

gc_csd_raw <- gc_csd_raw %>% dplyr::filter(`Size`<1000000)

# Specify markers to be excluded from clustering, NOTE that if the name of one marker is contained in another this will remove both!
failed <- c("DAPI", "NAKATPASE", "Ki67", "ADIPOQ", "CD56")

# Remove the size and x/y coords
gc_csd2 <- gc_csd_raw %>% select(-starts_with('x.coord'), -contains("Empty"), -starts_with('y.coord'),-contains('size'),-contains("Std.Dev"), -contains("indexes"),-contains('DAPI'), -contains('RowSum'), -contains('Blank'),-starts_with("Ch"), -starts_with("faile"), -contains(failed)) # Note removed failed markers
gc_csd2 <- gc_csd2 %>% select(-contains("Index"),-contains("CellID"), -matches("Python_Index"))
#Rename objects
gc_csd2 = gc_csd2 %>%
  rename_with(~str_remove_all(.x, 'Nucleus Intensity')) %>%
  rename_with(~str_remove_all(.x, '\\(.*\\) ')) %>%
  rename_with(~str_remove_all(.x, ' Mean'))
gc_csd2 <- na.omit(gc_csd2)

# Create a Seurat object
t_gc_csd = t(as.matrix(gc_csd2))
colnames(t_gc_csd) <- rownames(gc_csd2)
rownames(t_gc_csd) <- colnames(gc_csd2)
codex <- CreateAssayObject(counts = t_gc_csd, project = "Osteo_Fibro_MSC_Manual", assay="CODEX")
codex_seurat <- CreateSeuratObject(codex, assay = "CODEX", project = "Osteo_Fibro_MSC_Manual", meta.data = gc_csd_raw)
extract_parts <- lapply(gc_csd_raw$CellID, function(text) {
  parts <- strsplit(text, "_")[[1]]
  part1 <- parts[1]  
  part2 <- parts[2]  
  return(list(part1 = part1, part2 = part2))
})
sample_names <- sapply(extract_parts, function(x) x$part1)
cell_names <- sapply(extract_parts, function(x) x$part2)
cell_names <- paste0(cell_names, "_MSC")
codex_seurat$orig.ident <- sample_names
codex_seurat$cluster_anno_l2 <- cell_names
saveRDS(codex_seurat, "Osteo_Fibro_MSC_Manual_SeuratObj.RDS")

#Merge with existing Seurat object
immune.combined <- readRDS("immune.combined_final_040523.RDS")
added_cells <- readRDS("Osteo_Fibro_MSC_Manual_SeuratObj.RDS")
dim(immune.combined)
dim(added_cells)
metadata_names_before <- colnames(immune.combined@meta.data)
print(metadata_names_before)
colnames(immune.combined@meta.data) <- gsub("^Metadata_", "", colnames(immune.combined@meta.data))
metadata_names_after <- colnames(immune.combined@meta.data)
print(metadata_names_after)

merged_seurat <- merge(immune.combined, codex_seurat)
dim(merged_seurat)

new_assay <- merged_seurat@meta.data[,9:61]
new_assay$NAKATPASE <- NULL
new_assay$Ki67 <- NULL
new_assay$ADIPOQ <- NULL
dim(new_assay)
dim(merged_seurat)
new_assay <- t(new_assay)

merged_seurat[["CODEX2"]] <- CreateAssayObject(counts = new_assay)
DefaultAssay(merged_seurat) <- "CODEX2"
merged_seurat[["CODEX"]] <- NULL
merged_seurat[["integrated"]] <- NULL

table(merged_seurat$cluster_anno_l2)

#Renormalize merged data with CD56 now added
merged_seurat <- NormalizeData(object = merged_seurat, normalization.method = "CLR", margin = 1)
saveRDS(merged_seurat, "Immune_Combined_Added_MSC_CD56.RDS")

#Do some plots 
DefaultAssay(merged_seurat)
Stacked_VlnPlot(merged_seurat, features = c("CD45", "codex2_FOXC1", "codex2_CXCL12", "codex2_CD90", "codex2_CD56", "codex2_PDPN"), 
                group.by = "cluster_anno_l2", x_lab_rotate = TRUE, pt.size = 0, raster = F)

p1 <- VlnPlot(merged_seurat, "codex2_CD45", group.by = "cluster_anno_l2", raster = F, pt.size = 0) + theme(axis.text.x = element_blank(), legend.position = "none") 
p2 <- VlnPlot(merged_seurat, "codex2_FOXC1", group.by = "cluster_anno_l2", raster = F, pt.size = 0) + theme(axis.text.x = element_blank(), legend.position = "none") 
p3 <- VlnPlot(merged_seurat, "codex2_CXCL12", group.by = "cluster_anno_l2", raster = F, pt.size = 0) + theme(axis.text.x = element_blank(), legend.position = "none") 
p4 <- VlnPlot(merged_seurat, "codex2_CD90", group.by = "cluster_anno_l2", raster = F, pt.size = 0) + theme(axis.text.x = element_blank(), legend.position = "none") 
p5 <- VlnPlot(merged_seurat, "codex2_CD56", group.by = "cluster_anno_l2", raster = F, pt.size = 0) + theme(axis.text.x = element_blank(), legend.position = "none") 
p6 <- VlnPlot(merged_seurat, "codex2_PDPN", group.by = "cluster_anno_l2", raster = F, pt.size = 0) + theme(legend.position = "none")
p1+p2+p3+p4+p5+p6
