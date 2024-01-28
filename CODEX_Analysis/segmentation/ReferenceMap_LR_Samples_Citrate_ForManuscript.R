library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggrastr)
setwd("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand")

immune.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.combined_final_040523.RDS")

DefaultAssay(immune.combined) <- "integrated"

# read in LR objects
NBM67_H38_Citrate <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Seurat_Objects/NBM67_H38_Citrate_CODEX_SeuratObj.RDS")
NBM67_H41_Citrate <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Seurat_Objects/NBM67_H41_Citrate_CODEX_SeuratObj.RDS")
NBM67_H14_Citrate <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Seurat_Objects/NBM67_H14_Citrate_CODEX_SeuratObj.RDS")
NBM70_H10_Citrate <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Seurat_Objects/NBM70_H10_Citrate_CODEX_SeuratObj.RDS")
NBM70_H26_Citrate <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Seurat_Objects/NBM70_H26_Citrate_CODEX_SeuratObj.RDS")
NBM70_H32_Citrate <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Seurat_Objects/NBM70_H32_Citrate_CODEX_SeuratObj.RDS")

ob.list <- list(NBM67_H38_Citrate,NBM67_H41_Citrate,NBM67_H14_Citrate, 
NBM70_H10_Citrate,NBM70_H26_Citrate,NBM70_H32_Citrate)

for (i in 1:length(ob.list)) {
  print(i)
  ob.list[[i]] <- NormalizeData(object = ob.list[[i]], normalization.method = "CLR", margin = 1)
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  all_features = rownames(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], features = all_features, approx = FALSE)  
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:30, reduction.name = "UMAP_dim30", reduction.key = "UMAP_dim30_")
  ob.list[[i]] <- FindNeighbors(ob.list[[i]], features = all_features, k.param=30) #features = all_features
  ob.list[[i]] <- FindClusters(ob.list[[i]], algorithm = 2, resolution = 1.5) # Louvain
}

# ---- Work with Combined Data
NBM67_H38_Citrate_processed <- ob.list[[1]]
NBM67_H41_Citrate_processed <- ob.list[[2]]
NBM67_H14_Citrate_processed  <- ob.list[[3]]
NBM70_H10_Citrate_processed <- ob.list[[4]]
NBM70_H26_Citrate_processed <- ob.list[[5]]
NBM70_H32_Citrate_processed <- ob.list[[6]]

LR.combined_citrate <- merge(x=NBM67_H38_Citrate_processed, y = c(NBM67_H41_Citrate_processed,
                                                          NBM67_H14_Citrate_processed, 
                                                          NBM70_H10_Citrate_processed,
                                                          NBM70_H26_Citrate_processed,
                                                          NBM70_H32_Citrate_processed), add.cell.ids = c("NBM67_H38_Citrate", 
                                                                                                         "NBM67_H41_Citrate", 
                                                                                                         "NBM67_H14_Citrate",
                                                                                                         "NBM70_H10_Citrate", 
                                                                                                         "NBM70_H26_Citrate", 
                                                                                                         "NBM70_H32_Citrate"))


saveRDS(LR.combined_citrate, "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/LR.combined.mapped_unfiltered_citrate1.RDS")

anchors_citrate <- FindTransferAnchors(
  reference = immune.combined, k.filter = NA,
  query = LR.combined_citrate,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:30
)

saveRDS(anchors_citrate, "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/anchors_citrate.RDS")

Key(object = immune.combined[["UMAP_dim30"]]) <- "umapdim30_"

LR.combined_citrate_mapped <- MapQuery(
  anchorset = anchors_citrate,
  query = LR.combined_citrate,
  reference = immune.combined,
  refdata = list(
    refmap_fine = "cluster_anno_l2", 
    refmap_coarse = "cluster_anno_coarse"
  ),
  reference.reduction = "pca", 
  reduction.model = "UMAP_dim30"
)


LR.combined_citrate_mapped$classified_cluster_anno_l2 <- LR.combined_citrate_mapped$predicted.refmap_fine
LR.combined_citrate_mapped$classified_cluster_anno_l2_score <- LR.combined_citrate_mapped$predicted.refmap_fine.score
LR.combined_citrate_mapped$classified_cluster_anno_coarse <- LR.combined_citrate_mapped$predicted.refmap_coarse
LR.combined_citrate_mapped$classified_cluster_anno_coarse_score <- LR.combined_citrate_mapped$predicted.refmap_coarse.score

saveRDS(LR.combined_citrate_mapped, "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/LR.combined.mapped_unfiltered_citrate2.RDS")

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(as.factor(LR.combined_citrate_mapped$classified_cluster_anno_l2))))

DimPlot(LR.combined_citrate_mapped, reduction = "ref.umap", group.by = "classified_cluster_anno_l2", raster = F, 
        label = TRUE, label.size = 3, repel = TRUE, cols = cols) + NoLegend() + coord_fixed() + NoAxes()
VlnPlot(LR.combined_citrate_mapped, features = "classified_cluster_anno_l2_score", raster = F, pt.size = 0, cols = cols, 
        group.by = "classified_cluster_anno_l2")

ggplot(LR.combined_citrate_mapped@meta.data, aes(fill=classified_cluster_anno_l2, x=orig.ident)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  scale_fill_manual(values = cols) + 
  ylab("Percentage of Total Cells") + xlab("Sample") + guides(fill="none")

LR.combined_citrate_mapped@meta.data["Sample_Name"] <- LR.combined_citrate_mapped$orig.ident

#Generate masks for Citrate
annotations_col <- "classified_cluster_anno_l2"
samples_citrate <- unique(LR.combined_citrate_mapped$orig.ident)
cluster_factor <- as.numeric(as.factor(LR.combined_citrate_mapped$classified_cluster_anno_l2))
LR.combined_citrate_mapped$cluster_levels <- cluster_factor
conversion_table_citrate <- data.frame(
  Level = levels(as.factor(LR.combined_citrate_mapped$classified_cluster_anno_l2)),
  Numeric = levels(factor(as.numeric(factor(cluster_factor))))
)
conversion_table_citrate

for(name in samples_citrate){
  print(name)
  print("Subseting")
  sample_subset <- subset(LR.combined_citrate_mapped, subset=orig.ident==name)
  annotations_df <- as.data.frame(sample_subset$cluster_levels)
  rownames(annotations_df) <- sample_subset$CellID
  annotations_df$CellID <- sample_subset$CellID
  colnames(annotations_df) <- c("Annotation", "CellID")
  print("Writing Annotations")
  write.csv(annotations_df, file = paste0("Reference_Map_Masks/",name,"_ReferenceMap.csv"))
}
