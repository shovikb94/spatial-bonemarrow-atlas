### This script takes as input the Seurat Objects from each CODEX sample in the atlas and creates a combined CODEX Atlas Object

library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)

# Import RDS Files 
codex_H10 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H10_SeuratObj.RDS")
codex_H14 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H14_SeuratObj.RDS")
codex_H26 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H26_SeuratObj.RDS")
codex_H27 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H27_SeuratObj.RDS")
codex_H32 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H32_SeuratObj.RDS")
codex_H33 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H33_SeuratObj.RDS")
codex_H35 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H35_SeuratObj.RDS")
codex_H36 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H36_SeuratObj.RDS")
codex_H37 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H37_SeuratObj.RDS")
codex_H38 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H38_SeuratObj.RDS")
codex_H39 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H39_SeuratObj.RDS")
codex_H41 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Create_Seurat_Objs_Step1/Objects/codex_H41_SeuratObj.RDS")

ob.list <- list(codex_H10, codex_H14, codex_H26, codex_H27, codex_H32, codex_H33, codex_H35, codex_H36, codex_H37, codex_H38, codex_H39, codex_H41)
for (i in 1:length(ob.list)) {
  # ob.list[[i]] <- subset(ob.list[[i]], cells = cells.use)
  ob.list[[i]] <- NormalizeData(object = ob.list[[i]], normalization.method = "CLR", margin = 1)
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  all_features = rownames(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], features = all_features, approx = FALSE)  
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:30, reduction.name = "UMAP_dim30", reduction.key = "UMAP_dim30_")
  ob.list[[i]] <- FindNeighbors(ob.list[[i]], features = all_features, k.param=30) #features = all_features
  ob.list[[i]] <- FindClusters(ob.list[[i]], algorithm = 2, resolution = 1.5) # Louvain
  print('finished processing {i}th samples')
}

# ---- Work with Combined Data
codex_H10 <- ob.list[[1]]
codex_H14 <- ob.list[[2]]
codex_H26 <- ob.list[[3]]
codex_H27 <- ob.list[[4]]
codex_H32 <- ob.list[[5]]
codex_H33 <- ob.list[[6]]
codex_H35 <- ob.list[[7]]
codex_H36 <- ob.list[[8]]
codex_H37 <- ob.list[[9]]
codex_H38 <- ob.list[[10]]
codex_H39 <- ob.list[[11]]
codex_H41 <- ob.list[[12]]


combined <- merge(x = ob.list[[1]], y = c(ob.list[[2]], ob.list[[3]], ob.list[[4]], ob.list[[5]], ob.list[[6]], ob.list[[7]], ob.list[[8]],ob.list[[9]] ,ob.list[[10]], ob.list[[11]], ob.list[[12]]), add.cell.ids = c("H10", "H14","H26", "H27", "H32","H33", "H35","H36", "H37", "H38", "H39", "H41"))
immune.anchors <- FindIntegrationAnchors(object.list = ob.list, anchor.features = all_features, reduction = "rpca")
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
saveRDS(immune.combined, "CODEX_RPCA_Integrated.RDS")       