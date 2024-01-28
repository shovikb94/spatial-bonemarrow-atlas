### This script takes as input the Seurat Objects from each CODEX sample in the atlas and creates a combined CODEX Atlas Object

library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(MetBrewer)
library(cowplot)
library(circlize)


immune.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Integrate_Samples_RPCA_Step2/objects/CODEX_RPCA_Integrated.RDS")

# Output per sample metadata for checking the cluster mask

# H10_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM27_H10_CODEX_Mesmer')
# H14_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM28_H14_CODEX_Mesmer')
# H26_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM36_H26_CODEX_Mesmer')
# H27_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM41_H27_CODEX_Mesmer')
# H32_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM31_H32_CODEX_Mesmer')
# H33_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM38_H33_CODEX_Mesmer')
# H35_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM37_H35_CODEX_Mesmer')
# H36_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM33_H36_CODEX_Mesmer')
# H37_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM32_H37_CODEX_Mesmer')
# H38_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM34_H38_CODEX_Mesmer')
# H39_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM40_H39_CODEX_Mesmer')
# H41_md <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM39_H41_CODEX_Mesmer')
# 
# H10_md$cluster_anno_l2_n <- as.numeric(factor(H10_md$cluster_anno_l2))
# H14_md$cluster_anno_l2_n <- as.numeric(factor(H14_md$cluster_anno_l2))
# H26_md$cluster_anno_l2_n <- as.numeric(factor(H26_md$cluster_anno_l2))
# H27_md$cluster_anno_l2_n <- as.numeric(factor(H27_md$cluster_anno_l2))
# H32_md$cluster_anno_l2_n <- as.numeric(factor(H32_md$cluster_anno_l2))
# H33_md$cluster_anno_l2_n <- as.numeric(factor(H33_md$cluster_anno_l2))
# H35_md$cluster_anno_l2_n <- as.numeric(factor(H35_md$cluster_anno_l2))
# H36_md$cluster_anno_l2_n <- as.numeric(factor(H36_md$cluster_anno_l2))
# H37_md$cluster_anno_l2_n <- as.numeric(factor(H37_md$cluster_anno_l2))
# H38_md$cluster_anno_l2_n <- as.numeric(factor(H38_md$cluster_anno_l2))
# H39_md$cluster_anno_l2_n <- as.numeric(factor(H39_md$cluster_anno_l2))
# H41_md$cluster_anno_l2_n <- as.numeric(factor(H41_md$cluster_anno_l2))
# 
# write.csv(x= data.frame(H10_md$x.coord, H10_md$y.coord, H10_md$CellID, H10_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H10_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H14_md$x.coord, H14_md$y.coord, H14_md$CellID, H14_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H14_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H26_md$x.coord, H26_md$y.coord, H26_md$CellID, H26_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H26_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H27_md$x.coord, H27_md$y.coord, H27_md$CellID, H27_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H27_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H32_md$x.coord, H32_md$y.coord, H32_md$CellID, H32_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H32_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H33_md$x.coord, H33_md$y.coord, H33_md$CellID, H33_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H33_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H35_md$x.coord, H35_md$y.coord, H35_md$CellID, H35_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H35_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H36_md$x.coord, H36_md$y.coord, H36_md$CellID, H36_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H36_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H37_md$x.coord, H37_md$y.coord, H37_md$CellID, H37_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H37_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H38_md$x.coord, H38_md$y.coord, H38_md$CellID, H38_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H38_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H39_md$x.coord, H39_md$y.coord, H39_md$CellID, H39_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H39_mesmer_cluster_anno_l2_annotations_110322.csv")
# write.csv(x= data.frame(H41_md$x.coord, H41_md$y.coord, H41_md$CellID, H41_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H41_mesmer_cluster_anno_l2_annotations_110322.csv")
# 



# rename metadata
renamed_md_vect <- c(paste0("Metadata_",rownames(immune.combined), sep = ""))
colnames(immune.combined@meta.data)[35] <- "HLA-DR"
colnames(immune.combined@meta.data)[colnames(immune.combined@meta.data) %in% rownames(immune.combined)] <- renamed_md_vect

# Rename clusters based on marker expression from ridge plots
new.cluster.ids <- c("Erythroid", "Myeloid", "Myeloid", "Macrophages", "Erythroid", "B-Cells",
                     "Myeloid", "Monocytes", "CD8+ T-Cell", "Autofluorescent", "CD4+ T-Cell", "Stromal", "Plasma Cells",
                     "Endothelial", "Megakaryocyte", "Erythroid", "HSPC", "pDC", "VSMC", "Undetermined", "CD79A+ CD61+ Undetermined", "CD44+ Undetermined", "Undetermined", "Undetermined")


immune.combined <- SetIdent(immune.combined, value = "seurat_clusters")
names(new.cluster.ids) <- levels(immune.combined) 
immune.combined<- RenameIdents(immune.combined, new.cluster.ids)
immune.combined@meta.data["cluster_anno_l1"] <- immune.combined@active.ident  


# See if the other clusters are better resolved if filtering out the confident clusters
all_features <- rownames(immune.combined)
other.combined <- subset(immune.combined, subset = seurat_clusters %in% c(11,13, 16,18,19,20,22,23))
other.combined <- ScaleData(other.combined)
other.combined <- RunPCA(other.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
other.combined <- RunUMAP(other.combined, dims = 1:30)
other.combined <- FindNeighbors(other.combined, features = all_features, k.param=30)#features = all_features,
# other.combined <- FindClusters(other.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
other.combined <- FindClusters(other.combined, algorithm = 2, resolution = 1) # Louvain

library(RColorBrewer)
getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(other.combined@active.ident)))
DimPlot(other.combined, reduction = "umap", cols = cols, pt.size = 0.1, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()


# Change IDs based on reclustering results
immune.combined$cluster_anno_l2 <- immune.combined$cluster_anno_l1
AEC_ids <- colnames(subset(other.combined, subset = seurat_clusters == 0 | seurat_clusters == 14))
MSC_ids <- colnames(subset(other.combined, subset = seurat_clusters == 5))
SEC_ids <- colnames(subset(other.combined, subset = seurat_clusters == 3 | seurat_clusters == 12))
HSPC_ids <- colnames(subset(other.combined, subset = seurat_clusters == 9 | seurat_clusters == 7 | seurat_clusters == 15))
VSMC_ids <- colnames(subset(other.combined, subset = seurat_clusters == 16))


immune.combined$cluster_anno_l2 <- as.character(immune.combined$cluster_anno_l2)
immune.combined$cluster_anno_l2[AEC_ids] <- "AEC"
immune.combined$cluster_anno_l2[MSC_ids] <- "MSC"
immune.combined$cluster_anno_l2[SEC_ids] <- "SEC"
immune.combined$cluster_anno_l2[HSPC_ids] <- "HSPC"
immune.combined$cluster_anno_l2[VSMC_ids] <- "VSMC"

# Check the erythroid cells and find the erythroblasts

erythroid.combined <- subset(immune.combined, subset = cluster_anno_l1 == 'Erythroid' )
erythroid.combined <- ScaleData(erythroid.combined)
erythroid.combined <- RunPCA(erythroid.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
erythroid.combined <- RunUMAP(erythroid.combined, dims = 1:30)
erythroid.combined <- FindNeighbors(erythroid.combined, features = all_features, k.param=30)#features = all_features,
# erythroid.combined <- FindClusters(erythroid.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
erythroid.combined <- FindClusters(erythroid.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(erythroid.combined@active.ident)))
DimPlot(erythroid.combined, reduction = "umap", cols = cols, pt.size = 0.1, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

# Rename Clusters Based On Expression
Erythroblast_ids <- colnames(subset(erythroid.combined, subset = seurat_clusters == 9)) # GATA1 and CD71 expression
immune.combined$cluster_anno_l2[Erythroblast_ids] <- "Erythroblast"
EarlyMyeloid_Ery_correction_ids <- colnames(subset(erythroid.combined, subset = seurat_clusters == 5 | seurat_clusters == 8)) # MPO expression / lack of CD11B
immune.combined$cluster_anno_l2[EarlyMyeloid_Ery_correction_ids] <- "Early Myeloid Progenitor"
MatureMyeloid_Ery_correction_ids <- colnames(subset(erythroid.combined, subset = seurat_clusters == 3)) # CD15 + CD11B, lack of MPO
immune.combined$cluster_anno_l2[MatureMyeloid_Ery_correction_ids] <- "Mature Myeloid"

# Repeat for B-Cells
B_Cell.combined <- subset(immune.combined, subset = cluster_anno_l1 == 'B-Cells')
B_Cell.combined <- ScaleData(B_Cell.combined)
B_Cell.combined <- RunPCA(B_Cell.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
B_Cell.combined <- RunUMAP(B_Cell.combined, dims = 1:30)
B_Cell.combined <- FindNeighbors(B_Cell.combined, features = all_features, k.param=30)#features = all_features,
# B_Cell.combined <- FindClusters(B_Cell.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
B_Cell.combined <- FindClusters(B_Cell.combined, algorithm = 2, resolution = 1) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(B_Cell.combined@active.ident)))

DimPlot(B_Cell.combined, reduction = "umap", cols = cols, pt.size = 0.1, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

# Rename Immature-B cluster
ImmatureB_ids <- colnames(subset(B_Cell.combined, subset = seurat_clusters == 0))
immune.combined$cluster_anno_l2[ImmatureB_ids] <- "Immature_B_Cell"
EarlyMyeloid_B_correction_ids <- colnames(subset(B_Cell.combined, subset = seurat_clusters == 15)) # MPO expression / lack of CD11B
immune.combined$cluster_anno_l2[EarlyMyeloid_B_correction_ids] <- "Early Myeloid Progenitor"
MatureMyeloid_B_correction_ids <- colnames(subset(B_Cell.combined, subset = seurat_clusters == 4 | seurat_clusters == 12)) # CD15 + CD11B, lack of MPO
immune.combined$cluster_anno_l2[MatureMyeloid_B_correction_ids] <- "Mature Myeloid"

# Repeat for HSPCs
HSPC.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'HSPC')
HSPC.combined <- ScaleData(HSPC.combined)
HSPC.combined <- RunPCA(HSPC.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
HSPC.combined <- RunUMAP(HSPC.combined, dims = 1:30)
HSPC.combined <- FindNeighbors(HSPC.combined, features = all_features, k.param=30)#features = all_features,
# HSPC.combined <- FindClusters(HSPC.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
HSPC.combined <- FindClusters(HSPC.combined, algorithm = 2, resolution = 1) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(HSPC.combined@active.ident)))
DimPlot(HSPC.combined, group.by = 'seurat_clusters', reduction = "umap", cols = cols, pt.size = 1, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()


# Rename
GMP_ids <- colnames(subset(HSPC.combined, subset = seurat_clusters == 0))
immune.combined$cluster_anno_l2[GMP_ids] <- "GMP/Myeloblast"
SPINK2posHSPC_ids <- colnames(subset(HSPC.combined, subset = seurat_clusters == 1 | seurat_clusters == 5 | seurat_clusters == 11))
immune.combined$cluster_anno_l2[SPINK2posHSPC_ids] <- "SPINK2+ HSPC"
EarlyMyeloid_HSPC_correction_ids <- colnames(subset(HSPC.combined, subset = seurat_clusters == 9)) # MPO expression / lack of CD11B / high VECAD suggesting these are perivascular MPO+ progenitors mislabeled as CD34+ due to lateral spillover
immune.combined$cluster_anno_l2[EarlyMyeloid_HSPC_correction_ids] <- "Early Myeloid Progenitor"
MatureMyeloid_HSPC_correction_ids <- colnames(subset(HSPC.combined, subset = seurat_clusters == 2 | seurat_clusters == 3)) # CD15 + CD11B, lack of MPO
immune.combined$cluster_anno_l2[MatureMyeloid_HSPC_correction_ids] <- "Mature Myeloid"
Megakaryoblast_ids <- colnames(subset(HSPC.combined, subset = CD61 > 0.8)) # CD34 and CD61, checking on image a mix of spillover and true double positives
immune.combined$cluster_anno_l2[Megakaryoblast_ids] <- "CD34+ CD61+"
MEP_ids <- colnames(subset(HSPC.combined, subset = GATA1 > 0.8)) # GATA1 and CD34 co-expression. 
immune.combined$cluster_anno_l2[MEP_ids] <- "MEP/Early Erythroblast"

# Repeat for Macrophages
Macrophage.combined <- subset(immune.combined, subset = cluster_anno_l1 == 'Macrophages')
Macrophage.combined <- ScaleData(Macrophage.combined)
Macrophage.combined <- RunPCA(Macrophage.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
Macrophage.combined <- RunUMAP(Macrophage.combined, dims = 1:30)
Macrophage.combined <- FindNeighbors(Macrophage.combined, features = all_features, k.param=30)#features = all_features,
# Macrophage.combined <- FindClusters(Macrophage.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
Macrophage.combined <- FindClusters(Macrophage.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(Macrophage.combined@active.ident)))
DimPlot(Macrophage.combined, reduction = "umap", cols = cols, pt.size = 0.1, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()


# Rename Macrophage Clusters

Macrophage_ids <- colnames(subset(Macrophage.combined, subset = seurat_clusters == 1 | seurat_clusters == 5 | seurat_clusters == 7)) # Expression of CD68/CD163. Noted on image that some false positives are in cluster 7, but mostly these seem to be Cd68/CD163++ Macrophages
immune.combined$cluster_anno_l2[Macrophage_ids] <- "Macrophages"
EarlyMyeloid_Macrophage_correction_ids <- colnames(subset(Macrophage.combined, subset = seurat_clusters == 6)) # MPO expression, CD33 expression, lower CD11B and CD68/CD163
immune.combined$cluster_anno_l2[EarlyMyeloid_Macrophage_correction_ids] <- "Early Myeloid Progenitor"
Mature_Myeloid_Macrophage_correction_ids <- colnames(subset(Macrophage.combined, subset = seurat_clusters == 0 | seurat_clusters == 2 | seurat_clusters == 3 | seurat_clusters == 4)) # CD11B and CD15 expression, lower CD68/CD163
immune.combined$cluster_anno_l2[Mature_Myeloid_Macrophage_correction_ids] <- "Mature Myeloid"

# Repeat for Endothelial 
Endothelial.combined <- subset(immune.combined, subset = cluster_anno_l1 == 'Endothelial')
Endothelial.combined <- ScaleData(Endothelial.combined)
Endothelial.combined <- RunPCA(Endothelial.combined, features = all_features, approx = FALSE)
Endothelial.combined <- RunUMAP(Endothelial.combined, dims = 1:30)
Endothelial.combined <- FindNeighbors(Endothelial.combined, features = all_features, k.param=30)#features = all_features,
Endothelial.combined <- FindClusters(Endothelial.combined, algorithm = 2, resolution = 0.5) # Louvain

library(RColorBrewer)
getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(Endothelial.combined@active.ident)))
DimPlot(Endothelial.combined, reduction = "umap", cols = cols, pt.size = 0.1, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

# Rename
AEC_ids <- colnames(subset(Endothelial.combined, subset = seurat_clusters == 2))
immune.combined$cluster_anno_l2[AEC_ids] <- "AEC"
SEC_ids <- colnames(subset(Endothelial.combined, subset = seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 3 | seurat_clusters == 4)) # Expression of VECAD, no CXCL12 
immune.combined$cluster_anno_l2[SEC_ids] <- "SEC"
EarlyMyeloid_Endothelial_correction_ids <- colnames(subset(Endothelial.combined, subset = seurat_clusters == 6)) # MPO expression, CD33 expression, lower CD34 / VECAD 
immune.combined$cluster_anno_l2[EarlyMyeloid_Endothelial_correction_ids] <- "Early Myeloid Progenitor"
Erythroid_Endothelial_correction_ids <- colnames(subset(Endothelial.combined, subset = seurat_clusters == 5)) # GYPC/CD71, lower CD34/VECAD
immune.combined$cluster_anno_l2[Erythroid_Endothelial_correction_ids] <- "Erythroid"

# Repeat for Myeloid
Myeloid.combined <- subset(immune.combined, subset = cluster_anno_l1 == 'Myeloid')
Myeloid.combined <- ScaleData(Myeloid.combined)
Myeloid.combined <- RunPCA(Myeloid.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
Myeloid.combined <- RunUMAP(Myeloid.combined, dims = 1:30)
Myeloid.combined <- FindNeighbors(Myeloid.combined, features = all_features, k.param=30)#features = all_features,
# Myeloid.combined <- FindClusters(Myeloid.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
Myeloid.combined <- FindClusters(Myeloid.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(Myeloid.combined@active.ident)))
DimPlot(Myeloid.combined, reduction = "umap", cols = cols, pt.size = 0.1, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

# Rename Myeloid clusters
MyeloidProg_ids <- colnames(subset(Myeloid.combined, subset = seurat_clusters == 0))
immune.combined$cluster_anno_l2[MyeloidProg_ids] <- "Early Myeloid Progenitor"
IntermediateMyeloid_ids <- colnames(subset(Myeloid.combined, subset = seurat_clusters == 4 | seurat_clusters == 5 | seurat_clusters == 7))
immune.combined$cluster_anno_l2[IntermediateMyeloid_ids] <- "Intermediate Myeloid"
MatureMyeloid_ids <- colnames(subset(Myeloid.combined, subset = seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters == 3 | seurat_clusters == 6 |seurat_clusters == 8 |seurat_clusters == 9))
immune.combined$cluster_anno_l2[MatureMyeloid_ids] <- "Mature Myeloid"
Undetermined_Myeloid_ids <- colnames(subset(Myeloid.combined, subset = seurat_clusters == 10))
immune.combined$cluster_anno_l2[Undetermined_Myeloid_ids] <- "Undetermined"

# Repeat for MSC
MSC.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'MSC')
MSC.combined <- ScaleData(MSC.combined)
MSC.combined <- RunPCA(MSC.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
MSC.combined <- RunUMAP(MSC.combined, dims = 1:30)
MSC.combined <- FindNeighbors(MSC.combined, features = all_features, k.param=30)#features = all_features,
# MSC.combined <- FindClusters(MSC.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
MSC.combined <- FindClusters(MSC.combined, algorithm = 2, resolution = 1) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(MSC.combined@active.ident)))
DimPlot(MSC.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()



# Repeat for Mks
Megakaryocyte.combined <- subset(immune.combined, subset = cluster_anno_l1 == 'Megakaryocyte')
Megakaryocyte.combined <- ScaleData(Megakaryocyte.combined)
Megakaryocyte.combined <- RunPCA(Megakaryocyte.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
Megakaryocyte.combined <- RunUMAP(Megakaryocyte.combined, dims = 1:30)
Megakaryocyte.combined <- FindNeighbors(Megakaryocyte.combined, features = all_features, k.param=30)#features = all_features,
# Megakaryocyte.combined <- FindClusters(Megakaryocyte.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
Megakaryocyte.combined <- FindClusters(Megakaryocyte.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(Megakaryocyte.combined@active.ident)))
DimPlot(Megakaryocyte.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()


# Rename Mk clusters
GATA1neg_Mk_ids <- colnames(subset(Megakaryocyte.combined, subset = seurat_clusters == 3))
immune.combined$cluster_anno_l2[GATA1neg_Mk_ids] <- "GATA1neg_Mks"
GATA1pos_Mk_ids <- colnames(subset(Megakaryocyte.combined, subset = seurat_clusters == 8))
immune.combined$cluster_anno_l2[GATA1pos_Mk_ids] <- "GATA1pos_Mks"
EarlyMyeloid_Mk_correction_ids <- colnames(subset(Megakaryocyte.combined, subset = seurat_clusters == 5)) # MPO expression / lack of CD11B
immune.combined$cluster_anno_l2[EarlyMyeloid_Mk_correction_ids] <- "Early Myeloid Progenitor"
MatureMyeloid_Mk_correction_ids <- colnames(subset(Megakaryocyte.combined, subset = seurat_clusters == 0 | seurat_clusters == 2)) # CD15 + CD11B, lack of MPO
immune.combined$cluster_anno_l2[MatureMyeloid_Mk_correction_ids] <- "Mature Myeloid"
Erythroid_Mk_correction_ids <- colnames(subset(Megakaryocyte.combined, subset = seurat_clusters == 4)) # GYPC, lack of CD61 or myeloid markers
immune.combined$cluster_anno_l2[Erythroid_Mk_correction_ids] <- "Erythroid"
HSPC_Mk_correction_ids <- colnames(subset(Megakaryocyte.combined, subset = seurat_clusters == 7)) # GYPC, lack of CD61 or myeloid markers
immune.combined$cluster_anno_l2[HSPC_Mk_correction_ids] <- "HSPC"
Undetermined_Mk_correction_ids <- colnames(subset(Megakaryocyte.combined, subset = seurat_clusters == 1 | seurat_clusters == 6)) # No clear CD61 expression or other markers
immune.combined$cluster_anno_l2[Undetermined_Mk_correction_ids] <- "Undetermined"

# Repeat for Stromal
Stromal.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'Stromal')
Stromal.combined <- ScaleData(Stromal.combined)
Stromal.combined <- RunPCA(Stromal.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
Stromal.combined <- RunUMAP(Stromal.combined, dims = 1:30)
Stromal.combined <- FindNeighbors(Stromal.combined, features = all_features, k.param=30)#features = all_features,
# Stromal.combined <- FindClusters(Stromal.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
Stromal.combined <- FindClusters(Stromal.combined, algorithm = 2, resolution = 1) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(Stromal.combined@active.ident)))
DimPlot(Stromal.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()


# Rename Stromal Clusters 
Endosteal_ids <- colnames(subset(Stromal.combined, subset = seurat_clusters == 3 | seurat_clusters == 11)) # cells found on the image as being endosteally located, high expression of CD56 (*note CD56 not used for clustering)
immune.combined$cluster_anno_l2[Endosteal_ids] <- "Endosteal"
MSC_Stromal_ids <- colnames(subset(Stromal.combined, subset = seurat_clusters == 1 | seurat_clusters == 12)) # FOXC1+ CXCL12+ typical MSCs
immune.combined$cluster_anno_l2[MSC_Stromal_ids] <- "MSC"
Adipocyte_Stromal_ids <- colnames(subset(Stromal.combined, subset = seurat_clusters == 0 | seurat_clusters == 2 | seurat_clusters == 4 |seurat_clusters == 5 |seurat_clusters == 7 |seurat_clusters == 8 |seurat_clusters == 9 |seurat_clusters == 10 |seurat_clusters == 13)) # Cells lining the adipose checked on image, also see CD146/VIM and lack of CD45. Cluster 8 has kind of low CD146 on cluster but checking the image they are still lining the adipose
immune.combined$cluster_anno_l2[Adipocyte_Stromal_ids] <- "Adipocyte"
Undetermined_Stromal_ids <- colnames(subset(Stromal.combined, subset = seurat_clusters == 6)) # Cells not really lining fat, expression of hematopoietic markers, hard to tell what cell type
immune.combined$cluster_anno_l2[Undetermined_Stromal_ids] <- "Undetermined"

# Repeat for Plasma Cells
Plasma_Cells.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'Plasma Cells')
Plasma_Cells.combined <- ScaleData(Plasma_Cells.combined)
Plasma_Cells.combined <- RunPCA(Plasma_Cells.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
Plasma_Cells.combined <- RunUMAP(Plasma_Cells.combined, dims = 1:30)
Plasma_Cells.combined <- FindNeighbors(Plasma_Cells.combined, features = all_features, k.param=30)#features = all_features,
# Plasma_Cells.combined <- FindClusters(Plasma_Cells.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
Plasma_Cells.combined <- FindClusters(Plasma_Cells.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(Plasma_Cells.combined@active.ident)))
DimPlot(Plasma_Cells.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

# Rename Plasma Cell Clusters
Erythroid_PC_ids <- colnames(subset(Plasma_Cells.combined, subset = seurat_clusters == 6)) # GYPC / CD71 near plasma cells in the image
immune.combined$cluster_anno_l2[Erythroid_PC_ids] <- "Erythroid"
Myeloid_PC_ids <- colnames(subset(Plasma_Cells.combined, subset = seurat_clusters == 3)) # CD15, CD11B and confirmation on image. Note CD11B normalized values not that high but due to CD11B in PCs
immune.combined$cluster_anno_l2[Myeloid_PC_ids] <- "Mature Myeloid"

# Rename THY1+ MSCs
FeatureScatter(subset(immune.combined, cluster_anno_l2 != "Autofluorescent"), feature1 = "CXCL12", feature2 = "CD90", group.by = "cluster_anno_l2") + geom_density_2d()
THY1_MSC_ids <- colnames(subset(immune.combined, subset = CD90 > 1 & cluster_anno_l2 == "MSC")) # Gate set on all cells excluding the autofluorescent cells
immune.combined$cluster_anno_l2[THY1_MSC_ids] <- "THY1+ MSC"

# Repeat for Monocytes
Monocytes.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'Monocytes')
Monocytes.combined <- ScaleData(Monocytes.combined)
Monocytes.combined <- RunPCA(Monocytes.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
Monocytes.combined <- RunUMAP(Monocytes.combined, dims = 1:30)
Monocytes.combined <- FindNeighbors(Monocytes.combined, features = all_features, k.param=30)#features = all_features,
# Monocytes.combined <- FindClusters(Monocytes.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
Monocytes.combined <- FindClusters(Monocytes.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(Monocytes.combined@active.ident)))
DimPlot(Monocytes.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()


# Rename Monocyte Clusters
EMP_Mono_ids <- colnames(subset(Monocytes.combined, subset = seurat_clusters == 1)) # absent CD14, high MPO 
immune.combined$cluster_anno_l2[EMP_Mono_ids] <- "Early Myeloid Progenitor"
NonClassical_Mono_ids <- colnames(subset(Monocytes.combined, subset = seurat_clusters == 4)) # CD14 dim, CD11c high, very specifically HLA-DR high
immune.combined$cluster_anno_l2[NonClassical_Mono_ids] <- "Non-Classical Monocyte"

# Repeat for CD8+ T-Cell
CD8.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'CD8+ T-Cell')
CD8.combined <- ScaleData(CD8.combined)
CD8.combined <- RunPCA(CD8.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
CD8.combined <- RunUMAP(CD8.combined, dims = 1:30)
CD8.combined <- FindNeighbors(CD8.combined, features = all_features, k.param=30)#features = all_features,
# CD8.combined <- FindClusters(CD8.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
CD8.combined <- FindClusters(CD8.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(CD8.combined@active.ident)))
DimPlot(CD8.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

# Rename CD8 Clusters
EMP_CD8_ids <- colnames(subset(CD8.combined, subset = seurat_clusters == 6)) # High MPO/CD33, low CD15/CD11B but higher than T-cells, low CD3e
immune.combined$cluster_anno_l2[EMP_CD8_ids] <- "Early Myeloid Progenitor"
MatureMyeloid_CD8_ids <- colnames(subset(CD8.combined, subset = seurat_clusters == 3)) # CD3e low, CD11B and CD15 hi, low MPO/CD33
immune.combined$cluster_anno_l2[MatureMyeloid_CD8_ids] <- "Mature Myeloid"
Erythroid_CD8_ids <- colnames(subset(CD8.combined, subset = seurat_clusters == 4)) # GYPC much higher than rest, low CD3E
immune.combined$cluster_anno_l2[Erythroid_CD8_ids] <- "Erythroid"
Undetermined_CD8_ids <- colnames(subset(CD8.combined, subset = seurat_clusters == 7)) # CD19 and CD3E very dim positivity
immune.combined$cluster_anno_l2[Undetermined_CD8_ids] <- "Undetermined"

# Repeat for CD4+ T-Cell
CD4.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'CD4+ T-Cell')
CD4.combined <- ScaleData(CD4.combined)
CD4.combined <- RunPCA(CD4.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
CD4.combined <- RunUMAP(CD4.combined, dims = 1:30)
CD4.combined <- FindNeighbors(CD4.combined, features = all_features, k.param=30)#features = all_features,
# CD4.combined <- FindClusters(CD4.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
CD4.combined <- FindClusters(CD4.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(CD4.combined@active.ident)))
DimPlot(CD4.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

# Rename CD4 Clusters
MatureMyeloid_CD4_ids <- colnames(subset(CD4.combined, subset = seurat_clusters == 3)) # CD3e low, CD11B and CD15 hi, MPO and CD33 are more intermediate, confirmation on image
immune.combined$cluster_anno_l2[MatureMyeloid_CD4_ids] <- "Mature Myeloid"

# Repeat for VSMC
VSMC.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'VSMC')
VSMC.combined <- ScaleData(VSMC.combined)
VSMC.combined <- RunPCA(VSMC.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
VSMC.combined <- RunUMAP(VSMC.combined, dims = 1:30)
VSMC.combined <- FindNeighbors(VSMC.combined, features = all_features, k.param=30)#features = all_features,
# VSMC.combined <- FindClusters(VSMC.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
VSMC.combined <- FindClusters(VSMC.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(VSMC.combined@active.ident)))
DimPlot(VSMC.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

# Rename VSMC Clusters
MatureMyeloid_VSMC_ids <- colnames(subset(VSMC.combined, subset = seurat_clusters == 4)) # ASMA low, CD11B and CD15 hi, MPO and CD33 are more intermediate, confirmation on image
immune.combined$cluster_anno_l2[MatureMyeloid_VSMC_ids] <- "Mature Myeloid"
EMP_VSMC_ids <- colnames(subset(VSMC.combined, subset = seurat_clusters == 5)) # ASMA low, CD11B and CD15 low, MPO and CD33 high, confirmation on image
immune.combined$cluster_anno_l2[EMP_VSMC_ids] <- "Early Myeloid Progenitor"

# Rename Schwann Cells - noticed in both Arterial and VSMC combined there were some mislabeled cells
AEC.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'AEC')
AEC.combined <- ScaleData(AEC.combined)

FeatureScatter(AEC.combined, feature1 = "CD271", feature2 = "PLP1") + geom_density_2d()
FeatureScatter(VSMC.combined, feature1 = "CD271", feature2 = "PLP1", group.by = 'seurat_clusters') + geom_density_2d() # Note - positives all in the same cluster (seruat_clusters == 1 for VSMC.combined)

Schwann_VSMC_ids <- colnames(subset(VSMC.combined, subset = PLP1 > 0.9)) # Clear population of PLP1 hi CD271 (NGFR) high cells. The actual nuclei should be Schwann Cells lining myelinated axons which penetrate the BM (axons visible on image)
immune.combined$cluster_anno_l2[Schwann_VSMC_ids] <- "Schwann Cells"
Schwann_AEC_ids <- colnames(subset(AEC.combined, subset = PLP1 > 0.9)) # Clear population of PLP1 hi CD271 (NGFR) high cells. The actual nuclei should be Schwann Cells lining myelinated axons which penetrate the BM (axons visible on image)
immune.combined$cluster_anno_l2[Schwann_AEC_ids] <- "Schwann Cells"

# Rename CD79A+ CD61+ Undetermined to B-Cells after looking at image, just B-Cells near Mks
BCell_CD79ACD61_ids <- colnames(subset(immune.combined, cluster_anno_l2 == "CD79A+ CD61+ Undetermined"))
immune.combined$cluster_anno_l2[BCell_CD79ACD61_ids] <- "B-Cells"

#saveRDS(immune.combined, "immune.combined_untilRefining_112222.RDS")
#immune.combined_111622 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.combined_updated111622.RDS")

###Refine Annotations Which Have Been Renamed----
#immune.combined <- readRDS("immune.combined_untilRefining_112222.RDS")
# Refine Early Myeloid Progenitor Annotation Because observed significant heterogeneity in terms of MPO/CD33 Positivity
EarlyMyeloidProgenitor_Refined.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'Early Myeloid Progenitor')
EarlyMyeloidProgenitor_Refined.combined <- ScaleData(EarlyMyeloidProgenitor_Refined.combined)
EarlyMyeloidProgenitor_Refined.combined <- RunPCA(EarlyMyeloidProgenitor_Refined.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
EarlyMyeloidProgenitor_Refined.combined <- RunUMAP(EarlyMyeloidProgenitor_Refined.combined, dims = 1:30)
EarlyMyeloidProgenitor_Refined.combined <- FindNeighbors(EarlyMyeloidProgenitor_Refined.combined, features = all_features, k.param=30)#features = all_features,
# EarlyMyeloidProgenitor_Refined.combined <- FindClusters(EarlyMyeloidProgenitor_Refined.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
EarlyMyeloidProgenitor_Refined.combined <- FindClusters(EarlyMyeloidProgenitor_Refined.combined, algorithm = 2, resolution = 0.5) # Louvain

VlnPlot(EarlyMyeloidProgenitor_Refined.combined, features = c("GYPC", "CD34", "MPO", "CD33", "CD141", "ASMA", "CXCL12","CD45", "CD14", "CD3e", "CD15", "CD11B", "CD79A", "PAX5", "CD123", "CD34"), pt.size = 0)

#saveRDS(EarlyMyeloidProgenitor_Refined.combined, "EarlyMyeloidProgenitor_Refined.combined")

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(EarlyMyeloidProgenitor_Refined.combined@active.ident)))
DimPlot(EarlyMyeloidProgenitor_Refined.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

# generate cluster mask for H10 EarlyMyeloidProgenitor_Refined
H10_md_EarlyMyeloidProgenitor_Refined <- subset(EarlyMyeloidProgenitor_Refined.combined@meta.data, subset = orig.ident == 'SB67_NBM27_H10_CODEX_Mesmer')
# write.csv(x= data.frame(H10_md_EarlyMyeloidProgenitor_Refined$x.coord, H10_md_EarlyMyeloidProgenitor_Refined$y.coord, H10_md_EarlyMyeloidProgenitor_Refined$CellID, H10_md_EarlyMyeloidProgenitor_Refined$seurat_clusters), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H10_EarlyMyeloidProgenitor_Refined_refined_112222_mesmer_cluster_anno_l2_annotations_111022.csv")

# Rename EarlyMyeloidProgenitor_Refined Clusters
BCell_EarlyMyeloidProgenitor_Refined_ids <- colnames(subset(EarlyMyeloidProgenitor_Refined.combined, subset = seurat_clusters == 11)) # ASMA low, CD11B and CD15 hi, MPO and CD33 are more intermediate, confirmation on image
immune.combined$cluster_anno_l2[BCell_EarlyMyeloidProgenitor_Refined_ids] <- "B-Cells"
IntermediateMyeloid_EarlyMyeloidProgenitor_Refined_ids <- colnames(subset(EarlyMyeloidProgenitor_Refined.combined, subset = seurat_clusters == 8 | seurat_clusters == 9 | seurat_clusters == 10)) # ASMA low, CD11B and CD15 low, MPO and CD33 high, confirmation on image
immune.combined$cluster_anno_l2[IntermediateMyeloid_EarlyMyeloidProgenitor_Refined_ids] <- "Intermediate Myeloid"
Artifact_EarlyMyeloidProgenitor_Refined_ids <- colnames(subset(EarlyMyeloidProgenitor_Refined.combined, subset = seurat_clusters == 0 | seurat_clusters == 2 | seurat_clusters == 12)) # ASMA low, CD11B and CD15 low, MPO and CD33 high, confirmation on image
immune.combined$cluster_anno_l2[Artifact_EarlyMyeloidProgenitor_Refined_ids] <- "Artifact"



# Refine Mature Myeloid Annotation
MatureMyeloid.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'Mature Myeloid')
MatureMyeloid.combined <- ScaleData(MatureMyeloid.combined)
MatureMyeloid.combined <- RunPCA(MatureMyeloid.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
MatureMyeloid.combined <- RunUMAP(MatureMyeloid.combined, dims = 1:30)
MatureMyeloid.combined <- FindNeighbors(MatureMyeloid.combined, features = all_features, k.param=30)#features = all_features,
# MatureMyeloid.combined <- FindClusters(MatureMyeloid.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
MatureMyeloid.combined <- FindClusters(MatureMyeloid.combined, algorithm = 2, resolution = 0.5) # Louvain


getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(MatureMyeloid.combined@active.ident)))
DimPlot(MatureMyeloid.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

VlnPlot(MatureMyeloid.combined, features = c("GYPC", "CD34", "MPO", "CD33", "CD68", "ASMA", "CXCL12","CD45", "CD14", "CD3e", "CD15", "CD11B", "CD79A", "PAX5", "CD123", "CD34"), pt.size = 0)


# generate cluster mask for H10 MatureMyeloid
H10_md_MatureMyeloid <- subset(MatureMyeloid.combined@meta.data, subset = orig.ident == 'SB67_NBM27_H10_CODEX_Mesmer')
# write.csv(x= data.frame(H10_md_MatureMyeloid$x.coord, H10_md_MatureMyeloid$y.coord, H10_md_MatureMyeloid$CellID, H10_md_MatureMyeloid$seurat_clusters), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H10_MatureMyeloid_mesmer_cluster_anno_l2_annotations_112322.csv")

# Rename MatureMyeloid Clusters
BCell_MatureMyeloid_ids <- colnames(subset(MatureMyeloid.combined, subset = seurat_clusters == 9)) # CD11B positives are all CD79A positive, confirmation on image
immune.combined$cluster_anno_l2[BCell_MatureMyeloid_ids] <- "B-Cells"
IntermediateMyeloid_MatureMyeloid_ids <- colnames(subset(MatureMyeloid.combined, subset = seurat_clusters == 11)) # CD11B low CD15 intermediate, MPO and CD33 intermediate, confirmation on image
immune.combined$cluster_anno_l2[IntermediateMyeloid_MatureMyeloid_ids] <- "Intermediate Myeloid"
Erythroid_MatureMyeloid_ids <- colnames(subset(MatureMyeloid.combined, subset = seurat_clusters == 4)) # CD11B and CD15 low, GYPC high, no other lineage markers expressed confirmation on image
immune.combined$cluster_anno_l2[Erythroid_MatureMyeloid_ids] <- "Erythroid"

# Refine Erythroid Annotation
Erythroid_refined.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'Erythroid')
Erythroid_refined.combined <- ScaleData(Erythroid_refined.combined)
Erythroid_refined.combined <- RunPCA(Erythroid_refined.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
Erythroid_refined.combined <- RunUMAP(Erythroid_refined.combined, dims = 1:30)
Erythroid_refined.combined <- FindNeighbors(Erythroid_refined.combined, features = all_features, k.param=30)#features = all_features,
# Erythroid_refined.combined <- FindClusters(Erythroid_refined.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
Erythroid_refined.combined <- FindClusters(Erythroid_refined.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(Erythroid_refined.combined@active.ident)))
DimPlot(Erythroid_refined.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

H10_md_Erythroid_refined <- subset(Erythroid_refined.combined@meta.data, subset = orig.ident == 'SB67_NBM27_H10_CODEX_Mesmer')
#write.csv(x= data.frame(H10_md_Erythroid_refined$x.coord, H10_md_Erythroid_refined$y.coord, H10_md_Erythroid_refined$CellID, H10_md_Erythroid_refined$seurat_clusters), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H10_Erythroid_refined_mesmer_cluster_anno_l2_annotations_111022.csv")

VlnPlot(Erythroid_refined.combined, features = c("GYPC", "CD34", "MPO", "CD33", "CD68", "ASMA", "CXCL12","CD45", "CD14", "CD3e", "CD15", "CD11B", "CD79A", "PAX5", "CD123", "CD34", "GATA1", "CD61", "TGFB1"), pt.size = 0)

# generate cluster mask for H10 Erythroid_refined
H10_md_Erythroid_refined <- subset(Erythroid_refined.combined@meta.data, subset = orig.ident == 'SB67_NBM27_H10_CODEX_Mesmer')
#write.csv(x= data.frame(H10_md_Erythroid_refined$x.coord, H10_md_Erythroid_refined$y.coord, H10_md_Erythroid_refined$CellID, H10_md_Erythroid_refined$seurat_clusters), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H10_Erythroid_refined_mesmer_cluster_anno_l2_annotations_112322.csv")

# Refine B-Cell Annotations 
B_Cell_refined.combined <- subset(immune.combined, subset = cluster_anno_l2 == 'B-Cells')
B_Cell_refined.combined <- ScaleData(B_Cell_refined.combined)
B_Cell_refined.combined <- RunPCA(B_Cell_refined.combined, features = all_features, approx = FALSE)
# codex <- RunUMAP(codex, method = "umap-learn", metric = "correlation", features = all_features, dims = 1:33)
B_Cell_refined.combined <- RunUMAP(B_Cell_refined.combined, dims = 1:30)
B_Cell_refined.combined <- FindNeighbors(B_Cell_refined.combined, features = all_features, k.param=30)#features = all_features,
# B_Cell_refined.combined <- FindClusters(B_Cell_refined.combined, algorithm = 4,method = 'igraph', resolution = 0.5) # Leiden (takes a longgg time)
B_Cell_refined.combined <- FindClusters(B_Cell_refined.combined, algorithm = 2, resolution = 0.5) # Louvain

getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(B_Cell_refined.combined@active.ident)))
DimPlot(B_Cell_refined.combined, reduction = "umap", cols = cols, pt.size = 0.5, repel = TRUE, label = TRUE) + NoAxes() + coord_fixed()

VlnPlot(B_Cell_refined.combined, features = c("GYPC", "CD34", "MPO", "CD33", "CD68", "ASMA", "CXCL12","CD45", "CD14", "CD3e", "CD15", "CD11B", "CD79A", "PAX5", "CD123", "CD34", "GATA1", "CD61", "TGFB1"), pt.size = 0)

H10_md_B_Cell_refined <- subset(B_Cell_refined.combined@meta.data, subset = orig.ident == 'SB67_NBM27_H10_CODEX_Mesmer')
#write.csv(x= data.frame(H10_md_B_Cell_refined$x.coord, H10_md_B_Cell_refined$y.coord, H10_md_B_Cell_refined$CellID, H10_md_B_Cell_refined$seurat_clusters), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H10_B_Cell_refined_mesmer_cluster_anno_l2_annotations_112322.csv")

# B-Cell subsets did not need to be renamed after image inspection


# Rename HSC subsets based on classical immunophenotype
HSC_Spink2hspc_refined_ids <- colnames(subset(immune.combined, subset = cluster_anno_l2 == "SPINK2+ HSPC" & CD90 > 0.5 & CD45RA < 0.7 & CD38 < 0.6)) # Lin- CD34+ CD38- CD90+ CD45RA- ,scRNA-Seq shows SPINK2+ contains HSC/GMP/CLP
immune.combined$cluster_anno_l2[HSC_Spink2hspc_refined_ids] <- "HSC"
CLP_Spink2hspc_refined_ids <- colnames(subset(immune.combined, subset = cluster_anno_l2 == "SPINK2+ HSPC" & CD45RA > 0.7 & CD38 < 0.6)) # Lin- CD34+ CD38- CD45RA+ , scRNA-Seq shows SPINK2+ contains HSC/GMP/CLP
immune.combined$cluster_anno_l2[CLP_Spink2hspc_refined_ids] <- "CLP"
GMP_Spink2hspc_refined_ids <- colnames(subset(immune.combined, subset = cluster_anno_l2 == "SPINK2+ HSPC" & CD123 > 0.5 & CD45RA > 0.7 & CD38 > 0.6)) # Lin- CD34+ CD38+ CD45RA+ CD123+, scRNA-Seq shows SPINK2+ contains HSC/GMP/CLP
immune.combined$cluster_anno_l2[GMP_Spink2hspc_refined_ids] <- "GMP"

# saveRDS(immune.combined, "immune.combined_final_112322.RDS")
# immune.combined <- readRDS("immune")
# immune.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Test_Full_Script/immune.combined_final_112322.RDS")

table(immune.combined$cluster_anno_l2)
table(immune.combined_final$cluster_anno_l2)

# save 
saveRDS(immune.combined, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.combined_final_040523.RDS")


# Filter out artifacts and undetermined clusters 
immune.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Test_Full_Script/immune.combined_final_112322.RDS")




### Make Cell Annotation Figures ----


# Rename orig.idents 
levels(as.factor(immune.combined$orig.ident))
orig.ident_annotations <- c("H10", "H14", "H32", "H37", "H36", "H38", "H26", "H35", "H33", "H41", "H39", "H27")
immune.combined <- SetIdent(immune.combined, value = "orig.ident")
names(orig.ident_annotations) <- levels(as.factor(immune.combined$orig.ident))
immune.combined <- RenameIdents(immune.combined, orig.ident_annotations)
immune.combined@meta.data["Sample_Name"] <- immune.combined@active.ident 


# Give coarse definition for filtered combined object
levels(as.factor(immune.combined$cluster_anno_l2))
coarse_annotations <- c("Adipocyte", "AEC","Artifact","Artifact","Lymphoid", "HSPC", "Lymphoid","Artifact","Lymphoid", "HSPC", "Myeloid", "Endosteal", "Erythroid", "Erythroid", "Megakaryocyte", "Megakaryocyte", 
                        "HSPC", "Myeloid", "HSPC", "HSPC", "Lymphoid", "Myeloid", "Macrophages","Myeloid", "HSPC", "Monocytes", "MSC", "Monocytes", "pDC", "Lymphoid", "Schwann Cells", "SEC", "HSPC", "MSC","Artifact","VSMC")
immune.combined <- SetIdent(immune.combined, value = "cluster_anno_l2")
names(coarse_annotations) <- levels(as.factor(immune.combined$cluster_anno_l2))
immune.combined <- RenameIdents(immune.combined, coarse_annotations)
immune.combined@meta.data["cluster_anno_coarse"] <- immune.combined@active.ident 

# Rename MSC to Adipo-MSC based on scRNA-Seq
AdipoMSC_ids <- colnames(subset(immune.combined, subset = cluster_anno_l2 == "MSC")) 
immune.combined$cluster_anno_l2[AdipoMSC_ids] <- "Adipo-MSC"

# Filter out artifacts and undetermined clusters 
#immune.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Test_Full_Script/immune.combined_final_112322.RDS")
immune.filtered <- subset(immune.combined, subset = cluster_anno_l2 != "Autofluorescent" & cluster_anno_l2 != "Undetermined" &cluster_anno_l2 != "CD44+ Undetermined" & cluster_anno_l2 != "Artifact")

#saveRDS(immune.filtered, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.filtered_FINAL.RDS")

# Figure 4C CODEX Expression Heatmap ----
avg <- AverageExpression(object = immune.filtered, group.by = 'cluster_anno_l2', slot = 'data',features = rownames(immune.filtered)) # Return average expression values across cells in each cluster for the selected genes from top3 object
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")) 
cell_order <- c("HSC", "SPINK2+ HSPC", "HSPC", "GMP", "GMP/Myeloblast", "Early Myeloid Progenitor", "Intermediate Myeloid", "Mature Myeloid", "Monocytes", "Non-Classical Monocyte", "Macrophages", "pDC", "CLP", "Immature_B_Cell", "B-Cells", "CD4+ T-Cell", "CD8+ T-Cell", "Plasma Cells", "MEP/Early Erythroblast", "CD34+ CD61+", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "Adipo-MSC", "THY1+ MSC", "Adipocyte", "Endosteal", "AEC", "SEC", "VSMC", "Schwann Cells")
cell_lineages <- c("HSPC", "HSPC", "HSPC", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Meg/E", "Meg/E", "Meg/E", "Meg/E", "Meg/E", "Meg/E", "Mesenchymal", "Mesenchymal", "Mesenchymal", "Mesenchymal", "Endothelial", "Endothelial", "Muscle", "Schwann Cells")
cell_lineage_colors <- c("#E0B0FF", "#E0B0FF", "#E0B0FF", "#A7C7E7", "#A7C7E7", "#A7C7E7", "#A7C7E7", "#A7C7E7", "#A7C7E7", "#A7C7E7", "#A7C7E7", "#A7C7E7", "#AFE1AF", "#AFE1AF", "#AFE1AF", "#AFE1AF", "#AFE1AF", "#AFE1AF", "#BDB5D5", "#BDB5D5", "#BDB5D5", "#BDB5D5", "#BDB5D5", "#BDB5D5", "#FFB6C1", "#FFB6C1", "#FFB6C1", "#FFB6C1", "#F28C28", "#F28C28", "#DD3F4E", "#FF69B4")
names(cell_lineage_colors) <- cell_lineages
cal2_cols <- c("#CF9FFF", "#E7C7DC", "#CAA7DD", "#A8A2D2", "#B6D0E2", "#2874A6", "#5599C8", "#AEC6CF", "#6495ED", "#64b8ed", "#96C5D7", "#40B5AD", "#8FCACA", "#CCE2CB", "#63BA97", "#7DB954", "#64864A", "#019477", "#953553", "#A1045A", "#a15891", "#9C58A1", "#79127F", "#BF40BF", "#FFD580", "#FFC8A2", "#FDDA0D", "#F3B0C3", "#FFBF00", "#ff9d5c", "#DD3F4E", "#FF69B4")
names(cal2_cols) <- cell_order
cell_counts <- as.numeric(table(factor(immune.filtered$cluster_anno_l2, levels = cell_order)))
anno_df <- data.frame(cell_order, cell_lineages)
ha = HeatmapAnnotation(annotation_label = c("Cell Type", "Cell Lineage"), df = anno_df,
                       border = TRUE, which = 'col', col = list(cell_order = c(cal2_cols),cell_lineages = c(cell_lineage_colors)))
avg <- as.data.frame(avg$CODEX)
Heatmap(t(scale(t(avg[,cell_order]))), name = "mat", rect_gp = gpar(col = "black", lwd = 0), col = col_fun, column_order = cell_order,
        column_title = "scRNA Average Scaled Cell Type Marker Expression", clustering_method_rows = "single", show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 8), top_annotation =  ha)

# change levels for plotting order
# immune.filtered$cluster_anno_l2 <- factor(immune.filtered$cluster_anno_l2, levels = c("HSC", "SPINK2+ HSPC", "HSPC", "GMP", "GMP/Myeloblast", "Early Myeloid Progenitor", "Intermediate Myeloid", "Mature Myeloid", "Monocytes", "Non-Classical Monocyte", "Macrophages", "pDC", "CLP", "Immature_B_Cell", "B-Cells", "CD4+ T-Cell", "CD8+ T-Cell", "Plasma Cells", "MEP/Early Erythroblast", "CD34+ CD61+", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "MSC", "THY1+ MSC", "Adipocyte", "Endosteal", "AEC", "SEC", "VSMC", "Schwann Cells", "Artifact","Autofluorescent","CD44+ Undetermined", "Undetermined"))


# Supplemental Figure S5C - coarse frequency barplot ----
immune.filtered$cluster_anno_coarse <- droplevels(immune.filtered$cluster_anno_coarse)
immune.filtered$cluster_anno_coarse <- factor(immune.filtered$cluster_anno_coarse, levels = c("HSPC", "Myeloid", "Monocytes", "Macrophage", "pDC", "Erythroid", "Megakaryocyte", "Lymphoid", "Adipocyte", "Adipo-MSC", "Endosteal", "SEC", "AEC", "VSMC", "Schwann Cells"))

cell_counts_coarse <- as.data.frame(table(immune.filtered$cluster_anno_coarse,immune.filtered$Sample_Name))
col_coarse <- c("#FDDA0D", "#FFBF00", "#AFE1AF", "#E0B0FF", "#A7C7E7", "#c3cede", "#BDB5D5", "#79127F", "#96C5D7", "#6495ED", "#FFB6C1", "#40B5AD", "#FF69B4", "#ff9d5c", "#DD3F4E")
names(col_coarse) <- c("Adipocyte", "AEC", "Lymphoid", "HSPC", "Myeloid", "Endosteal", "Erythroid", "Megakaryocyte", "Macrophages", "Monocytes", "MSC", "pDC", "Schwann Cells", "SEC", "VSMC")
ggplot(data=cell_counts_coarse, aes(x=Var2 ,y=Freq, fill=Var1)) +
  geom_bar(position="fill",stat="identity") + 
  theme_bw() + 
  labs(x='Sample') + 
  labs(y='Cell Type') + scale_fill_manual(values=col_coarse) + theme(axis.text = element_text(size=12)) + ggtitle("CODEX Cell Type Frequencies Per Sample") + RotatedAxis()



