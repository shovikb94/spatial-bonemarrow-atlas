# Transcriptomic and Spatial Proteomic Profiling Reveals the Cellular Composition and Spatial Organization of the Human Bone Marrow Microenvironment
# Author: Shovik Bandyopadhyay
# This script documents the analysis of scRNA-Seq data in this paper for Figure 1,2 and 5F

# Load necessary libraries ----
library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(RColorBrewer)
library(ggplot2)
library(VisCello)
library(viridis)
library(forcats)
library(ggrastr)
library(cowplot)
library(GSEABase)
library(SeuratDisk)
library(ComplexHeatmap)
library(DoubletFinder)
library(irlba)
library(AUCell)

# Load Doublet Finder ------
FindDoublets <- function(seurat.rna, PCs = 1:50, exp_rate = 0.02, sct = FALSE){
  # sct--do SCTransform or not
  
  ## pK identification
  sweep.res.list <- paramSweep_v3(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## Homotypic Doublet proportion Estimate
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seurat.rna <- doubletFinder_v3(seurat.rna, PCs = PCs, pN = 0.25,
                                 pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, 
                                 sct = sct)
  
  seurat.rna <- doubletFinder_v3(seurat.rna, PCs = PCs, pN = 0.25, 
                                 pK = 0.09, nExp = nExp_poi.adj,
                                 reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                                 sct = sct)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seurat.rna[['Doublet_Singlet']] = seurat.rna[[doublet_var]]
  
  mnames = names(seurat.rna@meta.data)
  seurat.rna@meta.data[, grep(mnames, pattern = '0.25_0.09')] <- NULL
  #seurat.rna = subset(seurat.rna, Doublet_Singlet == 'Singlet')
  return(seurat.rna)
}

# Remove doublets, preprocess the data, compute dimension reductions and cluster ----

# Load the cellranger output files
MACS <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/SB51_H21andH23_scRNA/H14/With_Introns/MACS/filtered_feature_bc_matrix/")
H21 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/SB51_H21andH23_scRNA/H21/filtered_feature_bc_matrix/")
H23 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/SB51_H21andH23_scRNA/H23/filtered_feature_bc_matrix/")
H24 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/SB62_scRNASeq_S2/H24/filtered_feature_bc_matrix/")
H32 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/SB62_scRNASeq_S2/H32/filtered_feature_bc_matrix/")
H33 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/SB62_scRNASeq_S2/H33/filtered_feature_bc_matrix/")
H34 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/SB62_scRNASeq_S2/H34/filtered_feature_bc_matrix/")
H35 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/SB62_scRNASeq_S2/H35/filtered_feature_bc_matrix/")
H36 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/H36/filtered_feature_bc_matrix/")
H38 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/H38/filtered_feature_bc_matrix/")
H39 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/H39/filtered_feature_bc_matrix/")
H41 <- Read10X(data.dir = "~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/H41/filtered_feature_bc_matrix/")


# Initialize the Seurat object with the raw (non-normalized data).
MACS <- CreateSeuratObject(counts = MACS, project = "H14_MACS", min.cells = 3, min.features = 100)
MACS
H21 <- CreateSeuratObject(counts = H21, project = "H21", min.cells = 3, min.features = 100)
H21
H23 <- CreateSeuratObject(counts = H23, project = "H23", min.cells = 3, min.features = 100)
H23
H24 <- CreateSeuratObject(counts = H24, project = "H24", min.cells = 3, min.features = 100)
H24
H32 <- CreateSeuratObject(counts = H32, project = "H32", min.cells = 3, min.features = 100)
H32
H33 <- CreateSeuratObject(counts = H33, project = "H33", min.cells = 3, min.features = 100)
H33
H34 <- CreateSeuratObject(counts = H34, project = "H34", min.cells = 3, min.features = 100)
H34
H35 <- CreateSeuratObject(counts = H35, project = "H35", min.cells = 3, min.features = 100)
H35
H36 <- CreateSeuratObject(counts = H36, project = "H36", min.cells = 3, min.features = 100)
H36
H38 <- CreateSeuratObject(counts = H38, project = "H38", min.cells = 3, min.features = 100)
H38
H39 <- CreateSeuratObject(counts = H39, project = "H39", min.cells = 3, min.features = 100)
H39
H41 <- CreateSeuratObject(counts = H41, project = "H41", min.cells = 3, min.features = 100)
H41

# Preprocess samples individually to remove doublets -----
ob.list <- list(MACS, H21, H23,H24,H32,H33,H34,H35, H36, H38, H39, H41)
for (i in 1:length(ob.list)) {
  ob.list[[i]][["percent.mt"]] <- PercentageFeatureSet(ob.list[[i]], pattern = "^MT-")
  ob.list[[i]] <- subset(ob.list[[i]], subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 10)
  # ob.list[[i]] <- subset(ob.list[[i]], cells = cells.use)
  ob.list[[i]] <- NormalizeData(ob.list[[i]])
  ob.list[[i]] <- FindVariableFeatures(ob.list[[i]])
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], features = VariableFeatures(object = ob.list[[i]]))
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:30, reduction.name = "UMAP_dim30", reduction.key = "UMAP_dim30_")
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:50, reduction.name = "UMAP_dim50", reduction.key = "UMAP_dim50_")
  ob.list[[i]] <- FindNeighbors(ob.list[[i]], dims = 1:30)
  ob.list[[i]] <- FindClusters(ob.list[[i]], resolution = 1)
  ob.list[[i]] <- FindDoublets(ob.list[[i]], PCs = 1:30, sct = FALSE, exp_rate = (length(colnames(ob.list[[i]]))/125000))
}

MACS <- ob.list[[1]]
H21 <- ob.list[[2]]
H23 <- ob.list[[3]]
H24 <- ob.list[[4]]
H32 <- ob.list[[5]]
H33 <- ob.list[[6]]
H34 <- ob.list[[7]]
H35 <- ob.list[[8]]
H36 <- ob.list[[9]]
H38 <- ob.list[[10]]
H39 <- ob.list[[11]]
H41 <- ob.list[[12]]

# Merge all seurat objects together 
combined <- merge(x = MACS, y = c(H21, H23,H24, H32, H33,H34,H35,H36, H38, H39, H41), add.cell.ids = c("H14_MACS", "H21", "H23","H24", "H32", "H33", "H34", "H35", "H36", "H38", "H39","H41"), project = "SB62_NBM")

ob.list <- list(combined)

for (i in 1:length(ob.list)) {
  ob.list[[i]][["percent.mt"]] <- PercentageFeatureSet(ob.list[[i]], pattern = "^MT-")
  ob.list[[i]] <- subset(ob.list[[i]], subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 10)
  ob.list[[i]] <- subset(ob.list[[i]], subset = Doublet_Singlet == "Singlet")
  ob.list[[i]] <- NormalizeData(ob.list[[i]])
  ob.list[[i]] <- FindVariableFeatures(ob.list[[i]])
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], features = VariableFeatures(object = ob.list[[i]]))
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:30, reduction.name = "UMAP_dim30", reduction.key = "UMAP_dim30_")
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:50, reduction.name = "UMAP_dim50", reduction.key = "UMAP_dim50_")
  ob.list[[i]] <- FindNeighbors(ob.list[[i]], dims = 1:30)
  ob.list[[i]] <- FindClusters(ob.list[[i]], algorithm = 2, resolution = 1)
}

combined <- FindClusters(combined, algorithm = 2, resolution = 1.5)

combined <- ob.list[[1]] 

# Find and save all markers for combined UMAP
combined_markers <- FindAllMarkers(combined, max.cells.per.ident = 1000) # Downsample 1000 cells per cluster
write.csv(combined_markers, "SB66_combined_markers.csv")

# remove contaminating cells and doublets

# Cluster 46 - muscle
# Cluster 41 - CXCL14+ Fibroblasts
# Cluster 45 - Doublets Plasma + T
combined <- readRDS("SB66_combined.RDS")
combined <- subset(combined, subset = seurat_clusters != 46 & seurat_clusters != 41 & seurat_clusters != 45)

ob.list <- list(combined)

for (i in 1:length(ob.list)) {
  ob.list[[i]][["percent.mt"]] <- PercentageFeatureSet(ob.list[[i]], pattern = "^MT-")
  ob.list[[i]] <- subset(ob.list[[i]], nFeature_RNA > 300 & nCount_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
  ob.list[[i]] <- subset(ob.list[[i]], subset = Doublet_Singlet == "Singlet")
  ob.list[[i]] <- NormalizeData(ob.list[[i]])
  ob.list[[i]] <- FindVariableFeatures(ob.list[[i]])
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], features = VariableFeatures(object = ob.list[[i]]))
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:30, reduction.name = "UMAP_dim30", reduction.key = "UMAP_dim30_")
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:50, reduction.name = "UMAP_dim50", reduction.key = "UMAP_dim50_")
  ob.list[[i]] <- FindNeighbors(ob.list[[i]], dims = 1:30)
  ob.list[[i]] <- FindClusters(ob.list[[i]], algorithm = 2, resolution = 1.5)
}


combined <- ob.list[[1]] 

# Find and save all markers for combined UMAP
combined_markers <- FindAllMarkers(combined, max.cells.per.ident = 1000) # Downsample 1000 cells per cluster
write.csv(combined_markers, "SB66_combined_markers.csv")

# Annotate the cells in the final dataset ----

# Prepare seurat object for viscello
library(VisCello)
DefaultAssay(combined) <- "RNA"
fmeta <- data.frame(symbol = rownames(combined)) 
rownames(fmeta) <- fmeta$symbol
eset <- new("ExpressionSet",
            assayData = assayDataNew("environment", exprs=combined@assays$RNA@counts, 
                                     norm_exprs = combined@assays$RNA@data),
            phenoData =  new("AnnotatedDataFrame", data = combined@meta.data),
            featureData = new("AnnotatedDataFrame", data = fmeta))
saveRDS(eset, "VisCello/eset.rds") 
# Creating a cello for all the cells
cello <- new("Cello", name = "Combined All Cells SB66", idx = 1:ncol(eset)) # Index is basically the column index, here all cells are included 
# Code for computing dimension reduction, not all of them is necessary, and you can input your own dimension reduction result into the cello@proj list.
# It is also recommended that you first filter your matrix to remove low expression genes and cells, and input a matrix with variably expressed genes
cello@proj <- list('PCA' = combined@reductions$pca@cell.embeddings,
                   'UMAP_dim30' =combined@reductions$UMAP_dim30@cell.embeddings,
                   'UMAP_dim50' =combined@reductions$UMAP_dim50@cell.embeddings)
clist <- list()
clist[["SB66 Global dataset"]] <- cello
saveRDS(clist, "VisCello/clist.rds") 

# Rename clusters after manual inspection of viscello obj and markers csv
combined <- SetIdent(combined, value = "seurat_clusters")
new.cluster.ids.combined <- c("Adipo-MSC", "Plasma Cell", "Plasma Cell", "Plasma Cell", "THY1+ MSC", "Osteo-MSC",
                              "Neutrophil", "SEC", "Pro-B", "CD4+ T-Cell", "Osteoblast", "VSMC", "HSC",
                              "Late Myeloid", "Late Erythroid", "GMP", "MEP", "OsteoFibro-MSC",
                              "CD8+ T-Cell", "MPP", "Early Myeloid Progenitor", "GMP", "Fibro-MSC", "RNAlo MSC", "Monocyte", "Mature B",
                              "AEC", "Cycling HSPC", "Pre-Pro B", "Late Myeloid", "Erythroblast", "Plasma Cell", "VSMC",
                              "Pre-B", "pDC", "Ba/Eo/Ma", "AEC", "CLP", "Plasma Cell",
                              "Cycling DCs", "RBC", "Plasma Cell", "MPP", "Macrophages", "RBC", "NKT Cell")

names(new.cluster.ids.combined) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids.combined)
DimPlot(combined, label = TRUE, reduction = "UMAP_dim30") + coord_fixed() + NoAxes() + NoLegend()
combined$cluster_anno_l2 <- combined@active.ident


# Mks do not cluster separately, need to subcluster erythroblasts to properly annotate



# subset erythroblast cluster and rename Mks 
EBs <- subset(combined, cluster_anno_l2 == "Erythroblast")
EBs <- NormalizeData(EBs)
EBs <- FindVariableFeatures(EBs)
EBs <- ScaleData(EBs)
EBs <- RunPCA(EBs, reduction.name = "EB_pca", features = VariableFeatures(EBs))
EBs <- RunUMAP(EBs, dims = 1:10, reduction = "EB_pca", reduction.name = "EB_UMAP_dim10")
EBs <- RunUMAP(EBs, dims = 1:20, reduction = "EB_pca", reduction.name = "EB_UMAP_dim20")
EBs <- RunUMAP(EBs, dims = 1:30, reduction = "EB_pca", reduction.name = "EB_UMAP_dim30")
EBs <- RunUMAP(EBs, dims = 1:40, reduction = "EB_pca", reduction.name = "EB_UMAP_dim40")
EBs <- RunUMAP(EBs, dims = 1:50, reduction = "EB_pca", reduction.name = "EB_UMAP_dim50")
EBs <- FindNeighbors(EBs, dims = 1:20)
EBs <- FindClusters(EBs, algorithm = 2, resolution = 0.5, group.singletons = FALSE)
DimPlot(EBs, label = TRUE, reduction = "EB_UMAP_dim50") + coord_fixed() + NoAxes() + NoLegend()
FeaturePlot(EBs, features = c("ITGA2B", "GATA2"), reduction = "EB_UMAP_dim50", coord.fixed = TRUE) & NoAxes() 
VlnPlot(EBs, features = c("ITGA2B","ITGB3", "GATA2" ,"GYPC")) # can see that cluster 4 is most megakaryocytic
# Rename Mk cluster based on CD41 expression
Mk_ids <- colnames(subset(EBs, subset = seurat_clusters == 4)) # GATA2/CD41/CD61 expression
combined$cluster_anno_l2 <- as.character(combined$cluster_anno_l2)
combined$cluster_anno_l2[Mk_ids] <- "Megakaryocyte" 
combined$cluster_anno_l2 <- as.factor(combined$cluster_anno_l2)
levels(combined$cluster_anno_l2)
DimPlot(combined, label= TRUE, repel = TRUE, group.by = 'cluster_anno_l2') + NoAxes() + coord_fixed()

# remove NKT cluster because they are likely doublets (T-cell markers, Osteolineage markers co-expressed)
combined <- subset(combined, cluster_anno_l2 != "NKT Cell")

# Make coarse annotations
cluster_anno_coarse <- c("Mesenchymal", "Endothelial", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Mesenchymal", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Mesenchymal", "Mesenchymal", "Mesenchymal", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Hematopoietic", "Mesenchymal", "Endothelial", "Mesenchymal", "Muscle")
combined <- SetIdent(combined, value = "cluster_anno_l2")
names(cluster_anno_coarse) <- levels(combined)
combined <- RenameIdents(combined, cluster_anno_coarse)
DimPlot(combined, label = TRUE, reduction = "UMAP_dim30") + coord_fixed() + NoAxes() + NoLegend()
combined$cluster_anno_coarse <- combined@active.ident
coarse_cols <- c('#f88379','#FFA750','#6495ED','#B8242D')
DimPlot(combined, label = TRUE, reduction = "UMAP_dim30", cols=coarse_cols) + coord_fixed() + NoAxes() + NoLegend()

# save cluster_anno_l2 colors
cal2_cols <- c("#FFD580", "#FFBF00", "#B0DFE6", "#7DB954", "#64864A", "#8FCACA", "#4682B4", "#CAA7DD", "#B6D0E2", "#a15891", "#FED7C3", "#A8A2D2", "#CF9FFF", "#9C58A1", "#2874A6", "#96C5D7", "#63BA97", "#BF40BF", "#953553", "#6495ED", "#E7C7DC", "#5599C8", "#FA8072", "#F3B0C3", "#F89880", "#40B5AD", "#019477", "#97C1A9", "#C6DBDA", "#CCE2CB", "#79127F", "#FFC5BF", "#ff9d5c", "#FFC8A2", "#DD3F4E")
cal2_col_names <- c("Adipo-MSC", "AEC", "Ba/Eo/Ma", "CD4+ T-Cell", "CD8+ T-Cell", "CLP", "Cycling DCs", "Cycling HSPC", "Early Myeloid Progenitor", "Erythroblast", "Fibro-MSC", "GMP", "HSC", "Late Erythroid", "Late Myeloid", "Macrophages", "Mature B", "Megakaryocyte", "MEP", "Monocyte", "MPP", "Neutrophil", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "pDC", "Plasma Cell", "Pre-B", "Pre-Pro B", "Pro-B", "RBC", "RNAlo MSC", "SEC", "THY1+ MSC", "VSMC")
names(cal2_cols) <- cal2_col_names

combined <- SetIdent(combined, value = "cluster_anno_l2")
combined_markers_anno <- FindAllMarkers(combined, max.cells.per.ident = 500) # Downsample 1000 cells per cluster
write.csv(combined_markers_anno, "SB66_combined_markers_annotated.csv")


# Figure 1 Generation ----

# in text QC metrics
combined@meta.data %>% dplyr::summarise(c(mean(nCount_RNA), mean(nFeature_RNA), mean(percent.mt)))


# Supplemental Figure 1C QC Metrics -----
VlnPlot(combined, group.by ="cluster_anno_coarse", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0) & scale_fill_manual(values =c("#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896"))
VlnPlot(combined, group.by ="cluster_anno_coarse", features = c("nCount_RNA"), pt.size=0, y.max = quantile(combined$nCount_RNA,0.99)) & scale_fill_manual(values =c("#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896"))

combined@meta.data %>% dplyr::group_by(cluster_anno_coarse)  %>% summarise(median(nFeature_RNA), median(nCount_RNA), median(percent.mt)) %>% gt() -> supp_fig_1D
gtsave(supp_fig_1D, filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure1/Supp_Fig_1D.pdf")

# Supplemental Figure 1D Inflammation ----- 
#Run AUCell
combined.small <- subset(combined, downsample = 300) # note the results may be subtly different between runs because of downsampling step
fmeta <- data.frame(symbol = rownames(combined.small)) 
rownames(fmeta) <- fmeta$symbol
eset_combined.small <- new("ExpressionSet",
                           assayData = assayDataNew("environment", exprs=combined.small@assays$RNA@counts, 
                                                    norm_exprs = combined.small@assays$RNA@data),
                           phenoData =  new("AnnotatedDataFrame", data = combined.small@meta.data),
                           featureData = new("AnnotatedDataFrame", data = fmeta))
exprMatrix <- exprs(eset_combined.small)
geneSets <- getGmt("~/Downloads/h.all.v2022.1.Hs.symbols.gmt")
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets)))
# add random noise
# Random
set.seed(321)
extraGeneSets <- c(
  GeneSet(sample(rownames(exprMatrix), 50), setName="Random (50g)"),
  GeneSet(sample(rownames(exprMatrix), 500), setName="Random (500g)"))

countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
# Housekeeping-like
extraGeneSets <- c(extraGeneSets,
                   GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK-like (100g)"))

geneSets <- GeneSetCollection(c(geneSets,extraGeneSets))
names(geneSets)
# Start to build the cell rankings
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=FALSE)
#save(cells_rankings, file="AUCell/cells_rankings.RData")

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
#save(cells_AUC, file="AUCell/AUCecells_AUC.RData")

set.seed(123)
par(mfrow=c(3,3)) 
#cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=FALSE, assign=TRUE)
#save(cells_assignment, file="AUCell/cellassignment_AUC.RData")

auc_test <- as.data.frame(t(cells_AUC@assays@data$AUC))
combined_AUC <- AddMetaData(combined.small, auc_test)


FeaturePlot(object = combined_AUC, features = c("HALLMARK_INFLAMMATORY_RESPONSE"), coord.fixed = TRUE) & NoAxes()



# Supplemental Figure 1E UMAP Feature Plots -----
p1 <- FeaturePlot(object = combined, raster = TRUE, raster.dpi = c(1028,1028), features = c("CXCL12","NCAM1","CDH5", "PTPRC", "MZB1", "CSF3R"), cols = brewer.pal(n = 100, name = "Reds"),ncol = 6, max.cutoff = 'q99', coord.fixed = TRUE) & NoAxes() 
p1_raster <- rasterize(p1)
ggsave(p1_raster,file = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure1/PanelE_scRNA_MarkerFeaturePlots.pdf", device = "pdf", width = 15, height = 2)

# Supplemental Figure 1F Azimuth Comparison ----
azimuth <- readRDS("~/Documents/AML_Microenvironment/SB36_CibersortX/Azimuth_BM_Reference/ref.Rds")
DimPlot(azimuth, reduction = "refUMAP", group.by = "celltype.l2") + coord_fixed() + NoAxes()
table(azimuth$celltype.l2)

# Figure 1B - UMAP Showing all of the Cell Types Captured in Our Analysis ----
p1 <- DimPlot(combined,group.by = 'cluster_anno_l2',raster = TRUE, raster.dpi = c(1028,1028),label = TRUE, repel = TRUE, reduction = "UMAP_dim30", cols=cal2_cols) + coord_fixed() + NoAxes() + NoLegend()
p1_raster <- rasterize(p1)
ggsave(p1_raster,file = "~/Documents/Manuscripts/NBM_Atlas/Figures/Figure1/PanelB_scRNA_UMAP.pdf")

# Change levels 
cell_order <- c("HSC", "MPP","Cycling HSPC", "MEP", "Erythroblast", "Late Erythroid", "RBC","Megakaryocyte","GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil", "Monocyte","Macrophages","Ba/Eo/Ma", "Cycling DCs", "pDC", "CLP", "Pre-Pro B", "Pro-B", "Pre-B", "Mature B", "Plasma Cell","CD4+ T-Cell", "CD8+ T-Cell", "NKT Cell", "AEC", "SEC", "VSMC", "THY1+ MSC", "Adipo-MSC", "OsteoFibro-MSC", "Osteo-MSC", "Osteoblast", "Fibro-MSC", "RNAlo MSC")
combined$cluster_anno_l2 <- factor(combined$cluster_anno_l2, levels = cell_order)

# Figure 1C Create barplot showing lineage frequencies ----

## Create the level 1 (coarse) annotation based on broad lineages
# cluster_anno_l1 <- c("Mesenchymal", "Endothelial", "Myeloid", "Lymphoid", "Lymphoid", "HSPC","Myeloid","HSPC", "Myeloid", "Meg/E", "Mesenchymal","HSPC","HSPC","Meg/E","Myeloid","Myeloid","Lymphoid","Meg/E","HSPC","Myeloid","HSPC","Myeloid","Mesenchymal","Mesenchymal","Mesenchymal","Myeloid","Lymphoid","Lymphoid","Lymphoid","Lymphoid","Meg/E", "Mesenchymal", "Endothelial",  "Mesenchymal", "Muscle") # First make level 1 annotation resolved for the major hematopoietic lineages
cluster_anno_l1 <- c("HSPC", "HSPC", "HSPC", "Meg/E", "Meg/E", "Meg/E", "Meg/E", "Meg/E", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Endothelial", "Endothelial", "Muscle", "Mesenchymal", "Mesenchymal", "Mesenchymal", "Mesenchymal", "Mesenchymal", "Mesenchymal", "Mesenchymal")
combined <- SetIdent(combined, value = "cluster_anno_l2")
names(cluster_anno_l1) <- levels(combined)
combined <- RenameIdents(combined, cluster_anno_l1)
combined$cluster_anno_l1 <- combined@active.ident
combined$cluster_anno_l1 <- factor(combined$cluster_anno_l1, levels = c("HSPC", "Myeloid", "Lymphoid", "Meg/E", "Mesenchymal", "Endothelial", "Muscle"))
cal1_cols <- c("#E0B0FF", "#A7C7E7", "#AFE1AF", "#BDB5D5", "#FFB6C1", "#F28C28", "#DD3F4E")
DimPlot(combined, group.by = 'cluster_anno_l1',cols = cal1_cols, label = TRUE, reduction = "UMAP_dim30") + coord_fixed() + NoAxes() + NoLegend() # Just make sure reannotation has been done properly

combined$cluster_anno_l2 <- droplevels(combined$cluster_anno_l2)
col_coarse <- c("#FFB6C1", "#FFBF00", "#A7C7E7", "#AFE1AF", "#AFE1AF", "#AFE1AF", "#40B5AD", "#E0B0FF", "#A7C7E7", "#BDB5D5", "#FFB6C1", "#A7C7E7", "#E0B0FF", "#BDB5D5", "#A7C7E7", "#96C5D7", "#AFE1AF", "#CAA7DD", "#BDB5D5", "#6495ED", "#E0B0FF", "#A7C7E7", "#c3cede", "#c3cede", "#FFB6C1", "#40B5AD", "#AFE1AF", "#AFE1AF", "#AFE1AF", "#AFE1AF", "#BDB5D5", "#FFB6C1", "#ff9d5c", "#FFB6C1", "#DD3F4E")
names_coarse <- c("Mesenchymal", "AEC", "Myeloid", "Lymphoid", "Lymphoid", "Lymphoid", "pDC", "HSPC", "Myeloid", "Erythroid", "Mesenchymal", "Myeloid", "HSPC", "Erythroid", "Myeloid", "Macrophage", "Lymphoid", "Megakaryocyte", "Erythroid", "Monocyte", "HSPC", "Myeloid", "Osteolineage", "Osteolineage", "Mesenchymal", "pDC", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Erythroid", "Mesenchymal", "SEC", "Mesenchymal", "VSMC")
names(col_coarse) <- names_coarse

# Figure 1C - Bar chart showing the counts of each lineage assayed -----
## Figure 1C (left)
p2 <- combined@meta.data %>% ggplot(aes(y=forcats::fct_rev(forcats::fct_infreq(cluster_anno_l1)), fill = cluster_anno_l1)) + geom_bar(stat = 'count') +
  theme_bw() + 
  scale_fill_manual(values = cal1_cols, labels = c("HSPC", "Myeloid", "Lymphoid", "Meg/E", "Mesenchymal", "Endothelial", "Muscle")) + 
  theme(axis.text = element_text(size=12)) 

cell_counts <- as.data.frame(table(combined$cluster_anno_l1,combined$orig.ident))

## Figure 1D (right)
p3 <- ggplot(data=cell_counts, aes(x=forcats::fct_rev(Var2),y=Freq, fill=Var1)) +
  geom_bar(position="fill",stat="identity") + 
  coord_flip() +
  theme_bw() + 
  labs(y='Cell Type Frequency') + 
  labs(x='') + scale_fill_manual(values=cal1_cols) + theme(axis.text = element_text(size=16)) + ggtitle("CODEX Cell Type Frequencies Per Sample") + scale_x_discrete(labels = c("H41", "H39", "H38", "H36", "H35", "H34", "H33", "H32", "H24", "H23", "H21", "H14"))

# Figure 1D Heatmap ----
library(ComplexHeatmap)
library(circlize)
combined_markers_anno %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = 1/p_val_adj) -> top10 # Get top 10 most significant genes for each cell type
combined <- SetIdent(combined, value = "cluster_anno_l2")
avg <- AverageExpression(object = combined, group.by = 'cluster_anno_l2', slot = 'data',features = top10$gene) # Return average expression values across cells in each cluster for the selected genes from top3 object

col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")) 
cell_order <- c("HSC", "MPP","Cycling HSPC", "MEP", "Erythroblast", "Late Erythroid", "RBC","Megakaryocyte","GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil", "Monocyte","Macrophages","Ba/Eo/Ma", "Cycling DCs", "pDC", "CLP", "Pre-Pro B", "Pro-B", "Pre-B", "Mature B", "Plasma Cell","CD4+ T-Cell", "CD8+ T-Cell", "AEC", "SEC", "VSMC", "THY1+ MSC", "Adipo-MSC", "OsteoFibro-MSC", "Osteo-MSC", "Osteoblast", "Fibro-MSC", "RNAlo MSC")
cell_lineages <- c("HSPC", "HSPC", "HSPC", "Meg/E", "Meg/E", "Meg/E", "Meg/E","Meg/E","Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Endothelial", "Endothelial", "Muscle", "Mesenchymal", "Mesenchymal", "Mesenchymal", "Mesenchymal", "Mesenchymal", "Mesenchymal","Mesenchymal")
cell_lineage_colors <- c("HSPC" = "steelblue", "HSPC"= "steelblue", "HSPC"= "steelblue", "Meg/E" = 'darkorchid1', "Meg/E" = 'darkorchid1', "Meg/E" = 'darkorchid1', "Meg/E" = 'darkorchid1', "Myeloid" = 'salmon', "Myeloid" = 'salmon', "Myeloid" = 'salmon', "Myeloid" = 'salmon', "Myeloid" = 'salmon', "Myeloid" = 'salmon', "Myeloid" = 'salmon', "Myeloid" = 'salmon', "Myeloid" = 'salmon', "Lymphoid" = 'darkgoldenrod1', "Lymphoid"= 'darkgoldenrod1', "Lymphoid"= 'darkgoldenrod1', "Lymphoid" = 'darkgoldenrod1', "Lymphoid" = 'darkgoldenrod1', "Lymphoid" = 'darkgoldenrod1', "Lymphoid"= 'darkgoldenrod1', "Lymphoid"= 'darkgoldenrod1', "Endothelial" = 'lavender', "Endothelial" = 'lavender', "Muscle" = 'cornflowerblue', "Mesenchymal" = 'darkolivegreen2', "Mesenchymal" = 'darkolivegreen2', "Mesenchymal" = 'darkolivegreen2', "Mesenchymal" = 'darkolivegreen2', "Mesenchymal" = 'darkolivegreen2', "Mesenchymal" = 'darkolivegreen2',"Mesenchymal" = 'darkolivegreen2')
cell_lineage_colors <- c("HSPC" = "#E0B0FF", "HSPC"= "#E0B0FF", "HSPC"= "#E0B0FF", "Meg/E" = "#BDB5D5", "Meg/E" = "#BDB5D5", "Meg/E" = "#BDB5D5", "Meg/E" = "#BDB5D5", "Meg/E" = "#BDB5D5", "Myeloid" = "#A7C7E7", "Myeloid" = "#A7C7E7", "Myeloid" = "#A7C7E7", "Myeloid" = "#A7C7E7", "Myeloid" = "#A7C7E7", "Myeloid" = "#A7C7E7", "Myeloid" = "#A7C7E7", "Myeloid" = "#A7C7E7", "Myeloid" = "#A7C7E7", "Lymphoid" = "#AFE1AF", "Lymphoid"= "#AFE1AF", "Lymphoid"= "#AFE1AF", "Lymphoid" = "#AFE1AF", "Lymphoid" = "#AFE1AF", "Lymphoid" = "#AFE1AF", "Lymphoid"= "#AFE1AF", "Lymphoid"= "#AFE1AF", "Endothelial" = "#F28C28", "Endothelial" = "#F28C28", "Muscle" = "#DD3F4E", "Mesenchymal" = "#FFB6C1", "Mesenchymal" = "#FFB6C1", "Mesenchymal" = "#FFB6C1", "Mesenchymal" = "#FFB6C1", "Mesenchymal" = "#FFB6C1", "Mesenchymal" = "#FFB6C1", "Mesenchymal" ="#FFB6C1")
cell_counts <- as.numeric(table(factor(combined$cluster_anno_l2, levels = cell_order)))
names(cell_order) <- cal2_cols
anno_df <- data.frame(cell_order, cell_lineages)
names(cal2_cols) <- levels(droplevels(as.factor(combined$cluster_anno_l2)))
ha = HeatmapAnnotation(annotation_label = c("Cell Type", "Cell Lineage"), df = anno_df,
                       border = TRUE, which = 'col', col = list(cell_order = cal2_cols,cell_lineages = c(cell_lineage_colors)))
avg <- as.data.frame(avg$RNA)
genes_show = c("AVP","SPINK2","IL7R", "GZMA","DNTT","CD79B","PAX5","MZB1", "MPO", "SRGN","S100A9" ,"LYZ", "CD14", "C1Q","GATA1", "ITGA2B", "HBA1","HBA2", "CXCL12","LEPR", "PDGFRA", "NCAM1","FLT1","EMCN", "TAGLN", "ACTA2")
Heatmap(t(scale(t(avg[,cell_order]))), name = "mat", rect_gp = gpar(col = "black", lwd = 0), col = col_fun, column_order = cell_order,
        column_title = "scRNA Average Scaled Cell Type Marker Expression", clustering_method_rows = "single", show_row_names = FALSE, row_names_gp = grid::gpar(fontsize = 8), top_annotation =  ha) +
  rowAnnotation(link=anno_mark(at=which(rownames(avg) %in% genes_show), labels = rownames(avg)[rownames(avg) %in% genes_show], which = "rows"))


# Perform Analysis of MSC and Endothelial Cell Subsets -----

# Subset MSCs, Re-Normalize/Scale, Compute PCA and UMAP, leaving out the RNAlo MSCs ----
MSCs <- subset(combined, cluster_anno_l2 == "Fibro-MSC" | cluster_anno_l2 == "THY1+ MSC" | cluster_anno_l2 == "Adipo-MSC" | cluster_anno_l2 == "OsteoFibro-MSC" | cluster_anno_l2 == "Osteo-MSC" | cluster_anno_l2 == "Osteoblast")
MSCs <- NormalizeData(MSCs)
MSCs <- FindVariableFeatures(MSCs)
MSCs <- ScaleData(MSCs)
MSCs <- RunPCA(MSCs, reduction.name = "MSC_pca", features = VariableFeatures(MSCs))
MSCs <- RunUMAP(MSCs, dims = 1:10, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim10")
MSCs <- RunUMAP(MSCs, dims = 1:20, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim20")
MSCs <- RunUMAP(MSCs, dims = 1:30, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim30")
MSCs <- RunUMAP(MSCs, dims = 1:40, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim40")
MSCs <- RunUMAP(MSCs, dims = 1:50, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim50")

DimPlot(MSCs, label = TRUE, reduction = "MSC_UMAP_dim50") + coord_fixed() + NoAxes() + NoLegend()
VlnPlot(subset(combined, cluster_anno_l2 == "Fibro-MSC" | cluster_anno_l2 == "THY1+ MSC" | cluster_anno_l2 == "Adipo-MSC" | cluster_anno_l2 == "OsteoFibro-MSC" | cluster_anno_l2 == "Osteo-MSC" | cluster_anno_l2 == "Osteoblast"
               | cluster_anno_l2 == "RNAlo MSC"), features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), pt.size = 0, group.by = 'cluster_anno_l2')

MSCs$original_seurat_clusters <- MSCs@active.ident
MSCs <- FindNeighbors(MSCs, dims = 1:20)
MSCs <- FindClusters(MSCs, algorithm = 2, resolution = 0.25, group.singletons = FALSE)

# update viscello with annotated combined obj
DefaultAssay(combined) <- "RNA"
fmeta <- data.frame(symbol = rownames(combined)) 
rownames(fmeta) <- fmeta$symbol
eset <- new("ExpressionSet",
            assayData = assayDataNew("environment", exprs=combined@assays$RNA@counts, 
                                     norm_exprs = combined@assays$RNA@data),
            phenoData =  new("AnnotatedDataFrame", data = combined@meta.data),
            featureData = new("AnnotatedDataFrame", data = fmeta))
saveRDS(eset, "VisCello/eset.rds") 
# Creating a cello for all the cells
cello <- new("Cello", name = "Combined All Cells SB66", idx = 1:ncol(eset)) # Index is basically the column index, here all cells are included 
# Code for computing dimension reduction, not all of them is necessary, and you can input your own dimension reduction result into the cello@proj list.
# It is also recommended that you first filter your matrix to remove low expression genes and cells, and input a matrix with variably expressed genes
cello@proj <- list('PCA' = combined@reductions$pca@cell.embeddings,
                   'UMAP_dim30' =combined@reductions$UMAP_dim30@cell.embeddings,
                   'UMAP_dim50' =combined@reductions$UMAP_dim50@cell.embeddings)
clist <- list()
clist[["SB66 Global dataset"]] <- cello
saveRDS(clist, "VisCello/clist.rds") 


# Add zoom of MSC to Viscello 
zoom_type <- "MSC_Trajectory"
cur_idx <-  which(colnames(eset) %in% colnames(MSCs))
cello <- new("Cello", name = zoom_type, idx = cur_idx) 
cello@proj <- list('MSC_PCA' = MSCs@reductions$MSC_pca@cell.embeddings,
                   'MSC_UMAP_dim10'=MSCs@reductions$MSC_UMAP_dim10@cell.embeddings, 
                   'MSC_UMAP_dim20'=MSCs@reductions$MSC_UMAP_dim20@cell.embeddings, 
                   'MSC_UMAP_dim30'=MSCs@reductions$MSC_UMAP_dim30@cell.embeddings, 
                   'MSC_UMAP_dim40'=MSCs@reductions$MSC_UMAP_dim40@cell.embeddings,
                   'MSC_UMAP_dim50'=MSCs@reductions$MSC_UMAP_dim50@cell.embeddings)
clist[[zoom_type]] <- cello
saveRDS(clist, "VisCello/clist.rds") 



# remove doublets (PTPRC+ cluster) identified in viscello and recompute
doublet_ids <- read_csv("MSC_Doublets_Data_VisCello.csv")
MSCs <- subset(MSCs, invert = TRUE, cells = doublet_ids$...1)
# recompute umap
MSCs <- NormalizeData(MSCs)
MSCs <- FindVariableFeatures(MSCs)
MSCs <- ScaleData(MSCs)
MSCs <- RunPCA(MSCs, reduction.name = "MSC_pca", features = VariableFeatures(MSCs))
MSCs <- RunUMAP(MSCs, dims = 1:10, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim10")
MSCs <- RunUMAP(MSCs, dims = 1:20, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim20")
MSCs <- RunUMAP(MSCs, dims = 1:30, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim30")
MSCs <- RunUMAP(MSCs, dims = 1:40, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim40")
MSCs <- RunUMAP(MSCs, dims = 1:50, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim50")
MSCs <- FindNeighbors(MSCs, dims = 1:20)
MSCs <- FindClusters(MSCs, algorithm = 2, resolution = 0.25, group.singletons = FALSE)

# Repeat for Endothelial Cells ----
# Subset Endothelial Cells, Re-Normalize/Scale, Compute PCA and UMAP ----
Endo <- subset(combined, cluster_anno_l2 == "AEC" | cluster_anno_l2 == "SEC")
Endo <- NormalizeData(Endo)
Endo <- FindVariableFeatures(Endo)
Endo <- ScaleData(Endo)
Endo <- RunPCA(Endo, reduction.name = "Endo_pca", features = VariableFeatures(Endo))
Endo <- RunUMAP(Endo, dims = 1:10, reduction = "Endo_pca", reduction.name = "Endo_UMAP_dim10")
Endo <- RunUMAP(Endo, dims = 1:20, reduction = "Endo_pca", reduction.name = "Endo_UMAP_dim20")
Endo <- RunUMAP(Endo, dims = 1:30, reduction = "Endo_pca", reduction.name = "Endo_UMAP_dim30")
Endo <- RunUMAP(Endo, dims = 1:40, reduction = "Endo_pca", reduction.name = "Endo_UMAP_dim40")
Endo <- RunUMAP(Endo, dims = 1:50, reduction = "Endo_pca", reduction.name = "Endo_UMAP_dim50")
Endo <- FindNeighbors(Endo, dims = 1:20)
Endo <- FindClusters(Endo, algorithm = 2, resolution = 0.25, group.singletons = FALSE)
Endo_Markers <- FindAllMarkers(Endo)
write.csv(Endo_Markers, "Endo_Markers.csv")

DimPlot(Endo, label = TRUE, reduction = "Endo_UMAP_dim50") + coord_fixed() + NoAxes() + NoLegend()
VlnPlot(Endo, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), pt.size = 0, group.by = 'seurat_clusters') # note high nCount and nFeature for clusters 4 and 5, extremely low for Cluster 3
FeaturePlot(Endo, reduction = "Endo_UMAP_dim50", features= c("CDH5", "PTPRC", "CD34", "DCN", "CXCL12", "KITLG"))

# Add the new zoom of Endo with putative doublets removed to Viscello
zoom_type <- "Endo_Trajectory"
cur_idx <-  which(colnames(eset) %in% colnames(Endo))
cello <- new("Cello", name = zoom_type, idx = cur_idx) 
cello@proj <- list('Endo_PCA' = Endo@reductions$Endo_pca@cell.embeddings,
                   'Endo_UMAP_dim10'=Endo@reductions$Endo_UMAP_dim10@cell.embeddings, 
                   'Endo_UMAP_dim20'=Endo@reductions$Endo_UMAP_dim20@cell.embeddings, 
                   'Endo_UMAP_dim30'=Endo@reductions$Endo_UMAP_dim30@cell.embeddings, 
                   'Endo_UMAP_dim40'=Endo@reductions$Endo_UMAP_dim40@cell.embeddings,
                   'Endo_UMAP_dim50'=Endo@reductions$Endo_UMAP_dim50@cell.embeddings)
clist[[zoom_type]] <- cello
saveRDS(clist, "VisCello/clist.rds") 

# remove doublets (PTPRC+ and DCN+ cluster) identified in viscello and recompute
Endo_doublet_ids <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/Final_scRNA_Analysis/Endo_Doublets_LQcells_Removed_Data_VisCello.csv")
Endo <- subset(Endo,  cells = Endo_doublet_ids$...1) # remove doublets
Endo <- subset(Endo, subset = seurat_clusters != 3) # remove likely low quality cells

# recompute umap and clustering
Endo <- NormalizeData(Endo)
Endo <- FindVariableFeatures(Endo)
Endo <- ScaleData(Endo)
Endo <- RunPCA(Endo, reduction.name = "Endo_pca", features = VariableFeatures(Endo))
Endo <- RunUMAP(Endo, dims = 1:10, reduction = "Endo_pca", reduction.name = "Endo_UMAP_dim10")
Endo <- RunUMAP(Endo, dims = 1:20, reduction = "Endo_pca", reduction.name = "Endo_UMAP_dim20")
Endo <- RunUMAP(Endo, dims = 1:30, reduction = "Endo_pca", reduction.name = "Endo_UMAP_dim30")
Endo <- RunUMAP(Endo, dims = 1:40, reduction = "Endo_pca", reduction.name = "Endo_UMAP_dim40")
Endo <- RunUMAP(Endo, dims = 1:50, reduction = "Endo_pca", reduction.name = "Endo_UMAP_dim50")
Endo <- FindNeighbors(Endo, dims = 1:20)
Endo <- FindClusters(Endo, algorithm = 2, resolution = 1, group.singletons = FALSE)

# Figure 2A Make UMAP for MSC Subsets ----
p1 <- DimPlot(MSCs, label = TRUE, group.by = 'cluster_anno_l2',reduction = "MSC_UMAP_dim50", raster = TRUE, raster.dpi = c(4112,4112), pt.size = 8, cols = cal2_cols) + coord_fixed() + NoAxes() + NoLegend()
p1_raster <- rasterize(p1)
ggsave(p1_raster,filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Figure2/panel/PanelA_scRNA_MSC_UMAP.eps", device = 'eps', width = 5,height = 5)
ggsave(p1_raster,filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Figure2/panel/PanelA_scRNA_MSC_UMAP.pdf")


# Supplemental Figure 2A RNAlo MSC QC ----
VlnPlot(subset(combined, cluster_anno_l2 %in% c("Adipo-MSC", "Osteo-MSC", "THY1+ MSC", "OsteoFibro-MSC", "Osteoblast", "Fibro-MSC", "RNAlo MSC")), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, cols = cal2_cols)

# Supplemental Figure 2B - Various Canonical MSC genes ---- 

mesenchymal_cell_types <-  c("Adipo-MSC","THY1+ MSC","Fibro-MSC","OsteoFibro-MSC","Osteo-MSC","Osteoblast", "SEC", "AEC", "VSMC") # need to reload the original untransformed combined object
mes <- subset(combined, cluster_anno_l2 %in% mesenchymal_cell_types)
mes$cluster_anno_l2 <- factor(mes$cluster_anno_l2, levels =c("Fibro-MSC","OsteoFibro-MSC","Osteo-MSC","Osteoblast","Adipo-MSC","THY1+ MSC","AEC", "SEC", "VSMC"))
p1 <- DotPlot(mes,group.by = "cluster_anno_l2", scale = FALSE, features = c("NES","LBP","NCAM1", "SPP1", "DCN", "CSPG4"), cluster.idents = FALSE) & RotatedAxis() & scale_color_viridis(option = "plasma") & coord_flip()
ggsave(p1, filename = '~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure2/PanelS2B_DotPlot_ExtraMarkerGenes_v1.pdf', width = 8, height = 3.5)



DotPlot(subset(combined, cluster_anno_l2 %in% mesenchymal_cell_types),scale = FALSE, cols = "Inferno", features = c("NES","LBP","NCAM1", "SPP1", "DCN", "CSPG4")) & RotatedAxis() 
p1 <- VlnPlot(subset(combined, cluster_anno_l2 %in% mesenchymal_cell_types),slot = 'data' , pt.size = 0,  features = "NES", cols = cal2_cols) & RotatedAxis() & NoLegend()
ggsave(p1, filename = '~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure2/PanelB_NES_VlnPlots.pdf',
       device = 'pdf', width = 10, height = 3.5)

# Supplemental Figure 2C - Lymphatic Endothelial cell markers ----
p1 <- FeaturePlot(Endo, reduction = "Endo_UMAP_dim30", features = c("LYVE1"), pt.size = 0.1, max.cutoff = 'q99', cols=brewer.pal(n = 100, name = "Reds"), coord.fixed = TRUE, ncol = 3) & NoAxes()
p1_raster <- rasterize(p1, dpi = 300)
ggsave(p1_raster, filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure2/PanelC_LEC_markers_LYVE1.pdf", device = "pdf", height = 5, width = 7)
p2 <- FeaturePlot(Endo, reduction = "Endo_UMAP_dim30", features = c("PROX1"), pt.size = 0.1, max.cutoff = 'q99', cols=brewer.pal(n = 100, name = "Reds"), coord.fixed = TRUE, ncol = 3) & NoAxes()
p2_raster <- rasterize(p2, dpi = 300)
ggsave(p2_raster, filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure2/PanelC_LEC_markers_PROX1.pdf", device = "pdf", height = 5, width = 7)
p3 <- FeaturePlot(Endo, reduction = "Endo_UMAP_dim30", features = c("PDPN"), pt.size = 0.1, max.cutoff = 'q99', cols=brewer.pal(n = 100, name = "Reds"), coord.fixed = TRUE, ncol = 3) & NoAxes()
p3_raster <- rasterize(p3, dpi = 300)
ggsave(p3_raster, filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure2/PanelC_LEC_markers_PDPN.pdf", device = "pdf", height = 5, width = 7)

 
FeaturePlot(MSCs, reduction = "MSC_UMAP_dim30", features = c("LYVE1", "PROX1", "PDPN"), pt.size = 0.1, max.cutoff = 'q99', coord.fixed = TRUE, cols=brewer.pal(n = 100, name = "Reds")) & NoAxes()

# Supplemental Figure 2E - CSF3 -----
p1 <- VlnPlot(combined, features = c("CSF3"), pt.size = 1,  sort = TRUE) + NoLegend()
ggsave(p1, filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure2/PanelD_GCSF_markers.pdf", device = "pdf", height = 5, width = 12)

p1 <- FeaturePlot(combined, max.cutoff = 'q99', features = c("CSF3"), coord.fixed = TRUE, order = TRUE, cols=brewer.pal(n = 100, name = "Reds")) & NoAxes() # note use of order=TRUE to highlight the positive cells
p1_raster <- rasterize(p1, dpi = 300)
ggsave(p1_raster, filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure2/PanelD_CSF3_UMAP.pdf", device = "pdf", height = 5, width = 7)

# Figure 2B - Dot Plot showing key DEGs between MSC subsets -----
DotPlot(MSCs, group.by = 'cluster_anno_l2', features = c("NES", "CSPG4", "LEPR", "COL1A1", "WIF1", "SPP1", "MGP", "CDH2"), scale = FALSE, cluster.idents = TRUE) + scale_color_viridis(option = 'plasma') + RotatedAxis()
ggsave(p1, filename = '~/Documents/Manuscripts/NBM_Atlas/Figures/Figure2/panel/PanelB_DotPlot_MCS_DEGs.pdf',
       device = 'pdf', width = 10, height = 3.5)
VlnPlot(MSCs, group.by = 'cluster_anno_l2', features = c("NES", "CSPG4", "LEPR", "COL1A1", "WIF1", "SPP1", "MGP", "CDH2"), scale = FALSE, cluster.idents = TRUE) + scale_color_viridis(option = 'plasma') + RotatedAxis()

# Not a figure - but sorting strategy justification -----
DotPlot(MSCs, group.by = 'cluster_anno_l2', features = c("PTPRC","GYPA", "CD38","CDH5", "PDPN", "NCAM1", "LEPR", "THY1"), scale = FALSE, cluster.idents = TRUE) + scale_color_viridis(option = 'plasma') + RotatedAxis()

# Supplemental Figure 2D MSC and Endo frequency per sample ----
MSCs$cluster_anno_l2 <- droplevels(MSCs$cluster_anno_l2)
cell_counts_MSC <- as.data.frame(table(MSCs$cluster_anno_l2,MSCs$orig.ident))

cell_counts_MSC$Var2 <- droplevels(cell_counts_MSC$Var1)
MSC_cols <- MetBrewer::met.brewer("Austria", n = 6)
MSC_cols <- RColorBrewer::brewer.pal(n = 6, name = "Set2")


p3 <- ggplot(data=cell_counts_MSC, aes(x=forcats::fct_rev(Var2),y=Freq, fill=Var1)) +
  geom_bar(position="fill",stat="identity") + 
  coord_flip() +
  theme_bw() + 
  labs(y='Cell Type Frequency') + 
  labs(x='') + scale_fill_manual(values=MSC_cols) + theme(axis.text = element_text(size=16)) + ggtitle("CODEX Cell Type Frequencies Per Sample") + scale_x_discrete(labels = c("H41", "H39", "H38", "H36", "H35", "H34", "H33", "H32", "H24", "H23", "H21", "H14") )

Endo$cluster_anno_l2 <- droplevels(Endo$cluster_anno_l2)
cell_counts_Endo <- as.data.frame(table(Endo$cluster_anno_l2,Endo$orig.ident))



p4 <- ggplot(data=cell_counts_Endo, aes(x=forcats::fct_rev(Var2),y=Freq, fill=Var1)) +
  geom_bar(position="fill",stat="identity") + 
  coord_flip() +
  theme_bw() + 
  labs(y='Cell Type Frequency') + 
  labs(x='') + scale_fill_manual(values= c("#AA336A", "#93E9BE")) + theme(axis.text = element_text(size=16)) + ggtitle("CODEX Cell Type Frequencies Per Sample") + scale_x_discrete(labels = c("H41", "H39", "H38", "H36", "H35", "H34", "H33", "H32", "H24", "H23", "H21", "H14") + RotatedAxis())

ggsave(p3, device = "pdf", height = 3, width = 4, units = "in", filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure2/MSC_Frequency.pdf")
ggsave(p4, device = "pdf", height = 3, width = 4, units = "in", filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure2/Endo_Frequency.pdf")

## calculate relative frequency of each cell type
cell_counts_MSC %>% dplyr::group_by(Var1) %>% summarise(rel_freq = Freq/n())

cell_counts_MSC <- cell_counts_MSC %>%
  group_by(Var2) %>%
  mutate(total = sum(Freq))

cell_counts_MSC$rel_freq <-  cell_counts_MSC$Freq / cell_counts_MSC$total
cell_counts_MSC_summarystats <-  cell_counts_MSC %>% group_by(Var1) %>% summarise(median=median(rel_freq), std.dev = sd(rel_freq), min = min(rel_freq), max = max(rel_freq))
cell_counts_MSC %>% group_by(Var2) %>% summarise(sum(rel_freq))

cell_counts_MSC$is_thy1oradipo <- cell_counts_MSC$Var1 == "THY1+ MSC" | cell_counts_MSC$Var1 == "Adipo-MSC"
cell_counts_MSC_summarystats <-  cell_counts_MSC %>% group_by(is_thy1oradipo) %>% summarise(median=median(rel_freq), std.dev = sd(rel_freq), min = min(rel_freq), max = max(rel_freq))

cell_counts_Endo <- cell_counts_Endo %>%
  group_by(Var2) %>%
  mutate(total = sum(Freq))

cell_counts_Endo$rel_freq <-  cell_counts_Endo$Freq / cell_counts_Endo$total
cell_counts_Endo_summarystats <-  cell_counts_Endo %>% group_by(Var1) %>% summarise(median=median(rel_freq), std.dev = sd(rel_freq), min = min(rel_freq), max = max(rel_freq))
# test that SEC is more than AEC
prop.test(7,12, p=0.5, correct = FALSE)
# Supplemental Table 2 - calculate DEGs between MSC and Endo subsets ----
MSCs<- SetIdent(MSCs, value = "cluster_anno_l2")
MSC_markers <- FindAllMarkers(MSCs)
Endo_markers <- FindAllMarkers(Endo)

# Figure 2B - Dot Plot of canonical mesenchymal genes ----
p1 <- DotPlot(MSCs, group.by = 'cluster_anno_l2', features = c("PDGFRA", "VIM","CXCL12", "LEPR", "LPL", "APOE","THY1","RUNX2", "SP7", "IBSP", "BGLAP", "COL1A1","GSN", "APOD", "PDPN", "HAS1","DPT", "CD164", "NT5E"), scale = FALSE, cluster.idents = FALSE) + scale_color_viridis(option = 'plasma') + RotatedAxis()
ggsave(p1, filename = '~/Documents/Manuscripts/NBM_Atlas/Figures/Figure2/panel/PanelB_DotPlot_MCS_DEGs.eps',
       device = 'pdf', width = 10, height = 3.5)
# Figure 2C Dot Plot of ISCT genes and commonly used MSC markers ----
p1 <- DotPlot(MSCs, group.by = 'cluster_anno_l2', features = c("NGFR", "MCAM", "NT5E", "THY1", "ENG"), scale = FALSE, cluster.idents = FALSE) + scale_color_viridis(option = 'plasma') + RotatedAxis()
ggsave(p1, filename = '~/Documents/Manuscripts/NBM_Atlas/Figures/Figure2/panel/PanelB_DotPlot_MCS_DEGs.eps',
       device = 'pdf', width = 10, height = 3.5)

# Figure 2D - CytoTRACE Analysis  ----
## See Cytotrace script, run by Jon Sussman 
MSCs_CT <- readRDS("~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/Final_scRNA_Analysis/CytoTRACE/MSCs_CytoTRACE_computed.RDS")
p1 <- FeaturePlot(MSCs_CT, features = "CytoTRACE_score", reduction = "MSC_UMAP_dim50", pt.size = 1,coord.fixed = TRUE) + NoAxes() + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p1_raster <- rasterize(p1, dpi = 300)
ggsave(p1_raster, device = "pdf", filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Figure2/panel/CytoTRACE_Score_FeaturePlot.pdf", units = "in", width = 5, height = 5)

# Figure 2D Make dot plot of DEGs between SEC and AEC -----
Endo_markers <- FindMarkers(Endo, group.by = 'cluster_anno_l2', ident.1 = "AEC", ident.2 = "SEC")
Endo_markers$gene_name <- rownames(Endo_markers)
Endo_markers%>% 
  top_n(n = 10, wt = avg_log2FC) -> top10_AEC # Get 10 most significant genes for each cell type
Endo_markers  %>% 
  top_n(n = 10, wt = -avg_log2FC) -> top10_SEC # Get 10 most significant genes for each cell type
# note chose both some canonical markers and top DEGs
p1 <- DotPlot(Endo, group.by = 'cluster_anno_l2', features = c("CDH5", "CD34", "KDR","PODXL", "ICAM2", rownames(top10_AEC), "NR2F2", "EPHB4", rownames(top10_SEC)), scale = FALSE,  cluster.idents = TRUE) + scale_color_viridis(option = 'plasma') + RotatedAxis()
ggsave(p1, filename = '~/Documents/Manuscripts/NBM_Atlas/Figures/Figure2/panel/PanelD_DotPlot_BMEC_DEGs.pdf', width = 10, height = 3.5)

# Reference Map Other MSC Datasets to Our Reference -----
# Figure 2F - Reference Map DeJong MSCs ----- 
# data first published https://pubmed.ncbi.nlm.nih.gov/34017122/) 
MSCs <- RunUMAP(MSCs, dims = 1:50, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim50", return.model = TRUE) # rerun find umap to return the umap model

DeJong_MSCs <- readRDS(file = "~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/Final_scRNA_Analysis/Objects_To_RefMap/DeJong_MSCs.RDS")
anchors <- FindTransferAnchors(
  reference = MSCs,
  query = DeJong_MSCs,
  normalization.method = "LogNormalize",
  reference.reduction = "MSC_pca",
  dims = 1:50
)
DeJong_MSCs <- MapQuery(
  anchorset = anchors,
  query = DeJong_MSCs,
  reference = MSCs,
  refdata = list(
    MSC_refmap = "cluster_anno_l2"
  ),
  reference.reduction = "MSC_pca", 
  reduction.model = "MSC_UMAP_dim50"
)

DimPlot(DeJong_MSCs, reduction = "umap", group.by = "predicted.MSC_refmap", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
DimPlot(DeJong_MSCs, reduction = "ref.umap", group.by = "predicted.MSC_refmap", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

p1 <- ggplot(DeJong_MSCs@meta.data, aes(fill=predicted.MSC_refmap, x=orig.ident)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force) # limits=force drops unused levels from the legend

ggsave(p1, filename = '~/Documents/Manuscripts/NBM_Atlas/Figures/Figure2/panel/PanelF_StackedBar_DeJongMSCs.pdf',width = 5, height = 5)

# Repeat for Fetal BM MSCs
# Figure 2F Fetal MSC Reference Mapping ------
# With fetal bone marrow
Convert("fig1b_fbm_scaled_gex_updated_dr_20210104.h5ad", dest = "h5seurat", overwrite = TRUE) # Convert h5ad obj from Jardine paper to seurat
fbm <- LoadH5Seurat("~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/Final_scRNA_Analysis/Objects_To_RefMap/fig1b_fbm_scaled_gex_updated_dr_20210104.h5seurat")
fbm_stroma <- subset(fbm, broad_fig1_cell.labels == "stroma")
fbm_stroma_MSC <- subset(fbm_stroma, cell.labels == "adipo-CAR" | cell.labels == "osteoblast" | cell.labels == "osteoblast precursor" | cell.labels == "osteochondral precursor" | cell.labels =="endosteal fibroblast" | cell.labels =="arteriolar fibroblast") # chondrocytes and muscle/muscle stem cells not included because they aren't represented in our dataset (not in adult BM)

remove(fbm)
anchors <- FindTransferAnchors(
  reference = MSCs,
  query = fbm_stroma_MSC,
  normalization.method = "LogNormalize",
  reference.reduction = "MSC_pca",
  dims = 1:50
)
fbm_stroma_MSC <- MapQuery(
  anchorset = anchors,
  query = fbm_stroma_MSC,
  reference = MSCs,
  refdata = list(
    MSC_refmap = "cluster_anno_l2"
  ),
  reference.reduction = "MSC_pca", 
  reduction.model = "MSC_UMAP_dim50"
)

DimPlot(fbm_stroma_MSC, reduction = "umap", group.by = "predicted.MSC_refmap", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 <- ggplot(fbm_stroma_MSC@meta.data, aes(fill=predicted.MSC_refmap, x=orig.ident)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force)

plot_grid(p1+p2)
ggsave(p1, filename = '~/Documents/Manuscripts/NBM_Atlas/Figures/Figure2/panel/PanelG_StackedBar_Fetal_MSCs.pdf',width = 5, height = 5)


# Figure 3A - Plot Notable Putative Hematopoietic Supportive Factors ----
cytokines <- c("TGFB1","TNFRSF11B","CDH2","SPP1", "DLL1", "DLL4", "JAG1", "JAG2","SELE", "SELP","SELPLG","FLT3LG", "CSF1", "CXCL12", "KITLG","VCAM1", "TNFSF11", "IL6","IL7","PTN", "FGF1", "FGF2","NGF", "ANGPT1", "ANGPT2", "IGF1")
mesenchymal_cell_types <-  c("Adipo-MSC","THY1+ MSC","Fibro-MSC","OsteoFibro-MSC","Osteo-MSC","Osteoblast", "SEC", "AEC", "VSMC") # need to reload the original untransformed combined object
mes <- subset(combined, cluster_anno_l2 %in% mesenchymal_cell_types)
mes$cluster_anno_l2 <- factor(mes$cluster_anno_l2, levels =c("Fibro-MSC","OsteoFibro-MSC","Osteo-MSC","Osteoblast","Adipo-MSC","THY1+ MSC","AEC", "SEC", "VSMC"))
p1 <- DotPlot(mes,group.by = "cluster_anno_l2", scale = FALSE,  features = cytokines, cluster.idents = FALSE) & RotatedAxis() & scale_color_viridis(option = "plasma")
ggsave(p1, filename = '~/Documents/Manuscripts/NBM_Atlas/Figures/Figure3/panels/PanelA_DotPlot_SupportiveFactors_v3.pdf', width = 10, height = 3.5)



# Figure 5F - AUCell Hypoxia Score Analysis ----
myeloid_development <- c("HSC", "MPP","Cycling HSPC", "GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil")
FeaturePlot(combined_AUC, reduction = "UMAP_dim30", features = "HALLMARK_HYPOXIA", coord.fixed = TRUE, max.cutoff = 'q99') + 
  scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + NoAxes()

VlnPlot(subset(combined_AUC, subset = cluster_anno_l2 %in% myeloid_development), group.by = 'cluster_anno_l2', features = "HALLMARK_HYPOXIA", pt.size = 0, sort = TRUE, cols = cal2_cols) + NoLegend()

combined_AUC_myeloid <- combined_AUC %>% subset(cluster_anno_l2 %in% myeloid_development)

# Perform one vs. rest T.test with multiple hypothesis correction for each population
test <- c()
c <- c()
d <- c()
for (i in levels(droplevels(combined_AUC_myeloid$cluster_anno_l2))) {
  a <- subset(combined_AUC_myeloid@meta.data, cluster_anno_l2 == {i})$HALLMARK_HYPOXIA
  d <- append(x = d, values = i)
  b <- subset(combined_AUC_myeloid@meta.data, cluster_anno_l2 != {i})$HALLMARK_HYPOXIA
  test <- t.test(a,b, alternative = "less")
  c <- append(x = c, values = test$p.value)
}
#correct for multiple hypothesis 
c <- p.adjust(c,method = "BH")
names(c) <- d
c

subset(combined_AUC_myeloid@meta.data, cluster_anno_l2 == "Neutrophil")$HALLMARK_HYPOXIA

# Supplemental Figure S5B Distance to Bone ----
cn_ranks <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Non-Cell Microenvironment Analysis/combined_neighbor.csv")

nbs_names <- c("HSC / Mature Myeloid", "Erythroid/Myeloid", "PC/Arteriolar", "Erythroid", "Arteriolar", "Erythroid", "Lymphoid", "Erythroid/Myeloid/Lymphoid", "Early Myeloid / Endosteal", "Myeloid/Lymphoid", "HSPC/Intermediate Myeloid", "Erythroid/Myeloid/Lymphoid", "Erythroid/Myeloid", "Early Myeloid / Arteriolar", "Peri-Arterolar Lymphoid")
nbs_cols <- c("#A8D37F", "#51C8EB", "#FF00FF", "#FBED24", "#FF0000", "#FBED24", "#89919C", "#BD7CB5", "#4B409A", "#FFA500", "#A15A26", "#BD7CB5", "#51C8EB", "#4B409A", "#FF007F")
names(nbs_cols) <- nbs_names

p2 <- cn_ranks %>% dplyr::filter(cn_ranks$structure == "bone") %>% ggplot(aes(x = reorder(neighborhoods, -normalized_rank), y = normalized_rank, fill = neighborhoods)) + geom_boxplot() + theme_minimal() + NoLegend()+ RotatedAxis() + scale_fill_manual(values = nbs_cols)
ggsave(p2, filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure5/neighborhoods_bonedist_normrank.pdf", device = "pdf", units = "in", width = 8, height = 4)



