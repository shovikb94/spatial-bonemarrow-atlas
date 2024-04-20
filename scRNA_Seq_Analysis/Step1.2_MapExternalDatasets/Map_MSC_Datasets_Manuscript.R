# Script for Mapping all MSC Subsets

# External Datasets Reference Mapping


library(DoubletFinder)
library(cowplot)
library(Seurat)
library(SeuratDisk)
library(ggplot2)

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

# save cluster_anno_l2 colors
cal2_cols <- c("#FFD580", "#FFBF00", "#B0DFE6", "#7DB954", "#64864A", "#8FCACA", "#4682B4", "#CAA7DD", "#B6D0E2", "#a15891", "#FED7C3", "#A8A2D2", "#CF9FFF", "#9C58A1", "#2874A6", "#96C5D7", "#63BA97", "#BF40BF", "#953553", "#6495ED", "#E7C7DC", "#5599C8", "#FA8072", "#F3B0C3", "#F89880", "#40B5AD", "#019477", "#97C1A9", "#C6DBDA", "#CCE2CB", "#79127F", "#FFC5BF", "#ff9d5c", "#FFC8A2", "#DD3F4E")
cal2_col_names <- c("Adipo-MSC", "AEC", "Ba/Eo/Ma", "CD4+ T-Cell", "CD8+ T-Cell", "CLP", "Cycling DCs", "Cycling HSPC", "Early Myeloid Progenitor", "Erythroblast", "Fibro-MSC", "GMP", "HSC", "Late Erythroid", "Late Myeloid", "Macrophages", "Mature B", "Megakaryocyte", "MEP", "Monocyte", "MPP", "Neutrophil", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "pDC", "Plasma Cell", "Pre-B", "Pre-Pro B", "Pro-B", "RBC", "RNAlo MSC", "SEC", "THY1+ MSC", "VSMC")
names(cal2_cols) <- cal2_col_names


library(Seurat)
library(dplyr)

# import Atlas Obj and MSC-subsetted Atlas
combined <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Seurat/SB66_combined_corrected.RDS")
MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Seurat/MSCs_Final_091523.RDS")

# filter out rna lo msc and likely doublet so-called NKT population
combined <- subset(combined, subset = cluster_anno_l2 != "RNAlo MSC" & cluster_anno_l2 != "NKT Cell")
combined <- RunUMAP(combined, dims = 1:30, reduction = "pca", reduction.name = "UMAP_dim30", return.model = TRUE) # rerun find umap to return the umap model

# Reference Map Other MSC Datasets to Our Reference -----

MSCs <- RunUMAP(MSCs, dims = 1:50, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim50", return.model = TRUE) # rerun find umap to return the umap model


# Figure 2I - Reference Map DeJong MSCs ----- 

DeJong_MSCs <- readRDS(file = "~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/Final_scRNA_Analysis/Objects_To_RefMap/DeJong_MSCs.RDS")
anchors <- FindTransferAnchors(
  reference = combined,
  query = DeJong_MSCs,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

DeJong_MSCs <- MapQuery(
  anchorset = anchors,
  query = DeJong_MSCs,
  reference = combined,
  refdata = list(
    MSC_refmap = "cluster_anno_l1"
  ),
  reference.reduction = "pca", 
  reduction.model = "UMAP_dim50"
)

DeJong_MSCs$predicted_cluster_anno_l1 <- DeJong_MSCs$predicted.combined_refmap
DeJong_MSCs <- subset(x=DeJong_MSCs, subset = predicted_cluster_anno_l1 == "Mesenchymal")

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



# Figure 2I Fetal MSC Reference Mapping (Supplemental Figure S4B-C) ------
# With fetal bone marrow
library(SeuratDisk)
Convert("fig1b_fbm_scaled_gex_updated_dr_20210104.h5ad", dest = "h5seurat", overwrite = TRUE) # Convert h5ad obj from Jardine paper to seurat
fbm <- LoadH5Seurat("~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/Final_scRNA_Analysis/Objects_To_RefMap/fig1b_fbm_scaled_gex_updated_dr_20210104.h5seurat")
fbm <- LoadH5Seurat("/mnt/isilon/tan_lab/bandyopads/SB_CHOP_Comp_Backup/110923/NBM_Microenvironment/NBM_Atlas_scRNA/Final_scRNA_Analysis/Objects_To_RefMap/fig1b_fbm_scaled_gex_updated_dr_20210104.h5seurat")

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
p1 <- ggplot(fbm_stroma_MSC@meta.data, aes(fill=predicted.MSC_refmap, x=cell.labels)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force)

p2 <- ggplot(fbm_stroma_MSC@meta.data, aes(fill=predicted.MSC_refmap, x=orig.ident)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force)

plot_grid(p1+p2)
ggsave(plot_grid(p1+p2), filename = '/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/figures/PanelG_StackedBar_Fetal_MSCs.pdf',width = 7, height = 3)

# Map Triana MSCs (Figure 2I) ----

Triana_FullData <- readRDS(file = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/data/Triana_2021_NatCellBio/Healthy.rds")
table(levels(as.factor(Triana_FullData$CellTypes)))
Triana_MSCs <- subset(Triana_FullData, subset = CellTypes == "Mesenchymal cells_1" | CellTypes == "Mesenchymal cells_2")
remove(Triana_FullData)

anchors <- FindTransferAnchors(
  reference = MSCs,
  query = Triana_MSCs,
  normalization.method = "LogNormalize",
  reference.reduction = "MSC_pca", k.filter = 10, k.score = 17,
  dims = 1:50
)
Triana_MSCs <- MapQuery(
  transferdata.args = list(k.weight = 30), 
  anchorset = anchors,
  query = Triana_MSCs, 
  reference = MSCs,
  refdata = list(
    MSC_refmap = "cluster_anno_l2"
  ),
  reference.reduction = "MSC_pca", 
  reduction.model = "MSC_UMAP_dim50",
)

p1 <- ggplot(Triana_MSCs@meta.data, aes(fill=predicted.MSC_refmap, x=orig.ident)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force) # limits=force drops unused levels from the legend
DotPlot(Triana_MSCs, features = c("NES", "CSPG4", "LEPR", "COL1A1", "WIF1", "SPP1", "MGP", "CDH2"), scale = FALSE, cluster.idents = TRUE) + scale_color_viridis(option = 'plasma') + RotatedAxis()



# Map eLife Li MSCs (Figure 2I) ----

# Load the Matrix package
library(Matrix)

# Read the two .mtx files
#matrix1 <- readMM("ExternalDataset_Mapping/data/Li_2023_eLife/GSE190965_spliced.mtx.gz")
#matrix2 <- readMM("ExternalDataset_Mapping/data/Li_2023_eLife/GSE190965_unspliced.mtx.gz")

# Check if dimensions are the same
#if (dim(matrix1) != dim(matrix2)) {
#  stop("Matrices have different dimensions and cannot be summed.")
#}


# Sum the two matrices
#result_matrix <- matrix1 + matrix2
#result_matrix <- t(result_matrix)
# Save the result matrix as a .mtx file
#writeMM(result_matrix, "ExternalDataset_Mapping/data/Li_2023_eLife/matrix.mtx")

# read in 10x file
Li_FullData <- Read10X(data.dir = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/data/Li_2023_eLife/", gene.column = 2)

# create seurat obj and process
Li_FullData <- CreateSeuratObject(counts = Li_FullData, project = "Li_2023_data", min.cells = 3, min.features = 100)
Li_FullData

# add cluster information and apply their processing
li_md <- read_csv("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/data/Li_2023_eLife/Hongzhe_clustering_and_2dUMAP_louvain11_merged_95_renamed.csv")
rownames(li_md) <- li_md$...1

# 
Li_FullData <- AddMetaData(Li_FullData, metadata = li_md)
Li_Filtered <- subset(x = Li_FullData, subset = louvain11_merged.95_renamed != 'NA')  

ob.list <- list(Li_Filtered)
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

Li_Filtered <- ob.list[[1]]
VlnPlot(Li_Filtered, group.by = "louvain11_merged.95_renamed", features = c("PTPRC"), pt.size = 0) # confirm the cluster labels match the manuscript

# map using clusters annotated by authors
Li_MSCs <- subset(Li_Filtered, louvain11_merged.95_renamed %in% c(3,5,16,38,29,23,6,37,8))

anchors <- FindTransferAnchors(
  reference = MSCs,
  query = Li_MSCs,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

Li_MSCs<- MapQuery(
  anchorset = anchors,
  query = Li_MSCs,
  reference = MSCs,
  refdata = list(
    combined_refmap = "cluster_anno_l2"
  ),
  reference.reduction = "MSC_pca", 
  reduction.model = "MSC_UMAP_dim50"
)
Li_MSCs$predicted_cluster_anno_l2 <- Li_MSCs$predicted.combined_refmap


p1 <- ggplot(data = Li_MSCs@meta.data ,aes(fill=predicted_cluster_anno_l2, x = as.factor(louvain11_merged.95_renamed))) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force) # limits=force drops unused levels from the legend

table(Li_MSCs$louvain11_merged.95_renamed)
table(Li_MSCs$predicted_cluster_anno_l2)

# calculate percentage
6509/(6509+3+156+55+303)

DotPlot(Li_MSCs, group.by = "predicted.MSC_refmap", features = c("NES", "CSPG4", "LEPR", "COL1A1", "WIF1", "SPP1", "MGP", "CDH2"), scale = FALSE, cluster.idents = TRUE) + scale_color_viridis(option = 'plasma') + RotatedAxis()


# Map iScience Ennis MSCs (Figure 2I) ----

ennis_data <- readRDS("/mnt/isilon/tan_lab/sussmanj/Temp/BoneMarrow_scRNA/Datasets/Ennis_Atlas_Seurat.rds")
ennis_data@assays$RNA@counts <- ennis_data@assays$RNA@data

table(ennis_data$celltype)
ennis_MSCs <- subset(ennis_data, subset = celltype == "Mesenchymal")
remove(ennis_data)

ob.list <- list(ennis_MSCs)
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
  #ob.list[[i]] <- FindDoublets(ob.list[[i]], PCs = 1:30, sct = FALSE, exp_rate = (length(colnames(ob.list[[i]]))/125000))
}

ennis_MSCs <- ob.list[[1]]

anchors <- FindTransferAnchors(
  reference = MSCs,
  query = ennis_MSCs,
  normalization.method = "LogNormalize",
  reference.reduction = "MSC_pca", 
  dims = 1:50
)
ennis_MSCs <- MapQuery(
  anchorset = anchors,
  query = ennis_MSCs, 
  reference = MSCs,
  refdata = list(
    MSC_refmap = "cluster_anno_l2"
  ),
  reference.reduction = "MSC_pca", 
  reduction.model = "MSC_UMAP_dim50",
)

p1 <- ggplot(ennis_MSCs@meta.data, aes(fill=predicted.MSC_refmap, x=orig.ident)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force) # limits=force drops unused levels from the legend


table(ennis_MSCs$predicted.MSC_refmap)

DotPlot(ennis_MSCs, group.by = 'predicted.MSC_refmap', features = c("NES", "CSPG4", "LEPR", "COL1A1", "WIF1", "SPP1", "MGP", "CDH2"), scale = FALSE, cluster.idents = FALSE) + scale_color_viridis(option = 'plasma') + RotatedAxis()

#write.csv(ennis_data@assays$RNA@data, "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/data/Ennis_2023_iScience/ennis_msc_dataslot.csv")
#saveRDS(ennis_MSCs, "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/Ennis_MSCs.RDS")

# Map Int J Biol Sci Wang et al MSCs ----

Wang_oa <- Read10X(data.dir = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/data/Wang_2021_IntJBiolSci/GSE147287_RAW/oa_filtered_feature_bc_matrix/")
Wang_op <- Read10X(data.dir = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/data/Wang_2021_IntJBiolSci/GSE147287_RAW/op_filtered_feature_bc_matrix/")

# create seurat obj and process
Wang_oa <- CreateSeuratObject(counts = Wang_oa, project = "Wang_oa", min.cells = 3, min.features = 100)
Wang_oa
Wang_op <- CreateSeuratObject(counts = Wang_op, project = "Wang_op", min.cells = 3, min.features = 100)
Wang_op

ob.list <- list(Wang_oa, Wang_op)
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

Wang_oa <- ob.list[[1]]
Wang_op <- ob.list[[2]]



anchors <- FindTransferAnchors(
  reference = combined,
  query = Wang_oa,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)
Wang_oa <- MapQuery(
  anchorset = anchors,
  query = Wang_oa, 
  reference = combined,
  refdata = list(
    MSC_refmap = "cluster_anno_l2"
  ),
  reference.reduction = "pca", 
  reduction.model = "UMAP_dim30",
)

p1 <- ggplot(subset(Wang_oa@meta.data, subset = predicted.MSC_refmap %in% c("Adipo-MSC", "THY1+ MSC", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "Fibro-MSC")), aes(fill=predicted.MSC_refmap, x=orig.ident)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force) # limits=force drops unused levels from the legend

anchors <- FindTransferAnchors(
  reference = combined,
  query = Wang_op,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)
Wang_op <- MapQuery(
  anchorset = anchors,
  query = Wang_op, 
  reference = combined,
  refdata = list(
    MSC_refmap = "cluster_anno_l2"
  ),
  reference.reduction = "pca", 
  reduction.model = "UMAP_dim30",
)

p1 <- ggplot(subset(Wang_op@meta.data, subset = predicted.MSC_refmap %in% c("Adipo-MSC", "THY1+ MSC", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "Fibro-MSC")), aes(fill=predicted.MSC_refmap, x=orig.ident)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force) # limits=force drops unused levels from the legend

# calculate percentages
table(subset(Wang_op@meta.data, subset = predicted.MSC_refmap %in% c("Adipo-MSC", "THY1+ MSC", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "Fibro-MSC"))$predicted.MSC_refmap)
938/(938+118+7+52) # 0.84 adipo-msc percentage
table(subset(Wang_oa@meta.data, subset = predicted.MSC_refmap %in% c("Adipo-MSC", "THY1+ MSC", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "Fibro-MSC"))$predicted.MSC_refmap)
752/(752+1+4) #0.99 adipo-msc

(938+752) / (938+118+7+52 + 752 + 1 +4) # total dataset 0.90 adipo-MSC



