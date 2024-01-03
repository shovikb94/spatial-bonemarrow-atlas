# script to perform integrative analysis of fetal BM dataset Jardine et al and our data

library(Seurat)

fbm <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/fbm_MSCs.RDS")
MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Seurat/MSCs_Final_091523.RDS")
MSCs$Dataset <- "Bandyopadhyay_adult_BM"
MSCs$predicted.MSC_refmap <- MSCs$cluster_anno_l2
fbm$Dataset <- "Jardine_fetal_BM"

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- list(MSCs, fbm)

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
stroma.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(stroma.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
stroma.combined <- ScaleData(stroma.combined, verbose = FALSE)
stroma.combined <- RunPCA(stroma.combined, npcs = 30, verbose = FALSE)
stroma.combined <- RunUMAP(stroma.combined, reduction = "pca", dims = 1:30)
stroma.combined <- FindNeighbors(stroma.combined, reduction = "pca", dims = 1:30)
stroma.combined <- FindClusters(stroma.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(stroma.combined, reduction = "umap", group.by = "Dataset")
p2 <- DimPlot(stroma.combined, reduction = "umap", group.by = "predicted.MSC_refmap", label = TRUE,
              repel = TRUE)
p3 <- DimPlot(stroma.combined, reduction = "umap", group.by = "cell.labels", label = TRUE,
              repel = TRUE)
p1 + p2

stroma.combined <- SetIdent(stroma.combined, value = "Dataset")
integrated_markers_apodmsc <- FindMarkers(stroma.combined, ident.1 = "Bandyopadhyay_adult_BM", ident.2 = "Jardine_fetal_BM", group.by = "Dataset")      
integrated_markers_fibroMSC <- FindMarkers(subset(stroma.combined, subset = predicted.MSC_refmap == "Fibro-MSC"), ident.1 = "Bandyopadhyay_adult_BM", ident.2 = "Jardine_fetal_BM", group.by = "Dataset")      
write.csv(integrated_markers_fibroMSC, "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_FetalIntegration/integrated_markers_fibroMSC_bandyopadhyay_vs_jardine.csv")


VlnPlot(stroma.combined, features = "FOXD3", pt.size = 0, group.by = "predicted.MSC_refmap", split.by = "Dataset")

# Plot the composition of the fetal atlas with their labels comapred to ref mapped labels
ggplot(fbm@meta.data, aes(fill=predicted.MSC_refmap, x=cell.labels)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force) # limits=force drops unused levels from the legend


# Supplemental Figure S2H Gene Ontology Metascape Analysis ----
GO_Top10 <- read_csv("Top10_FibroMSC_FetalvsAdult_CellType_GO.csv")
GO_Top10$DescriptionID <- paste0(GO_Top10$Description,"_", GO_Top10$GO)
GO_Top10 <- GO_Top10 %>% mutate(Log10_qval_direction=replace(`Log10(q)`, Direction=="Fetal", -`Log10(q)`))
p1 <- GO_Top10 %>% ggplot(aes(x= reorder(DescriptionID, -Log10_qval_direction), y = Log10_qval_direction, fill = Direction)) + geom_bar(stat = "identity") + scale_fill_manual(values = c('#89CFF0', '#c1aed4')) + coord_flip() + theme_minimal()
ggsave(p1, filename = '/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_FetalIntegration/figures/FibroMSC_fetalvsadult_GO_metascape.pdf',width = 15, height = 5)

