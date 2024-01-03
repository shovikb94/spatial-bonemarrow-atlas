library(Seurat)
library(ggplot2)
library(rhdf5)

setwd("/mnt/isilon/tan_lab/sussmanj/Temp/BoneMarrow_scRNA/Datasets")
data <- h5read("Ennis_Counts.h5", "data")
cells = read.table("Ennis_Cells.txt", sep = '\t')
genes = read.table("Ennis_Genes.txt")

rownames(data) <- genes$V1
colnames(data) <- cells$V1
metacolumns <- c("dataset", "sample", "donor_status", "donor", "timepoint", "celltype")

seurat_assay <- CreateAssayObject(data = data, assay = "RNA")
metadata <- h5read("Ennis_Metadata.h5", "metadata")
rownames(metadata) <- metacolumns
colnames(metadata) <- cells$V1

umap <- h5read("Ennis_UMAP.h5", "umap_coordinates")
rownames(umap) <- c("UMAP1", "UMAP2")
umap <- as(t(umap), "matrix")
rownames(umap) <- cells$V1

seurat <- CreateSeuratObject(seurat_assay, project = "EnnisAtlas")
seurat <- AddMetaData(seurat, metadata = t(metadata))
seurat[["umap_original"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_")

DimPlot(seurat, group.by = "celltype", raster = F, alpha = 0.2) + coord_fixed()
FeaturePlot(seurat, features = "PTPRC", raster = F, alpha = 0.2, max.cutoff = 'q95') + coord_fixed()                      
saveRDS(seurat, "Ennis_Atlas_Seurat.rds")
