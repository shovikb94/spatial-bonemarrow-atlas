
library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(ggrastr)

immune.filtered <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.filtered_FINAL.RDS")
# immune.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.combined_final_112322.RDS")
immune.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.combined_final_040523.RDS")


DefaultAssay(immune.filtered) <- "integrated"
DefaultAssay(immune.combined) <- "integrated"


immune.filtered <- RunUMAP(immune.filtered, dims = 1:30, reduction = "pca", reduction.name = "UMAP_dim30", return.model = TRUE)
immune.combined <- RunUMAP(immune.combined, dims = 1:30, reduction = "pca", reduction.name = "UMAP_dim30", return.model = TRUE)

# read in NSM objects
NSM1 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/NSM_seurat_objs/codex_NSM1_1720_SeuratObj.RDS")
NSM2 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/NSM_seurat_objs/codex_NSM2_1086_SeuratObj.RDS")
NSM3 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/NSM_seurat_objs/codex_NSM3_1996_SeuratObj.RDS")

ob.list <- list(NSM1,NSM2,NSM3)
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
NSM1 <- ob.list[[1]]
NSM2 <- ob.list[[2]]
NSM3 <- ob.list[[3]]

NSM.combined <- merge(x=NSM1, y = c(NSM2,NSM3), add.cell.ids = c("NSM1", "NSM2", "NSM3") )

anchors <- FindTransferAnchors(
  reference = immune.combined, k.filter = NA,
  query = NSM.combined,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:30
)
#saveRDS(anchors, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/anchors_unfiltered.RDS")
NSM.combined <- MapQuery(
  anchorset = anchors,
  query = NSM.combined,
  reference = immune.combined,
  refdata = list(
    MSC_refmap = "cluster_anno_l2"
  ),
  reference.reduction = "pca", 
  reduction.model = "UMAP_dim30"
)

#saveRDS(NSM.combined, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/NSM.combined.mapped_unfiltered.RDS")
NSM.combined$classified_cluster_anno_l2 <- NSM.combined$predicted.MSC_refmap
NSM.combined$classified_cluster_anno_l2_score <- NSM.combined$predicted.MSC_refmap.score

# anchors <- readRDS("ReferenceMap_NSM_Step5/anchors_unfiltered.RDS")
NSM.combined <- MapQuery(
anchorset = anchors,
query = NSM.combined,
reference = immune.combined,
refdata = list(
  MSC_refmap = "cluster_anno_coarse"
),
reference.reduction = "pca", 
reduction.model = "UMAP_dim30"
)
NSM.combined$classified_cluster_anno_coarse <- NSM.combined$predicted.MSC_refmap
NSM.combined$classified_cluster_anno_coarse_score <- NSM.combined$predicted.MSC_refmap.score
#saveRDS(NSM.combined, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/NSM.combined.mapped_unfiltered.RDS")






getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(as.factor(NSM.combined$predicted.MSC_refmap))))

DimPlot(NSM.combined, reduction = "ref.umap", group.by = "predicted.MSC_refmap", label = TRUE, label.size = 3, repel = TRUE, cols = cols) + NoLegend() + coord_fixed() + NoAxes()
ggplot(NSM.combined@meta.data, aes(fill=predicted.MSC_refmap, x=orig.ident)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample")

# Rename orig.idents 
levels(as.factor(NSM.combined$orig.ident))
orig.ident_annotations <- c("NSM_1996","NSM_1720", "NSM_1086")
NSM.combined <- SetIdent(NSM.combined, value = "orig.ident")
names(orig.ident_annotations) <- levels(as.factor(NSM.combined$orig.ident))
NSM.combined <- RenameIdents(NSM.combined, orig.ident_annotations)
NSM.combined@meta.data["Sample_Name"] <- NSM.combined@active.ident 



# Give coarse definition for combined object
levels(as.factor(NSM.combined$predicted.MSC_refmap))
coarse_annotations <- c("Adipocyte", "AEC","Artifact","Artifact", "Lymphoid", "HSPC", "Lymphoid","Artifact", "Lymphoid", "HSPC", "Myeloid", "Endosteal", "Erythroid", "Erythroid", "Megakaryocyte", "Megakaryocyte", 
                        "HSPC", "Myeloid", "HSPC", "HSPC", "Lymphoid", "Myeloid", "Macrophages","Myeloid", "HSPC", "Monocytes", "MSC", "Monocytes", "pDC", "Lymphoid", "SEC", "HSPC", "MSC","Artifact", "VSMC")
NSM.combined <- SetIdent(NSM.combined, value = "predicted.MSC_refmap")
names(coarse_annotations) <- levels(as.factor(NSM.combined$predicted.MSC_refmap))
NSM.combined <- RenameIdents(NSM.combined, coarse_annotations)
NSM.combined@meta.data["predicted_anno_coarse"] <- NSM.combined@active.ident 



cell_counts <- as.data.frame(table(NSM.combined$predicted.MSC_refmap,NSM.combined$orig.ident))
getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(cell_counts$Var1)))

cell_counts_coarse <- as.data.frame(table(NSM.combined$predicted_anno_coarse,NSM.combined$Sample_Name))
getPalette <- colorRampPalette(brewer.pal(15,"Paired"))
cols <- getPalette(length(levels(cell_counts_coarse$Var1)))
cols <- met.brewer(name="Signac", n=16)


ggplot(data=cell_counts, aes(x=reorder(Var2, Freq),y=Freq, fill=Var1)) +
  geom_bar(position="fill",stat="identity") + 
  coord_flip() +
  theme_bw() + 
  labs(x='Count') + 
  labs(y='Cell Type') + scale_fill_manual(values=cols) + theme(axis.text = element_text(size=16)) + ggtitle("CODEX Cell Type Frequencies Per Sample")

ggplot(data=cell_counts_coarse, aes(x=reorder(Var2, Freq),y=Freq, fill=Var1)) +
  geom_bar(position="fill",stat="identity") + 
  theme_bw() + 
  labs(x='Sample') + 
  labs(y='Cell Type') + scale_fill_manual(values=cols) + theme(axis.text = element_text(size=12)) + ggtitle("Negative Staging Marrow Predicted Cell Frequency") + RotatedAxis()

# export labeled cluster mask

NSM1_md <- subset(NSM.combined@meta.data, subset = orig.ident == 'SB67_NBM48_NSM1_1720_CODEX_Mesmer')
NSM1_md$cluster_anno_l2_n <- as.numeric(factor(NSM1_md$classified_cluster_anno_l2))
write.csv(x= data.frame(NSM1_md$x.coord, NSM1_md$y.coord, NSM1_md$CellID, NSM1_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/cluster_mask_csv/SB67_NBM48_NSM1_1720_mesmer_cluster_anno_l2_predictions.csv")

NSM1_md <- subset(NSM.combined@meta.data, subset = orig.ident == 'SB67_NBM48_NSM1_1720_CODEX_Mesmer')
NSM1_md$cluster_anno_coarse_n <- as.numeric(factor(NSM1_md$classified_cluster_anno_coarse))
write.csv(x= data.frame(NSM1_md$x.coord, NSM1_md$y.coord, NSM1_md$CellID, NSM1_md$cluster_anno_coarse_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/cluster_mask_csv/SB67_NBM48_NSM1_1720_mesmer_cluster_anno_coarse_predictions.csv")

NSM2_md <- subset(NSM.combined@meta.data, subset = orig.ident == 'SB67_NBM49_NSM2_1086_CODEX_Mesmer')
NSM2_md$cluster_anno_l2_n <- as.numeric(factor(NSM2_md$classified_cluster_anno_l2))
write.csv(x= data.frame(NSM2_md$x.coord, NSM2_md$y.coord, NSM2_md$CellID, NSM2_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/cluster_mask_csv/SB67_NBM49_NSM2_1086_mesmer_cluster_anno_l2_predictions.csv")

NSM2_md <- subset(NSM.combined@meta.data, subset = orig.ident == 'SB67_NBM49_NSM2_1086_CODEX_Mesmer')
NSM2_md$cluster_anno_coarse_n <- as.numeric(factor(NSM2_md$classified_cluster_anno_coarse))
write.csv(x= data.frame(NSM2_md$x.coord, NSM2_md$y.coord, NSM2_md$CellID, NSM2_md$cluster_anno_coarse_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/cluster_mask_csv/SB67_NBM49_NSM2_1086_mesmer_cluster_anno_coarse_predictions.csv")

NSM3_md <- subset(NSM.combined@meta.data, subset = orig.ident == 'SB67_NBM47_NSM3_1996_CODEX_Mesmer')
NSM3_md$cluster_anno_l2_n <- as.numeric(factor(NSM3_md$classified_cluster_anno_l2))
write.csv(x= data.frame(NSM3_md$x.coord, NSM3_md$y.coord, NSM3_md$CellID, NSM3_md$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/cluster_mask_csv/SB67_NBM47_NSM3_1996_mesmer_cluster_anno_l2_predictions.csv")

NSM3_md <- subset(NSM.combined@meta.data, subset = orig.ident == 'SB67_NBM47_NSM3_1996_CODEX_Mesmer')
NSM3_md$cluster_anno_coarse_n <- as.numeric(factor(NSM3_md$classified_cluster_anno_coarse))
write.csv(x= data.frame(NSM3_md$x.coord, NSM3_md$y.coord, NSM3_md$CellID, NSM3_md$cluster_anno_coarse_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/cluster_mask_csv/SB67_NBM47_NSM3_1996_mesmer_cluster_anno_coarse_predictions.csv")



H10_md_refined <- subset(immune.combined@meta.data, subset = orig.ident == 'SB67_NBM27_H10_CODEX_Mesmer')
H10_md_refined$cluster_anno_l2_n <- as.numeric(factor(H10_md_refined$cluster_anno_l2))
write.csv(x= data.frame(H10_md_refined$x.coord, H10_md_refined$y.coord, H10_md_refined$CellID, H10_md_refined$cluster_anno_l2_n), file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/Cluster_CSVs/SB67_H10_refined_mesmer_cluster_anno_l2_annotations_112322.csv")

# plot the score by cell type 
ggplot(NSM.combined@meta.data, aes(x=NSM.combined$predicted.MSC_refmap, y = NSM.combined$predicted.MSC_refmap.score, fill = predicted.MSC_refmap)) + geom_boxplot(width=0.7) + RotatedAxis() 

# plot the histogram of the scores
ggplot(NSM.combined@meta.data, aes(x = NSM.combined$predicted.MSC_refmap.score, fill = predicted.MSC_refmap)) + geom_histogram(bins=100) + RotatedAxis() 


# Compare the Frequency of cell types NSM staging marrow from the manually annotated NBM samples ---- 

NSM_ct_freqs <- NSM.combined@meta.data %>%
  group_by(Sample_Name, predicted_anno_coarse) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

NBM_ct_freqs <- immune.combined@meta.data %>%
  group_by(Sample_Name, cluster_anno_coarse) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

colnames(NSM_ct_freqs)[2] <- "cluster_anno_coarse"

median_NBM_ct_freqs <- NBM_ct_freqs %>% group_by(cluster_anno_coarse) %>% summarise(median = median(freq), sd(freq))
median_NBM_ct_freqs <- median_NBM_ct_freqs %>% dplyr::filter(median_NBM_ct_freqs$cluster_anno_coarse != "Schwann Cells")
median_NSM_ct_freqs <- NSM_ct_freqs %>% group_by(cluster_anno_coarse) %>% summarise(median = median(freq), sd(freq))

median_NSM_ct_freqs$difference <- median_NSM_ct_freqs$median - median_NBM_ct_freqs$median
median_NSM_ct_freqs$FC <- median_NSM_ct_freqs / median_NBM_ct_freqs

median_NSM_ct_freqs %>% ggplot(aes(x=cluster_anno_coarse, y =difference, fill = cluster_anno_coarse)) +
  geom_bar(stat='identity') + ylim(-0.5,0.5) + 
  theme_bw() + 
  labs(x='Sample') + 
  labs(y='Cell Type') + scale_fill_manual(values=cols) + 
  theme(axis.text = element_text(size=12)) + 
  ggtitle("Deviation from NBM Median for Negative Staging Marrows") + 
  RotatedAxis()


# Figure 7A illustrate ref mapping ----
p1 <- DimPlot(NSM.combined, raster = TRUE, raster.dpi = c(32,32), split.by = "orig.ident") & coord_fixed() & NoAxes() 
p1_raster <- rasterize(p1)
ggsave(p1_raster,width = 6, height = 6, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/Figures/PanelA_NSM_unlabeled_UMAP.pdf")

p1 <- DimPlot(NSM.combined, raster = TRUE, raster.dpi = c(300,300), group.by = "classified_cluster_anno_l2", pt.size = , split.by = "orig.ident", cols = cal2_cols) & coord_fixed() & NoAxes() & NoLegend()
p1_raster <- rasterize(p1)
ggsave(p1_raster,width = 6, height = 6, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/Figures/PanelA_NSM_labeled_UMAP.pdf")

p1 <- DimPlot(AML1_183_mapped, raster = TRUE, raster.dpi = c(300,300), split.by = "orig.ident", pt.size = 3) & coord_fixed() & NoAxes() 
p1_raster <- rasterize(p1)
ggsave(p1_raster,width = 6, height = 6, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/Figures/PanelA_AML_183_unlabeled_UMAP.pdf")



p1 <- DimPlot(AML1_183_mapped, raster = TRUE, raster.dpi = c(300,300), group.by = "classified_cluster_anno_l2", pt.size = 3, split.by = "orig.ident", cols = cal2_cols) & coord_fixed() & NoAxes() & NoLegend()
p1_raster <- rasterize(p1)
ggsave(p1_raster,width = 6, height = 6, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/Figures/PanelA_AML_183_labeled_UMAP.pdf")

# Supplemental Figure S10A CODEX NPM1C of NSM Samples ----
VlnPlot(NSM.combined, features = "codex_NPM1C", slot = "data", pt.size = 0, group.by = "classified_cluster_anno_l2", cols = cal2_cols, sort = TRUE) + NoLegend()


