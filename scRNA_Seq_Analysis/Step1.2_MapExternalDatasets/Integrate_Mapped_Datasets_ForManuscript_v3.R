# Integrate MSC Mapped Samples

library(Seurat)
library(dplyr)
library(viridis)

# save cluster_anno_l2 colors
cal2_cols <- c("#FFD580", "#FFBF00", "#B0DFE6", "#7DB954", "#64864A", "#8FCACA", "#4682B4", "#CAA7DD", "#B6D0E2", "#a15891", "#FED7C3", "#A8A2D2", "#CF9FFF", "#9C58A1", "#2874A6", "#96C5D7", "#63BA97", "#BF40BF", "#953553", "#6495ED", "#E7C7DC", "#5599C8", "#FA8072", "#F3B0C3", "#F89880", "#40B5AD", "#019477", "#97C1A9", "#C6DBDA", "#CCE2CB", "#79127F", "#FFC5BF", "#ff9d5c", "#FFC8A2", "#DD3F4E")
cal2_col_names <- c("Adipo-MSC", "AEC", "Ba/Eo/Ma", "CD4+ T-Cell", "CD8+ T-Cell", "CLP", "Cycling DCs", "Cycling HSPC", "Early Myeloid Progenitor", "Erythroblast", "Fibro-MSC", "GMP", "HSC", "Late Erythroid", "Late Myeloid", "Macrophages", "Mature B", "Megakaryocyte", "MEP", "Monocyte", "MPP", "Neutrophil", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "pDC", "Plasma Cell", "Pre-B", "Pre-Pro B", "Pro-B", "RBC", "RNAlo MSC", "SEC", "THY1+ MSC", "VSMC")
names(cal2_cols) <- cal2_col_names

# read in seurat obj for nbm atlas
combined <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Seurat/SB66_combined_corrected.RDS")

# filter out rna lo msc and likely doublet so-called NKT population
combined <- subset(combined, subset = cluster_anno_l2 != "RNAlo MSC" & cluster_anno_l2 != "NKT Cell")
combined <- RunUMAP(combined, dims = 1:30, reduction = "pca", reduction.name = "UMAP_dim30", return.model = TRUE) # rerun find umap to return the umap model

#MSCs <- subset(x = combined, subset = cluster_anno_l2 %in% c("Adipo-MSC", "THY1+ MSC", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "Fibro-MSC"))
#MSCs <- RunUMAP(combined, dims = 1:50, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim50", return.model = TRUE) # rerun find umap to return the umap model
MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Seurat/MSCs_Final_091523.RDS")
MSCs <- RunUMAP(MSCs, dims = 1:50, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim50", return.model = TRUE) # rerun find umap to return the umap model

# read in mapped datasets
DeJong_MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/DeJong_MSCs.RDS")
DeJong_MSCs$Dataset <- "DeJong_MSCs"
Jardine_Fetal_MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/fbm_MSCs.RDS")
Jardine_Fetal_MSCs$Dataset <- "Jardine_Fetal_MSCs"
Li_MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/Li_MSCs.RDS")
Li_MSCs$Dataset <- "Li_MSCs"
Wang_OA_MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/Wang_OA_MSCs.RDS")
Wang_OA_MSCs$Dataset <- "Wang_OA_MSCs"
Wang_OP_MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/Wang_op.RDS")
Wang_OP_MSCs$Dataset <- "Wang_OP_MSCs"
Ennis_MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/Ennis_MSCs.RDS")
Ennis_MSCs$Dataset <- "Ennis_MSCs"
Triana_MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/Triana_MSCs.RDS")
Triana_MSCs$Dataset <- "Triana_MSCs"

# merge all datasets
MSCs$predicted.MSC_refmap <- MSCs$cluster_anno_l2
MSCs$Dataset <- "Bandyopadhyay_MSCs"
Li_MSCs$predicted.MSC_refmap <- Li_MSCs$predicted.combined_refmap

Wang_OA_MSCs <- subset(Wang_OA_MSCs, subset = predicted.MSC_refmap %in% c("Adipo-MSC", "THY1+ MSC", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "Fibro-MSC"))
Wang_OP_MSCs <- subset(Wang_OP_MSCs, subset = predicted.MSC_refmap %in% c("Adipo-MSC", "THY1+ MSC", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "Fibro-MSC"))


Combined_merge <- merge(MSCs, y = c(DeJong_MSCs, Jardine_Fetal_MSCs, Li_MSCs, Wang_OA_MSCs, Wang_OP_MSCs, Ennis_MSCs, Triana_MSCs))

anchors <- FindTransferAnchors(
  reference = MSCs,
  query = Combined_merge,
  normalization.method = "LogNormalize",
  reference.reduction = "MSC_pca",
  dims = 1:50
)

Combined_merge <- MapQuery(
  anchorset = anchors,
  query = Combined_merge,
  reference = MSCs,
  refdata = list(
    MSC_refmap = "cluster_anno_l2"
  ),
  reference.reduction = "MSC_pca", 
  reduction.model = "MSC_UMAP_dim50"
)

ggplot(Combined_merge@meta.data, aes(fill=predicted.MSC_refmap, x=Dataset)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
 ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force) # limits=force drops unused levels from the legend

# Make Figure 2I -----

p1 <- ggplot(Combined_merge@meta.data, aes(fill=predicted.MSC_refmap, x=Dataset)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force) # limits=force drops unused levels from the legend

saveRDS(Combined_merge, "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/Combined_MappedMSCs.RDS")

DimPlot(Combined_merge, group.by = "Dataset", reduction = "ref.umap", split.by = "Dataset") & coord_fixed() & NoAxes()
DimPlot(subset(Combined_merge, subset = Dataset == "Bandyopadhyay_MSCs"), group.by = "cluster_anno_l2", label = TRUE, repel = TRUE)


ggsave(p1, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/figures/MSC_Composition_ByDataset.pdf")

p2 <- ggplot(Combined_merge@meta.data, aes(fill=predicted.MSC_refmap, x=Dataset)) + 
           geom_bar(stat="count") + theme_minimal() + RotatedAxis() +
           ylab("Total MSCs") + xlab("Sample") + scale_fill_manual(values=cal2_cols, limits=force) # limits=force drops unused levels from the legend
         
ggsave(p2, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/figures/MSC_Counts_ByDataset.pdf")

# Calculate differences between subsets across samples -----
rel_freq <- data.frame(table(Combined_merge$predicted.MSC_refmap, Combined_merge$Dataset))

rel_freq <- rel_freq %>% group_by(Var2) %>% mutate(TotalFreq = sum(Freq), RelFreq = Freq / TotalFreq) %>%
  select(-TotalFreq) # manually calculated Bandyopadhyay_MSCs vs. Jardine for text fold change numbers

# Calculate fold difference between our atlas and other aspirate-derived atlases -----
calculate_fold_difference <- function(data, main_dataset, compare_datasets) {
  # Filter out 'Adipo-MSC' and calculate total frequency for main_dataset
  main_freq <- data %>%
    filter(Var2 == main_dataset, Var1 != 'Adipo-MSC') %>%
    summarise(TotalFreq = sum(Freq)) %>%
    pull(TotalFreq)
  
  # Filter for other datasets, sum their frequencies, and calculate fold difference
  other_datasets_freq <- data %>%
    filter(Var2 %in% compare_datasets, Var1 != 'Adipo-MSC') %>%
    summarise(TotalFreq = sum(Freq)) %>%
    pull(TotalFreq)
  
  fold_difference <- main_freq / sum(other_datasets_freq)
  
  return(fold_difference)
}
calculate_fold_difference(rel_freq, "Bandyopadhyay_MSCs", c("DeJong_MSCs", "Ennis_MSCs", "Li_MSCs", "Triana_MSCs", "Wang_OA_MSCs", "Wang_OP_MSCs"))
