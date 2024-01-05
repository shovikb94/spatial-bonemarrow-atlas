library(Seurat)
library(ggplot2)

# load combined scRNA-Seq object and prepare CellChat Obj
scrna_cal2_cols <- c("#FFD580", "#FFBF00", "#B0DFE6", "#7DB954", "#64864A", "#8FCACA", "#4682B4", "#CAA7DD", "#B6D0E2", "#a15891", "#FED7C3", "#A8A2D2", "#CF9FFF", "#9C58A1", "#2874A6", "#96C5D7", "#63BA97", "#BF40BF", "#953553", "#6495ED", "#E7C7DC", "#5599C8", "#FA8072", "#F3B0C3", "#F89880", "#40B5AD", "#019477", "#97C1A9", "#C6DBDA", "#CCE2CB", "#79127F", "#ff9d5c", "#FFC8A2", "#DD3F4E")
cal2_col_names <- c("Adipo-MSC", "AEC", "Ba/Eo/Ma", "CD4+ T-Cell", "CD8+ T-Cell", "CLP", "Cycling DCs", "Cycling HSPC", "Early Myeloid Progenitor", "Erythroblast", "Fibro-MSC", "GMP", "HSC", "Late Erythroid", "Late Myeloid", "Macrophages", "Mature B", "Megakaryocyte", "MEP", "Monocyte", "MPP", "Neutrophil", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "pDC", "Plasma Cell", "Pre-B", "Pre-Pro B", "Pro-B", "RBC", "SEC", "THY1+ MSC", "VSMC")
names(scrna_cal2_cols) <- cal2_col_names
cell_order <- c("HSC", "MPP", "Cycling HSPC", "MEP", "Erythroblast", "Late Erythroid", "RBC", "Megakaryocyte", "GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil", "Monocyte", "Macrophages", "Ba/Eo/Ma", "Cycling DCs", "pDC", "CLP", "Pre-Pro B", "Pro-B", "Pre-B", "Mature B", "Plasma Cell", "CD4+ T-Cell", "CD8+ T-Cell", "AEC", "SEC", "VSMC", "THY1+ MSC", "Adipo-MSC", "OsteoFibro-MSC", "Osteo-MSC", "Osteoblast", "Fibro-MSC")


combined <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Seurat/SB66_combined_corrected.RDS")
combined <- subset(combined, subset = cluster_anno_l2 != "NKT Cell" & cluster_anno_l2 != "RNAlo MSC") # remove RNAlo MSCs and likely artifactual NKT cell population
combined$cluster_anno_l2 <- droplevels(combined$cluster_anno_l2)
combined$cluster_anno_l2 <- factor(combined$cluster_anno_l2, levels = c("HSC", "MPP", "Cycling HSPC", "MEP", "Erythroblast", "Late Erythroid", "RBC", "Megakaryocyte", "GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil", "Monocyte", "Macrophages", "Ba/Eo/Ma", "Cycling DCs", "pDC", "CLP", "Pre-Pro B", "Pro-B", "Pre-B", "Mature B", "Plasma Cell", "CD4+ T-Cell", "CD8+ T-Cell", "AEC", "SEC", "VSMC", "THY1+ MSC", "Adipo-MSC", "OsteoFibro-MSC", "Osteo-MSC", "Osteoblast", "Fibro-MSC"))
cell_order <- c("HSC", "MPP", "Cycling HSPC", "MEP", "Erythroblast", "Late Erythroid", "RBC", "Megakaryocyte", "GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil", "Monocyte", "Macrophages", "Ba/Eo/Ma", "Cycling DCs", "pDC", "CLP", "Pre-Pro B", "Pro-B", "Pre-B", "Mature B", "Plasma Cell", "CD4+ T-Cell", "CD8+ T-Cell", "AEC", "SEC", "VSMC", "THY1+ MSC", "Adipo-MSC", "OsteoFibro-MSC", "Osteo-MSC", "Osteoblast", "Fibro-MSC")


# violin plots 

p1 <- VlnPlot(combined, features = c("CD34", "CD38", "SPINK2", "THY1", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p1, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/HSC_stacked_box_violin.pdf")

p2 <- VlnPlot(combined, features = c("CD34", "CD38", "SPINK2", "IL3RA", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p2, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/GMP_stacked_box_violin.pdf")

p3 <- VlnPlot(combined, features = c("CD34", "CD33", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p3, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/GMP_Myeloblast_stacked_box_violin.pdf")

p4 <- VlnPlot(combined, features = c("MPO", "CD33", "ITGAM"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p4, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/EarlyMyeloid_stacked_box_violin.pdf")

p5 <- VlnPlot(combined, features = c("MPO", "FUT4", "ITGAM"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p5, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/IntermediateMyeloid_stacked_box_violin.pdf")

p6 <- VlnPlot(combined, features = c("MPO", "CD33", "THBD", "ITGAM"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p6, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/MatureMyeloid_stacked_box_violin.pdf")

p7 <- VlnPlot(combined, features = c("CXCL12", "FOXC1", "THY1", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p7, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/AdipoMSC_stacked_box_violin.pdf")

p8 <- VlnPlot(combined, features = c("GATA1", "TFRC", "GYPC", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p8, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/Erythroblast_stacked_box_violin.pdf")

p9 <- VlnPlot(combined, features = c("GATA1", "ITGB3", "TGFB1", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p9, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/GATA1pos_Mk_stacked_box_violin.pdf")

p10 <- VlnPlot(combined, features = c("GATA1", "ITGB3", "TGFB1", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p10, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/GATA1pos_Mk_stacked_box_violin.pdf")

p11 <- VlnPlot(combined, features = c("GATA1", "CD34", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p11, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/MEP_EarlyEB_Mk_stacked_box_violin.pdf")

p12 <- VlnPlot(combined, features = c("PLP1", "NGFR"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p12, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/SchwannCell_stacked_box_violin.pdf")

p13 <- VlnPlot(combined, features = c("ACTA2", "CDH5"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p13, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/VSMC_stacked_box_violin.pdf")

p14 <- VlnPlot(combined, features = c("CD14","HLA-DRA","TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p14, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/Monocytes_stacked_box_violin.pdf")

p15 <- VlnPlot(combined, features = c("CD14","HLA-DRA","ITGAX", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p15, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/NonClassicalMono_stacked_box_violin.pdf")

p16 <- VlnPlot(combined, features = c("CD68","CD163","VCAM1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p16, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/Macrophage_stacked_box_violin.pdf")

p17 <- VlnPlot(combined, features = c("CD34","IL3RA","TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p17, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/pDC_stacked_box_violin.pdf")

p18 <- VlnPlot(combined, features = c("CD3E","CD4"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p18, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/CD4T_stacked_box_violin.pdf")

p19 <- VlnPlot(combined, features = c("CD3E","CD8A"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p19, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/CD8T_stacked_box_violin.pdf")

p20 <- VlnPlot(combined, features = c("PAX5","CD79A", "CD38", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p20, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/BCell_stacked_box_violin.pdf")

p21 <- VlnPlot(combined, features = c("CD79A","CD38", "SDC1", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p21, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/PlasmaCell_stacked_box_violin.pdf")

p22 <- VlnPlot(combined, features = c("NCAM1","VIM"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p22, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/Endosteal_stacked_box_violin.pdf")

p23 <- VlnPlot(combined, features = c("CDH5","CXCL12", "ACTA2"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p23, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/AEC_stacked_box_violin.pdf")

p24 <- VlnPlot(combined, features = c("CDH5","CXCL12", "ACTA2"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p24, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/SEC_stacked_box_violin.pdf")

p25 <- VlnPlot(combined, features = c("CD34","CD38", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p25, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/CLP_stacked_box_violin.pdf")

p26 <- VlnPlot(combined, features = c("CD34","SPINK2","CDH5", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p26, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/HSPC_stacked_box_violin.pdf")

p27 <- VlnPlot(combined, features = c("CD34","SPINK2", "CDH5", "TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p27, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/SPINK2+_HSPC_stacked_box_violin.pdf")

p28 <- VlnPlot(combined, features = c("CD34","ITGB3","TPSAB1"), pt.size = 0, fill.by = "ident", stack= TRUE, flip = TRUE, cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p28, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/CD34_CD61pos_stacked_box_violin.pdf")

p29 <- VlnPlot(combined, features = c("MCAM"), pt.size = 0, fill.by = "ident",cols = scrna_cal2_cols, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 8))
ggsave(p29, units = "in", device = "pdf", width = 6, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/Adpiocyte_stacked_box_violin.pdf")




p1 <- VlnPlot(combined, features = c("CD34"), pt.size = 0, cols = scrna_cal2_cols, sort = TRUE, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 6))
ggsave(p1, units = "in", device = "pdf", width = 4, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/CD34_box_violin.pdf")

p2 <- VlnPlot(combined, features = c("CD38"), pt.size = 0, cols = scrna_cal2_cols, sort = TRUE, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 6))
ggsave(p2, units = "in", device = "pdf",width = 4, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/CD38_box_violin.pdf")

p3 <- VlnPlot(combined, features = c("THY1"), pt.size = 0, cols = scrna_cal2_cols, sort = TRUE, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 6))
ggsave(p3, units = "in", device = "pdf",width = 4, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/CD90_box_violin.pdf")

p4 <- VlnPlot(combined, features = c("SPINK2"), pt.size = 0, cols = scrna_cal2_cols, sort = TRUE, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 6))
ggsave(p4, units = "in", device = "pdf",width = 4, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/SPINK2_box_violin.pdf")

p5 <- VlnPlot(combined, features = c("IL3RA"), pt.size = 0, cols = scrna_cal2_cols, sort = TRUE, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 6))
ggsave(p5, units = "in", device = "pdf",width = 4, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/CD123_box_violin.pdf")

p6 <- VlnPlot(combined, features = c("TPSAB1"), pt.size = 0, cols = scrna_cal2_cols, sort = TRUE, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 6))
ggsave(p6, units = "in", device = "pdf",width = 4, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/MCT_box_violin.pdf")

p7 <- VlnPlot(combined, features = c("CD33"), pt.size = 0, cols = scrna_cal2_cols, sort = TRUE, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 6))
ggsave(p7, units = "in", device = "pdf",width = 4, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/CD33_box_violin.pdf")

p8 <- VlnPlot(combined, features = c("MPO"), pt.size = 0, cols = scrna_cal2_cols, sort = TRUE, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 6))
ggsave(p8, units = "in",device = "pdf", width = 4, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/MPO_box_violin.pdf")

p9 <-VlnPlot(combined, features = c("THBD"), pt.size = 0, cols = scrna_cal2_cols, sort = TRUE, group.by = "cluster_anno_l2") + geom_boxplot(outlier.shape = NA) + NoLegend() + theme(axis.text = element_text(size = 6))
ggsave(p9, units = "in",device = "pdf", width = 4, height=3, filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_CellTypeAnnoGuide/figures/THBD_box_violin.pdf")







