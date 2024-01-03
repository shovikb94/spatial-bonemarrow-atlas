# This script is Step 4 of the NBM CODEX analysis pipeline and takes as input the immune_filtered_FINAL.RDS object from Step3
# from which we will perform neighborhood analysis


library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)

immune.filtered <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.filtered_FINAL.RDS")

# Add metadata columns and write annotated CSV for neighborhood analysis
immune.filtered@meta.data$x.coord_rev <- -(immune.filtered$x.coord)
immune.filtered@meta.data$y.coord_rev <- -(immune.filtered$y.coord)
immune.filtered@meta.data <- mutate(immune.filtered@meta.data, Region = "reg001")
write.csv(x = immune.filtered@meta.data, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/csv/NBM_CODEX_Atlas_annotated_ForNeighborhoods.csv")


### Go to Python script to generate neighborhoods

# import neighborhood csvs
fc <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/output/neighborhood_fc.csv")
neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/output/neighborhood.csv")

# Figure 5B ComplexHeatmap for neighborhoods ----
neighborhood_mat <- as.matrix(fc[2:33])
rownames(neighborhood_mat) <- rownames(fc)

# Annotate Neighborhoods (keep python zero index in mind)
neighborhood_names_zeroindex <- c('0' = "HSC / Mature Myeloid",
                                  "1" = "Erythroid/Myeloid",
                                  "2" = "PC/Arteriolar",
                                  "3" = "Erythroid",
                                  "4" = "Arteriolar",
                                  "5" = "Erythroid",
                                  "6" = "Lymphoid",
                                  "7" = "Erythroid/Myeloid/Lymphoid",
                                  "8" = "Early Myeloid / Endosteal",
                                  "9" = "Myeloid/Lymphoid",
                                  "10" = "HSPC/Intermediate Myeloid",
                                  "11" = "Erythroid/Myeloid/Lymphoid",
                                  "12" = "Erythroid/Myeloid",
                                  "13" = "Early Myeloid / Arteriolar",
                                  "14" = "Peri-Arterolar Lymphoid")

neighborhood_names <- c('1' = "HSC / Mature Myeloid",
                        "2" = "Erythroid/Myeloid",
                        "3" = "PC/Arteriolar",
                        "4" = "Erythroid",
                        "5" = "Arteriolar",
                        "6" = "Erythroid",
                        "7" = "Lymphoid",
                        "8" = "Erythroid/Myeloid/Lymphoid",
                        "9" = "Early Myeloid / Endosteal",
                        "10" = "Myeloid/Lymphoid",
                        "11" = "HSPC/Intermediate Myeloid",
                        "12" = "Erythroid/Myeloid/Lymphoid",
                        "13" = "Erythroid/Myeloid",
                        "14" = "Early Myeloid / Arteriolar",
                        "15" = "Peri-Arterolar Lymphoid")

cols <- met.brewer(name="Signac", n=12)

rownames(neighborhood_mat) <- neighborhood_names

neighborhood_order <- c(5,15,3,14,10,9,8,12,7,11,2,1,4,13,6)
neighborhood_counts <- as.numeric(table(factor(neighborhoods$neighborhood10, levels = neighborhood_order-1)))

ha = HeatmapAnnotation(cell_counts = anno_barplot(x = neighborhood_counts,
                                                  bar_width = 1, 
                                                  gp = gpar(col = "white", fill = "grey"), 
                                                  border = TRUE,
                                                  axis_param = list(at = c(0, 2e4, 4e4, 6e4, 8e4,1e5),
                                                                    labels = c("0", "20k", "40k", "60k", "80k","100k")),
                                                  height = unit(4, "cm"), width = unit(1,"cm"), show_annotation_name = FALSE), annotation_label = c("Counts"),annotation_name_side = 'top', annotation_name_rot = 360, annotation_name_align = TRUE,
                       border = TRUE, which = 'row')

htmp <- Heatmap(neighborhood_mat, name = "mat", rect_gp = gpar(col = "black", lwd = 2), 
                column_title = "Bone Marrow Neighborhood Enrichment",right_annotation =  ha, column_names_gp = grid::gpar(fontsize = 14), row_title_gp = grid::gpar(fontsize = 13))
draw(htmp, heatmap_legend_side="left")


# save each neighborhood mask 
neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/output/neighborhood.csv")
neighborhoods %>% dplyr::filter(orig.ident == "SB67_NBM27_H10_CODEX_Mesmer") -> H10_nbs
write.csv(H10_nbs, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/sample_mask_input/H10_nbs.csv")

neighborhoods %>% dplyr::filter(orig.ident == "SB67_NBM32_H37_CODEX_Mesmer ") -> H37_nbs
write.csv(H37_nbs, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/sample_mask_input/H37_nbs.csv")

neighborhoods %>% dplyr::filter(orig.ident == "SB67_NBM33_H36_CODEX_Mesmer") -> H36_nbs
write.csv(H36_nbs, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/sample_mask_input/H36_nbs.csv")

neighborhoods %>% dplyr::filter(orig.ident == "SB67_NBM39_H41_CODEX_Mesmer") -> H41_nbs
write.csv(H41_nbs, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/sample_mask_input/H41_nbs.csv")

neighborhoods %>% dplyr::filter(orig.ident == "SB67_NBM34_H38_CODEX_Mesmer") -> H38_nbs
write.csv(H38_nbs, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/sample_mask_input/H38_nbs.csv")

neighborhoods %>% dplyr::filter(orig.ident == "SB67_NBM31_H32_CODEX_Mesmer") -> H32_nbs
write.csv(H32_nbs, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/sample_mask_input/H32_nbs.csv")


# Figure 5E HIF1a analysis by neighborhood for EMP ----
rownames(neighborhoods) <- neighborhoods$`Unnamed: 0`
immune.filtered <- AddMetaData(object = immune.filtered, metadata = neighborhoods$neighborhood10, col.name = "Neighborhood")
# Rename neighborhoods based on above manual annotation
immune.filtered <- SetIdent(immune.filtered, value = "Neighborhood")
new.cluster.ids <- neighborhood_names_zeroindex[levels(immune.filtered)]
immune.filtered<- RenameIdents(immune.filtered, new.cluster.ids)
immune.filtered@meta.data["Neighborhoods_Annotated"] <- immune.filtered@active.ident  

nbs_cols = c("0"= "#A8D37F", "1" ="#51C8EB", "2" ="#AF7AB5", "3" ="#FBED24", "4" ="#FF0000", "5"="#FBED24","6"="#89919C", "7"="#BD7CB5", "8"="#4B409A", "9"="#FFA500","10"="#A15A26", "11"="#BD7CB5","12"="#BD7CB5", "13"="#51C8EB", "14"= "#4B409A")	
immune.filtered <- SetIdent(immune.filtered, value = "Neighborhood")
nbs_cols = nbs_cols[levels(immune.filtered)]
names(nbs_cols) = neighborhood_names_zeroindex[levels(immune.filtered)]

nbs_names <- c("HSC / Mature Myeloid", "Erythroid/Myeloid", "PC/Arteriolar", "Erythroid", "Arteriolar", "Erythroid", "Lymphoid", "Erythroid/Myeloid/Lymphoid", "Early Myeloid / Endosteal", "Myeloid/Lymphoid", "HSPC/Intermediate Myeloid", "Erythroid/Myeloid/Lymphoid", "Erythroid/Myeloid", "Early Myeloid / Arteriolar", "Peri-Arterolar Lymphoid")
nbs_cols <- c("#A8D37F", "#51C8EB", "#FF00FF", "#FBED24", "#FF0000", "#FBED24", "#89919C", "#BD7CB5", "#4B409A", "#FFA500", "#A15A26", "#BD7CB5", "#51C8EB", "#4B409A", "#FF007F")
names(nbs_cols) <- nbs_names
# emp
emp_df <- as.data.frame(t(as.data.frame(GetAssayData(subset(immune.filtered, cluster_anno_l2 == "Early Myeloid Progenitor"),slot = 'data'))))
emp_df$Neighborhoods <- subset(immune.filtered, cluster_anno_l2 == "Early Myeloid Progenitor")$Neighborhoods_Annotated
emp_df$Neighborhoods <- as.factor(emp_df$Neighborhoods)
p1 <- emp_df %>% ggplot(aes(x=fct_reorder(Neighborhoods, desc(HIF1A)), y = HIF1A, fill = Neighborhoods)) + geom_boxplot(width=0.5, outlier.shape = NA) + scale_fill_manual(values = nbs_cols)  + theme_classic() + 
  scale_y_continuous(limits = quantile(emp_df$HIF1A, c(0.1, 0.95))) + xlab("Neighborhood") + ylab("HIF1a CLR Normalized") + ggtitle("Early Myeloid Progenitor HIF1a Expression") & RotatedAxis()

ggsave(p1, device = "pdf", height = 3, width = 5, filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/Figures/HIF1a_expression_short.pdf")

# t.test for each emp neighborhood against all cells used in Figure 5E
test <- c()
c <- c()
d <- c()
for (i in levels(emp_df$Neighborhoods)) {
  a <- subset(emp_df, Neighborhoods == {i})$HIF1A
  d <- append(x = d, values = i)
  b <- subset(emp_df, Neighborhoods != {i})$HIF1A
  test <- t.test(a,b, alternative = "less")
  c <- append(x = c, values = test$p.value)
}
#correct for multiple hypothesis 
c <- p.adjust(c,method = "BH")
names(c) <- d
c # contains adjusted p-values for all neighborhoods 

# Supplemental Figure 5D - CODEX HIF1a expression per cell type -----
cal2_cols <- c("#CF9FFF", "#E7C7DC", "#CAA7DD", "#A8A2D2", "#B6D0E2", "#2874A6", "#5599C8", "#AEC6CF", "#6495ED", "#64b8ed", "#96C5D7", "#40B5AD", "#8FCACA", "#CCE2CB", "#63BA97", "#7DB954", "#64864A", "#019477", "#953553", "#A1045A", "#a15891", "#9C58A1", "#79127F", "#BF40BF", "#FFD580", "#FFC8A2", "#FDDA0D", "#F3B0C3", "#FFBF00", "#ff9d5c", "#DD3F4E", "#FF69B4")
names(cal2_cols) <- c("HSC", "SPINK2+ HSPC", "HSPC", "GMP", "GMP/Myeloblast", "Early Myeloid Progenitor", "Intermediate Myeloid", "Mature Myeloid", "Monocytes", "Non-Classical Monocyte", "Macrophages", "pDC", "CLP", "Immature_B_Cell", "B-Cells", "CD4+ T-Cell", "CD8+ T-Cell", "Plasma Cells", "MEP/Early Erythroblast", "CD34+ CD61+", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "Adipo-MSC", "THY1+ MSC", "Adipocyte", "Endosteal", "AEC", "SEC", "VSMC", "Schwann Cells")
VlnPlot(immune.filtered, slot = 'data', cols = cal2_cols, group.by = "cluster_anno_l2",features = "HIF1A", pt.size = 0, sort=TRUE) + NoLegend()

# Supplemental Figure X - Per Sample Nb Composition ----

neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/output/neighborhood.csv")
neighborhood_mat <- as.matrix(fc[2:33])

neighborhood_names_zeroindex_unique <- c('0' = "Mature Myeloid",
                                  "1" = "Erythroid/Myeloid_1",
                                  "2" = "PC/Arteriolar",
                                  "3" = "Erythroid_1",
                                  "4" = "Arteriolar",
                                  "5" = "Erythroid_2",
                                  "6" = "Lymphoid",
                                  "7" = "Erythroid/Myeloid/Lymphoid_1",
                                  "8" = "Early Myeloid / Endosteal",
                                  "9" = "Myeloid/Lymphoid",
                                  "10" = "Intermediate Myeloid",
                                  "11" = "Erythroid/Myeloid/Lymphoid_2",
                                  "12" = "Erythroid/Myeloid_2",
                                  "13" = "Early Myeloid / Arteriolar",
                                  "14" = "Peri-Arterolar Lymphoid")

neighborhood_names_zeroindex <- c('0' = "Mature Myeloid",
                                  "1" = "Erythroid/Myeloid",
                                  "2" = "PC/Arteriolar",
                                  "3" = "Erythroid",
                                  "4" = "Arteriolar",
                                  "5" = "Erythroid",
                                  "6" = "Lymphoid",
                                  "7" = "Erythroid/Myeloid/Lymphoid",
                                  "8" = "Early Myeloid",
                                  "9" = "Myeloid/Lymphoid",
                                  "10" = "Intermediate Myeloid",
                                  "11" = "Erythroid/Myeloid/Lymphoid",
                                  "12" = "Erythroid/Myeloid",
                                  "13" = "Early Myeloid",
                                  "14" = "Peri-Arterolar Lymphoid")

neighborhoods$neighborhood_named <- neighborhood_names_zeroindex[as.character(neighborhoods$neighborhood10)]

nbs_names <- c("Mature Myeloid", "Erythroid/Myeloid", "PC/Arteriolar", "Erythroid", "Arteriolar", "Erythroid", "Lymphoid", "Erythroid/Myeloid/Lymphoid", "Early Myeloid", "Myeloid/Lymphoid", "Intermediate Myeloid", "Erythroid/Myeloid/Lymphoid", "Erythroid/Myeloid", "Early Myeloid", "Peri-Arterolar Lymphoid")
nbs_cols <- c("#A8D37F", "#51C8EB", "#FF00FF", "#FBED24", "#FF0000", "#FBED24", "#89919C", "#BD7CB5", "#4B409A", "#FFA500", "#A15A26", "#BD7CB5", "#51C8EB", "#4B409A", "#FF007F")
names(nbs_cols) <- nbs_names

# add age information and simplify sample names 
neighborhoods$Sample_Name <- gsub(".*_(H[0-9]+)_.*", "\\1", neighborhoods$orig.ident)
age <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Revisions/Age Analysis/additional_metadata_CODEX.csv")
colnames(age)[1] <- "Sample_Name"
neighborhoods <- neighborhoods %>% left_join(age, by = "Sample_Name")

p1 <- ggplot(neighborhoods, aes(fill=neighborhood_named, x=reorder(as.factor(Sample_Name), ~(Age)))) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") + scale_fill_manual(values=nbs_cols, limits=force) # limits=force drops unused levels from the legend

ggsave(p1, device = "pdf", width = 5, height = 5, filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/Figures/PerSample_NeighborhoodComposition.pdf")


