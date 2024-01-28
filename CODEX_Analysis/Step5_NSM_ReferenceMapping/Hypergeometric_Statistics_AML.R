#Neighborhood hypergeometric test 
#Hypergeometric enrichment test
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(pheatmap)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(reshape2)

setwd("/mnt/isilon/tan_lab/sussmanj/CODEX/Seurat_Analysis/Bone_Marrow/")
###################################
#For AML - Dx
neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only_combined.neighborhood_blasts_separated_DxOnly.csv")
fc <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only.combined_neighborhood_fc_blasts_separated_DxOnly.csv")
neighborhood_mat <- as.matrix(fc[2:31])
rownames(neighborhood_mat) <- rownames(fc)
neighborhood_order <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)

neighborhood_counts <- as.numeric(table(factor(neighborhoods$neighborhood10, levels = neighborhood_order)))

neighborhoods_names <- unique(neighborhoods$neighborhood10)
annotations <- unique(neighborhoods$classified_cluster_anno_l2)
p_value_matrix <- matrix(NA, nrow = length(neighborhoods_names), ncol = length(annotations),
                         dimnames = list(neighborhoods_names, annotations))
n = nrow(neighborhoods)
for (i in 1:length(neighborhoods_names)) {
  print(i)
  for (j in 1:length(annotations)) {
    print(j)
    # Get counts for each combination of categories
    obs_count <- sum(neighborhoods$neighborhood10 == neighborhoods_names[i] & 
                       neighborhoods$classified_cluster_anno_l2 == annotations[j])
    
    # Perform hypergeometric test
    p_value <- phyper(obs_count - 1, 
                      sum(neighborhoods$classified_cluster_anno_l2 == annotations[j]),
                      n-sum(neighborhoods$classified_cluster_anno_l2 == annotations[j]),
                      sum(neighborhoods$neighborhood10 == neighborhoods_names[i]),
                      lower.tail = FALSE)
    print(p_value)
    p_value_matrix[i, j] <- p_value
  }
}
write.table(p_value_matrix, file = "AML_Dx_HG_Pvalue.txt", sep = '\t', quote = F)
long_df <- melt(p_value_matrix, variable.name = "Variable", value.name = "Value")
long_df$Var1 = long_df$Var1+1
colnames(long_df) <- c("CN", "Cell", "pval")
long_df$fdr = p.adjust(long_df$pval, method = "BH")
write.table(long_df, file = "AML_Dx_HG_Pvalue_BH_Adjusted.txt", sep = '\t', quote = F)
long_df$pval = NULL
bh_matrix <- pivot_wider(long_df, names_from = Cell, values_from = fdr)
bh_matrix <- bh_matrix[order(bh_matrix$CN), ]
rownames(bh_matrix) = bh_matrix$CN
bh_matrix$CN = NULL
bh_matrix <- as.data.frame(bh_matrix)

ha = HeatmapAnnotation(cell_counts = anno_barplot(x = neighborhood_counts,
                                                  bar_width = 1, 
                                                  gp = gpar(col = "white", fill = "grey"), 
                                                  border = TRUE,
                                                  axis_param = list(at = c(0, 1e4),
                                                                    labels = c("0", "10k")),
                                                  height = unit(4, "cm"), width = unit(1,"cm"), show_annotation_name = TRUE), 
                       annotation_label = c("Counts"),annotation_name_side = 'top',
                       annotation_name_rot = 360, annotation_name_align = TRUE,
                       border = TRUE, which = 'row')
htmp <- Heatmap(neighborhood_mat, name = "mat", #cluster_rows = F, cluster_columns = F, 
                cell_fun = function(j, i, x, y, w, h, fill){
                  if(bh_matrix[i, j] < 0.05) {
                    grid.text("*", x, y)
                  } 
                },
                rect_gp = gpar(col = "black", lwd = 2), column_names_rot = 45, 
                column_title = "Dx AML Neighborhood Enrichment", right_annotation =  ha, 
                column_names_gp = grid::gpar(fontsize = 14), row_title_gp = grid::gpar(fontsize = 13))
draw(htmp, heatmap_legend_side="left")

#For AML - PostTx
neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only_combined.neighborhood_blasts_separated_PostTx.csv")
fc <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only.combined_neighborhood_fc_blasts_separated_PostTx.csv")
neighborhood_mat <- as.matrix(fc[2:31])
rownames(neighborhood_mat) <- rownames(fc)
neighborhood_order <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)

neighborhood_counts <- as.numeric(table(factor(neighborhoods$neighborhood10, levels = neighborhood_order)))

neighborhoods_names <- unique(neighborhoods$neighborhood10)
annotations <- unique(neighborhoods$classified_cluster_anno_l2)
p_value_matrix <- matrix(NA, nrow = length(neighborhoods_names), ncol = length(annotations),
                         dimnames = list(neighborhoods_names, annotations))
n = nrow(neighborhoods)
for (i in 1:length(neighborhoods_names)) {
  print(i)
  for (j in 1:length(annotations)) {
    print(j)
    # Get counts for each combination of categories
    obs_count <- sum(neighborhoods$neighborhood10 == neighborhoods_names[i] & 
                       neighborhoods$classified_cluster_anno_l2 == annotations[j])
    
    # Perform hypergeometric test
    p_value <- phyper(obs_count - 1, 
                      sum(neighborhoods$classified_cluster_anno_l2 == annotations[j]),
                      n-sum(neighborhoods$classified_cluster_anno_l2 == annotations[j]),
                      sum(neighborhoods$neighborhood10 == neighborhoods_names[i]),
                      lower.tail = FALSE)
    print(p_value)
    p_value_matrix[i, j] <- p_value
  }
}
write.table(p_value_matrix, file = "AML_PostTx_HG_Pvalue.txt", sep = '\t', quote = F)
long_df <- melt(p_value_matrix, variable.name = "Variable", value.name = "Value")
long_df$Var1 = long_df$Var1+1
colnames(long_df) <- c("CN", "Cell", "pval")
long_df$fdr = p.adjust(long_df$pval, method = "BH")
write.table(long_df, file = "AML_PostTx_HG_Pvalue_BH_Adjusted.txt", sep = '\t', quote = F)
long_df$pval = NULL
bh_matrix <- pivot_wider(long_df, names_from = Cell, values_from = fdr)
bh_matrix <- bh_matrix[order(bh_matrix$CN), ]
rownames(bh_matrix) = bh_matrix$CN
bh_matrix$CN = NULL
bh_matrix <- as.data.frame(bh_matrix)

ha = HeatmapAnnotation(cell_counts = anno_barplot(x = neighborhood_counts,
                                                  bar_width = 1, 
                                                  gp = gpar(col = "white", fill = "grey"), 
                                                  border = TRUE,
                                                  axis_param = list(at = c(0, 1e4),
                                                                    labels = c("0", "10k")),
                                                  height = unit(4, "cm"), width = unit(1,"cm"), show_annotation_name = TRUE), 
                       annotation_label = c("Counts"),annotation_name_side = 'top',
                       annotation_name_rot = 360, annotation_name_align = TRUE,
                       border = TRUE, which = 'row')
htmp <- Heatmap(neighborhood_mat, name = "mat", #cluster_rows = F, cluster_columns = F, 
                cell_fun = function(j, i, x, y, w, h, fill){
                  if(bh_matrix[i, j] < 0.05) {
                    grid.text("*", x, y)
                  } 
                },
                rect_gp = gpar(col = "black", lwd = 2), column_names_rot = 45, 
                column_title = "Post-Tx AML Neighborhood Enrichment", right_annotation =  ha, 
                column_names_gp = grid::gpar(fontsize = 14), row_title_gp = grid::gpar(fontsize = 13))
draw(htmp, heatmap_legend_side="left")
