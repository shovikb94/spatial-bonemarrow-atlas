# Script for plotting heatmap of non-cell analysis

library(readr)
library(dplyr)
library(ComplexHeatmap)
library(tidyr)
library(reshape2)
library(circlize)
library(metap)
library(ggplot2)


# Figure 6A complex heatmap for structures ----
ranks <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Non-Cell Microenvironment Analysis/Most_Uptodate/combined_structure_v3.csv")
summarized_ranks <- ranks %>% group_by(structure_1, structure_2) %>% summarise(median_norm_rank = median(normalized_rank))
sum_ranks_st_mat <- dcast(summarized_ranks, structure_1 ~ structure_2, value.var = "median_norm_rank")
rownames(sum_ranks_st_mat) <- sum_ranks_st_mat$structure_1
sum_ranks_st_mat$structure_1 <- NULL
col_fun = colorRamp2(c(0, 0.5, 1), c("red", "white", "blue")) 

Heatmap(as.matrix(sum_ranks_st_mat), name = "mat", rect_gp = gpar(col = "black", lwd = 2), col = col_fun,
        column_title = "Cell Type Structural Proximity - Normalized Rank", clustering_method_rows = "single", row_names_gp = grid::gpar(fontsize = 8))


# note structure1 is "structure_to" and structure2 is "structure_from"

# Figure 6B - line length scaling by norm rank -----
to_bone_rank <- c(0.78986175, 0.24285714, 0.07142857,0.5857143,0.75714286,0.5267281)
names(to_bone_rank) <- colnames(sum_ranks_st_mat)
scalar1 <- function(x) {x / sqrt(sum(x^2))}
to_bone_rank_scaled <- scalar1(to_bone_rank) * 50
 
# Calculate using Stouffer's Method ----
stouffers <- function(pvec) {
  pvec[pvec == 0] <- pvec[pvec == 0] + 1E-16 # handle extreme values 
  pvec[pvec == 1] <- pvec[pvec == 1] - 1E-16 # handle extreme values
  zvec <- qnorm(pvec)
  #hist(zvec, breaks=50)
  sumz <- sum(zvec)
  zmeta <- sumz / sqrt(length(zvec))
  pmeta <- pnorm(zmeta)
  return(pmeta)
}

metap_sumz <- ranks %>% group_by(structure_1, structure_2) %>% summarise(pvec = stouffers(pvalue))
#metap_sumz <- ranks %>% group_by(celltype, structure) %>% summarise(meta_p = stouffers(p = ranks$pvalue))
metap_wide_sumz <- dcast(metap_sumz, structure_2 ~ structure_1, value.var = "pvec")
metap_wide_sumz_corplot <- metap_wide_sumz[,-1]
rownames(metap_wide_sumz_corplot) <- metap_wide_sumz$structure_2
metap_wide_sumz_corplot

# Figure 6C complex heatmap for neighborhoods ----

ranks <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Non-Cell Microenvironment Analysis/Most_Uptodate/combined_v3_health_neighbourhood.csv")
ranks$neighbourhood <- as.factor(ranks$neighbourhood)
ranks$neighbourhood <- gsub(ranks$neighbourhood, pattern = "HSPC/Intermediate Myeloid", replacement = "Intermediate Myeloid") # rename based on updated annotation 
ranks$neighbourhood <- gsub(ranks$neighbourhood, pattern = "HSC / Mature Myeloid", replacement = "Mature Myeloid") # rename based on updated annotation 
col_fun = colorRamp2(c(0, 0.5, 1), c("red", "white", "blue")) 

ranks$neighbourhood <- factor(ranks$neighbourhood)
summarized_ranks <- ranks %>% group_by(neighbourhood, structure) %>% summarise(median_norm_rank = median(normalized_rank))
sum_ranks_cn_mat <- dcast(summarized_ranks, neighbourhood ~ structure, value.var = "median_norm_rank")
rownames(sum_ranks_cn_mat) <- sum_ranks_cn_mat$neighbourhood
sum_ranks_cn_mat$neighbourhood <- NULL

Heatmap(as.matrix(sum_ranks_cn_mat), name = "mat", rect_gp = gpar(col = "black", lwd = 2), col = col_fun, 
        column_title = "Neighborhood Structural Proximity - Normalized Rank", clustering_method_rows = "single", row_names_gp = grid::gpar(fontsize = 8))

# Figure 6D complex heatmap for celltypes -----
ranks <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Non-Cell Microenvironment Analysis/Most_Uptodate/combined_v3_health_celltype.csv")
ranks$celltype <- as.factor(ranks$celltype)
# cell_order <- c("HSC", "SPINK2+ HSPC", "HSPC", "GMP", "GMP/Myeloblast", "Early Myeloid Progenitor", "Intermediate Myeloid", "Mature Myeloid", "Monocytes", "Non-Classical Monocyte", "Macrophages", "pDC", "CLP", "Immature_B_Cell", "B-Cells", "CD4+ T-Cell", "CD8+ T-Cell", "Plasma Cells", "MEP/Early Erythroblast", "CD34+ CD61+", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "MSC", "THY1+ MSC", "Adipocyte", "Endosteal", "AEC", "SEC", "VSMC", "Schwann Cells")
cell_order <- c("SPINK2+ HSPC", "HSPC","GMP","GMP/Myeloblast", "Early Myeloid Progenitor", "Intermediate Myeloid", "Mature Myeloid", "Monocytes", "Non-Classical Monocyte", "Macrophages", "pDC", "Immature_B_Cell", "B-Cells", "CD4+ T-Cell", "CD8+ T-Cell", "Plasma Cells", "MEP/Early Erythroblast", "CD34+ CD61+", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "MSC", "THY1+ MSC", "Adipocyte", "Endosteal", "AEC", "SEC", "VSMC")

ranks$celltype <- factor(ranks$celltype, levels = cell_order)
levels(ranks$celltype) <- cell_order
summarized_ranks <- ranks %>% group_by(celltype, structure) %>% summarise(median_norm_rank = median(normalized_rank))
sum_ranks_ct_mat <- dcast(summarized_ranks, celltype ~ structure, value.var = "median_norm_rank")
rownames(sum_ranks_ct_mat) <- sum_ranks_ct_mat$celltype
sum_ranks_ct_mat$celltype <- NULL

Heatmap(as.matrix(sum_ranks_ct_mat), name = "mat", rect_gp = gpar(col = "black", lwd = 2), col = col_fun, 
        column_title = "Cell Type Structural Proximity - Normalized Rank", clustering_method_rows = "single", row_names_gp = grid::gpar(fontsize = 8))

# Figure 6E Zoom in on HSPCs Structural Proximity ----
hspc_circ <- ranks %>% dplyr::filter(celltype == "SPINK2+ HSPC" | celltype == "HSPC")
# plot boxplot
hspc_circ %>% ggplot(aes(x = structure, y = normalized_rank, fill = structure, linetype = celltype))  + geom_boxplot() + scale_fill_manual(values = c("#FAA59E","#CADCEB","#C2E7B9","#E8DBEC","#FED194","#ffe2db")) + scale_linetype_manual(values=c("solid", "longdash")) + theme_minimal()
hspc_circ %>% ggplot(aes(x = structure, y = normalized_rank, fill = structure, linetype = celltype))  + geom_boxplot() + scale_fill_manual(values = c("#FAA59E","#CADCEB","#C2E7B9","#E8DBEC","#FED194","#ffe2db")) + scale_linetype_manual(values=c("solid", "longdash")) + theme_minimal() + coord_polar()

# Calculate using Stouffer's Method ----
stouffers <- function(pvec) {
pvec[pvec == 0] <- pvec[pvec == 0] + 1E-16 # handle extreme values 
pvec[pvec == 1] <- pvec[pvec == 1] - 1E-16 # handle extreme values
zvec <- qnorm(pvec)
#hist(zvec, breaks=50)
sumz <- sum(zvec)
zmeta <- sumz / sqrt(length(zvec))
pmeta <- pnorm(zmeta)
return(pmeta)
}

metap_sumz <- ranks %>% group_by(celltype, structure) %>% summarise(pvec = stouffers(pvalue))
#metap_sumz <- ranks %>% group_by(celltype, structure) %>% summarise(meta_p = stouffers(p = ranks$pvalue))
metap_wide_sumz <- dcast(metap_sumz, celltype ~ structure, value.var = "pvec")
metap_wide_sumz_corplot <- metap_wide_sumz[,-1]
rownames(metap_wide_sumz_corplot) <- metap_wide_sumz$celltype
metap_wide_sumz_corplot
## can use this to_plot table to check significane of any to_from comparison 

# Supplemental Figure 6B chord diagram ----
sum_ranks_ct_mat_trunc <- sum_ranks_ct_mat[c(1,2,4,5,7,10,13,20,23,26,27,28),]
summarized_ranks_sig <- summarized_ranks %>% left_join(metap,by = c("celltype" = "celltype", "structure" = "structure") )
summarized_ranks_sig_trunc <- summarized_ranks_sig %>% dplyr::filter(celltype %in% c("SPINK2+ HSPC", "HSPC", "GMP_Myeloblast", "Early Myeloid Progenitor", "Mature Myeloid", "Macrophages", "B-Cells", "Erythroid", "MSC", "AEC", "SEC") & meta_p <0.05)
summarized_ranks_sig_trunc$median_norm_rank <- -log(summarized_ranks_sig_trunc$median_norm_rank + 0.000000000001)
structure_cols <- c("#FAA59E","#CADCEB","#C2E7B9","#E8DBEC","#FED194","#ffe2db")
names(structure_cols) <- levels(as.factor(ranks$structure))
chord_to_plot <- data.frame(summarized_ranks_sig_trunc$structure, summarized_ranks_sig_trunc$celltype, summarized_ranks_sig_trunc$median_norm_rank)
chordDiagram(chord_to_plot, grid.col = c(cal2_cols,structure_cols), transparency = 0.2)

# Supplemental Figure S7G (new) -----
ranks_idx <- read_csv("~/Documents/Manuscripts/NBM_Atlas/ReviewerComments/AML_StructuralAnalysis_SuppFigS7G/combined_aml_idxneighbourhood.csv")
ranks_idx$neighbourhood <- ranks_idx$neighbourhood+1
ranks_idx$neighbourhood <- as.factor(ranks_idx$neighbourhood)

idx_circ <- ranks_idx %>% dplyr::filter(structure == "bone")
# plot boxplot
p1 <- idx_circ %>% ggplot(aes(x = reorder(neighbourhood,-normalized_rank), y = normalized_rank, fill = neighbourhood))  + geom_boxplot() + NoLegend() + theme_minimal() 


stouffers <- function(pvec) {
  pvec[pvec == 0] <- pvec[pvec == 0] + 1E-16 # handle extreme values 
  pvec[pvec == 1] <- pvec[pvec == 1] - 1E-16 # handle extreme values
  zvec <- qnorm(pvec)
  #hist(zvec, breaks=50)
  sumz <- sum(zvec)
  zmeta <- sumz / sqrt(length(zvec))
  pmeta <- pnorm(zmeta)
  return(pmeta)
}

metap_sumz <- ranks_idx %>% group_by(neighbourhood, structure) %>% summarise(pvec = stouffers(pvalue))
#metap_sumz <- ranks %>% group_by(celltype, structure) %>% summarise(meta_p = stouffers(p = ranks$pvalue))
metap_wide_sumz <- dcast(metap_sumz, neighbourhood ~ structure, value.var = "pvec")
metap_wide_sumz_corplot <- metap_wide_sumz[,-1]
rownames(metap_wide_sumz_corplot) <- metap_wide_sumz$neighbourhood


# repeat for post therapy
ranks_posttx <- read_csv("~/Documents/Manuscripts/NBM_Atlas/ReviewerComments/AML_StructuralAnalysis_SuppFigS7G/combined_aml_pt_neighbourhood.csv")
ranks_posttx$neighbourhood <- ranks_posttx$neighbourhood+1
ranks_posttx$neighbourhood <- as.factor(ranks_posttx$neighbourhood)

posttx_circ <- ranks_posttx %>% dplyr::filter(structure == "bone")
# plot boxplot
p2 <- posttx_circ %>% ggplot(aes(x = reorder(neighbourhood,-normalized_rank), y = normalized_rank, fill = neighbourhood))  + geom_boxplot() + NoLegend() + theme_minimal() 
hspc_circ %>% ggplot(aes(x = structure, y = normalized_rank, fill = structure, linetype = celltype))  + geom_boxplot() + scale_fill_manual(values = c("#FAA59E","#CADCEB","#C2E7B9","#E8DBEC","#FED194","#ffe2db")) + scale_linetype_manual(values=c("solid", "longdash")) + theme_minimal() + coord_polar()


stouffers <- function(pvec) {
  pvec[pvec == 0] <- pvec[pvec == 0] + 1E-16 # handle extreme values 
  pvec[pvec == 1] <- pvec[pvec == 1] - 1E-16 # handle extreme values
  zvec <- qnorm(pvec)
  #hist(zvec, breaks=50)
  sumz <- sum(zvec)
  zmeta <- sumz / sqrt(length(zvec))
  pmeta <- pnorm(zmeta)
  return(pmeta)
}

metap_sumz <- ranks_posttx %>% group_by(neighbourhood, structure) %>% summarise(pvec = stouffers(pvalue))
#metap_sumz <- ranks %>% group_by(celltype, structure) %>% summarise(meta_p = stouffers(p = ranks$pvalue))
metap_wide_sumz <- dcast(metap_sumz, neighbourhood ~ structure, value.var = "pvec")
metap_wide_sumz_corplot <- metap_wide_sumz[,-1]
rownames(metap_wide_sumz_corplot) <- metap_wide_sumz$neighbourhood

plot_grid(p1,p2, ncol = 1)

