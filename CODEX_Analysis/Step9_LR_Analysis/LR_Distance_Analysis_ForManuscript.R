# Script for performing distance calculation and figure making in Figure S8C

library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggrastr)
library(imcRtools)
library(SingleCellExperiment)
library(ggpubr)
set.seed(69420)
setwd("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand")

LR.combined_citrate_mapped = readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/LR.combined.mapped_unfiltered_citrate2.RDS")

#Compute distances - run with different cell types and save results for different cell types of interest 
NBM_LR_Citrate_SCE  <- SingleCellExperiment(LR.combined_citrate_mapped@assays$CODEX@data, colData=LR.combined_citrate_mapped@meta.data)
NBM_LR_Citrate_SCE <- minDistToCells(NBM_LR_Citrate_SCE, 
                                     x_cells = NBM_LR_Citrate_SCE$classified_cluster_anno_l2 == "AEC",
                                     coords = c("x.coord","y.coord"),
                                     img_id = "orig.ident", return_neg=F)
LR.combined_citrate_mapped$dist_to_MSC <- colData(NBM_LR_Citrate_SCE)$distToCells
saveRDS(LR.combined_citrate_mapped$dist_to_MSC, file = "dist_to_AEC_citrate.rds")

NBM_LR_Citrate_SCE <- minDistToCells(NBM_LR_Citrate_SCE, 
                                     x_cells = NBM_LR_Citrate_SCE$classified_cluster_anno_l2 == "Erythroid",
                                     coords = c("x.coord","y.coord"),
                                     img_id = "orig.ident", return_neg=F)
LR.combined_citrate_mapped$dist_to_MSC <- colData(NBM_LR_Citrate_SCE)$distToCells
saveRDS(LR.combined_citrate_mapped$dist_to_MSC, file = "dist_to_Erythroid_citrate.rds")

NBM_LR_Citrate_SCE <- minDistToCells(NBM_LR_Citrate_SCE, 
                                     x_cells = NBM_LR_Citrate_SCE$classified_cluster_anno_l2 == "THY1+ MSC",
                                     coords = c("x.coord","y.coord"),
                                     img_id = "orig.ident", return_neg=F)
LR.combined_citrate_mapped$dist_to_MSC <- colData(NBM_LR_Citrate_SCE)$distToCells
saveRDS(LR.combined_citrate_mapped$dist_to_MSC, file = "dist_to_THY1MSC_citrate.rds")

#Distance analysis with LR expression based on percentile citrate, based on expression level of target --> distance to source cell type 
##############################################################################################
###Plot1
##############################################################################################
target_cell <- subset(LR.combined_citrate_mapped, subset = classified_cluster_anno_l2_refined == "Monocytes") 
expression <- target_cell@assays$CODEX@data
expression_cell <- as.vector(expression["CD44",])
top_percentile <- quantile(expression_cell, probs = 0.75)
bottom_percentile <- quantile(expression_cell, probs = 0.25)

LR.combined_citrate_mapped$Level <- "Intermediate"
Pos_ids <- WhichCells(LR.combined_citrate_mapped, expression = codex_CD44 > top_percentile, slot = "data")
LR.combined_citrate_mapped$Level[Pos_ids] <- "Top Quartile"
Neg_ids <- WhichCells(LR.combined_citrate_mapped, expression = codex_CD44 < bottom_percentile, slot = "data")
LR.combined_citrate_mapped$Level[Neg_ids] <- "Bottom Quartile"

#Load in distance data
LR.combined_citrate_mapped$dist_to_THY1MSC <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/dist_to_THY1MSC_citrate.rds")

target_cell_distance <- subset(LR.combined_citrate_mapped, subset = classified_cluster_anno_l2_refined == "Monocytes")
target_cell_distance$Level = factor(target_cell_distance$Level, levels = c("Bottom Quartile", "Intermediate", "Top Quartile"))
table(target_cell_distance$Level)
saveRDS(target_cell_distance, file = "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Plots/Monocyte_CD44_to_THY1MSC.rds")

#Plot 1 Monocyte to THY1+ MSC ----
target_cell_distance <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Plots/Monocyte_CD44_to_THY1MSC.rds")

# Calculate the 99th percentile of dist_to_THY1MSC for each Level
percentile_99 <- target_cell_distance@meta.data %>%
  group_by(Level) %>%
  summarise(dist_to_THY1MSC_99 = quantile(dist_to_THY1MSC*0.5069, 0.99))

# Filter out data above the 99th percentile for each Level
filtered_data <- target_cell_distance@meta.data %>%
  left_join(percentile_99, by = "Level") %>%
  filter(dist_to_THY1MSC*0.5069 <= dist_to_THY1MSC_99)

# Plotting the filtered data
p1 <- ggplot(filtered_data, aes(x = Level, y = dist_to_THY1MSC*0.5069, fill = Level)) +
  geom_violin(scale = "width", color = NA) + 
  geom_boxplot(width=0.2, fill = "NA", outlier.shape = NA) +
  theme_bw() + 
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means(comparisons = list(c("Bottom Quartile", "Top Quartile")),
                     method.args = list(alternative = "greater")) +
  ggtitle("Monocyte to THY1+ MSC") +
  labs(x = "CD44 Level", y = "Min Distance to THY1+ MSC (um)") + scale_fill_manual(values = c("steelblue", "#F3D4FE", "#DC143C")) + coord_cartesian(ylim = c(0,max(filtered_data$dist_to_THY1MSC_99)+35))       

ggsave(p1, filename = "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Monocyte_to_THY1MSC_Final.pdf", device = "pdf", width = 3, height = 3)


##############################################################################################
#Plot 2 
##############################################################################################
target_cell <- subset(LR.combined_citrate_mapped, subset = classified_cluster_anno_l2_refined == "VSMC") 
expression <- target_cell@assays$CODEX@data
expression_cell <- as.vector(expression["NOTCH3",])
top_percentile <- quantile(expression_cell, probs = 0.75)
bottom_percentile <- quantile(expression_cell, probs = 0.25)

LR.combined_citrate_mapped$Level <- "Intermediate"
Pos_ids <- WhichCells(LR.combined_citrate_mapped, expression = codex_NOTCH3 > top_percentile, slot = "data")
LR.combined_citrate_mapped$Level[Pos_ids] <- "Top Quartile"
Neg_ids <- WhichCells(LR.combined_citrate_mapped, expression = codex_NOTCH3 < bottom_percentile, slot = "data")
LR.combined_citrate_mapped$Level[Neg_ids] <- "Bottom Quartile"

#Load in distance data
LR.combined_citrate_mapped$dist_to_AEC <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/dist_to_AEC_citrate.rds")

target_cell_distance <- subset(LR.combined_citrate_mapped, subset = classified_cluster_anno_l2_refined == "VSMC")
target_cell_distance$Level = factor(target_cell_distance$Level, levels = c("Bottom Quartile", "Intermediate", "Top Quartile"))
table(target_cell_distance$Level)
saveRDS(target_cell_distance, file = "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Plots/VSMC_NOTCH3_to_AEC.rds")

#Plot 2 VSMC to AEC ----
target_cell_distance <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Plots/VSMC_NOTCH3_to_AEC.rds")

# Calculate the 99th percentile of dist_to_AEC for each Level
percentile_99 <- target_cell_distance@meta.data %>%
  group_by(Level) %>%
  summarise(dist_to_AEC_99 = quantile(dist_to_AEC*0.5069, 0.99))

# Filter out data above the 99th percentile for each Level
filtered_data <- target_cell_distance@meta.data %>%
  left_join(percentile_99, by = "Level") %>%
  filter(dist_to_AEC*0.5069 <= dist_to_AEC_99)

# Plotting the filtered data
p2 <- ggplot(filtered_data, aes(x = Level, y = dist_to_AEC*0.5069, fill = Level)) +
  geom_violin(scale = "width", color = NA) + 
  geom_boxplot(width=0.2, fill = "NA", outlier.shape = NA) +
  theme_bw() + 
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means(comparisons = list(c("Bottom Quartile", "Top Quartile")),
                     method.args = list(alternative = "greater")) +
  ggtitle("VSMC to AEC") +
  labs(x = "NOTCH3 Level", y = "Min Distance to AEC (um)") + scale_fill_manual(values = c("steelblue", "#F3D4FE", "#DC143C")) + coord_cartesian(ylim = c(0,max(filtered_data$dist_to_AEC_99+30)))       

ggsave(p2, filename = "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/VSMC_to_AEC_Final.pdf", device = "pdf", width = 3, height = 3)

##############################################################################################
#Plot 3
##############################################################################################
target_cell <- subset(LR.combined_citrate_mapped, subset = classified_cluster_anno_l2_refined == "SEC") 
expression <- target_cell@assays$CODEX@data
expression_cell <- as.vector(expression["TIE2",])
top_percentile <- quantile(expression_cell, probs = 0.75)
bottom_percentile <- quantile(expression_cell, probs = 0.25)

LR.combined_citrate_mapped$Level <- "Intermediate"
Pos_ids <- WhichCells(LR.combined_citrate_mapped, expression = codex_TIE2 > top_percentile, slot = "data")
LR.combined_citrate_mapped$Level[Pos_ids] <- "Top Quartile"
Neg_ids <- WhichCells(LR.combined_citrate_mapped, expression = codex_TIE2 < bottom_percentile, slot = "data")
LR.combined_citrate_mapped$Level[Neg_ids] <- "Bottom Quartile"

#Load in distance data
LR.combined_citrate_mapped$dist_to_Erythroid <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/dist_to_Erythroid_citrate.rds")

target_cell_distance <- subset(LR.combined_citrate_mapped, subset = classified_cluster_anno_l2_refined == "SEC")
target_cell_distance$Level = factor(target_cell_distance$Level, levels = c("Bottom Quartile", "Intermediate", "Top Quartile"))
table(target_cell_distance$Level)
saveRDS(target_cell_distance, file = "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Plots/SEC_TIE_to_Erythroid.rds")

#Plot 3 SEC to Erythroid ----
target_cell_distance <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Plots/SEC_TIE_to_Erythroid.rds")

# Calculate the 99th percentile of dist_to_Erythroid for each Level
percentile_99 <- target_cell_distance@meta.data %>%
  group_by(Level) %>%
  summarise(dist_to_Erythroid_99 = quantile(dist_to_Erythroid*0.5069, 0.99))

# Filter out data above the 99th percentile for each Level
filtered_data <- target_cell_distance@meta.data %>%
  left_join(percentile_99, by = "Level") %>%
  filter(dist_to_Erythroid*0.5069 <= dist_to_Erythroid_99)

# Plotting the filtered data
p3 <- ggplot(filtered_data, aes(x = Level, y = dist_to_Erythroid*0.5069, fill = Level)) +
  geom_violin(scale = "width", color = NA) + 
  geom_boxplot(width=0.2, fill = "NA", outlier.shape = NA) +
  theme_bw() + 
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means(comparisons = list(c("Bottom Quartile", "Top Quartile")),
                     method.args = list(alternative = "greater")) +
  ggtitle("SEC to Erythroid") +
  labs(x = "TIE2 Level", y = "Min Distance to Erythroid (um)") + scale_fill_manual(values = c("steelblue", "#F3D4FE", "#DC143C")) + coord_cartesian(ylim = c(0,max(filtered_data$dist_to_Erythroid_99+10)))       

ggsave(p3, filename = "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/SEC_to_Erythroid_Final.pdf", device = "pdf", width = 3, height = 3)


##############################################################################################
#Plot 4
##############################################################################################
#Load in distance data
LR.combined_citrate_mapped$dist_to_THY1MSC <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/dist_to_THY1MSC_citrate.rds")

#Annotate CLP using manual annotation of a subset by hand, then read in the matching CellIDs that are found int the Seurat object
NBM67 <- subset(LR.combined_citrate_mapped, subset = orig.ident %in% c("NBM67_H14_Citrate_CODEX", "NBM67_H38_Citrate_CODEX", "NBM67_H341_Citrate_CODEX"))
NBM70 <- subset(LR.combined_citrate_mapped, subset = orig.ident %in% c("NBM70_H10_Citrate_CODEX", "NBM70_H26_Citrate_CODEX", "NBM70_H32_Citrate_CODEX"))
NBM67file = read.table(file = "NBM67_CLP_CellID.txt", sep = "\t", header = T)
NBM70file = read.table(file = "NBM70_CLP_CellID.txt", sep = "\t", header = T)
NBM67IDs = NBM67file$Pixel_Value
NBM70IDs = NBM70file$Pixel_Value

length(unique(NBM67IDs))-1
length(unique(NBM70IDs))-1
NBM67_CLR <- subset(NBM67, subset = CellID %in% NBM67IDs)
NBM70_CLR <- subset(NBM70, subset = CellID %in% NBM70IDs)
dim(NBM67_CLR)
dim(NBM70_CLR)
CLR_only = merge(NBM67_CLR, NBM70_CLR)
CLR_only$classified_cluster_anno_l2 = "CLP_Manual"
dim(CLR_only)

target_cell <- subset(CLR_only, subset = classified_cluster_anno_l2 == "CLP_Manual") 
expression <- target_cell@assays$CODEX@data
expression_cell <- as.vector(expression["IL7RA",])
top_percentile <- quantile(expression_cell, probs = 0.75)
bottom_percentile <- quantile(expression_cell, probs = 0.25)

CLR_only$Level <- "Intermediate"
Pos_ids <- WhichCells(CLR_only, expression = codex_IL7RA > top_percentile, slot = "data")
CLR_only$Level[Pos_ids] <- "Top Quartile"
Neg_ids <- WhichCells(CLR_only, expression = codex_IL7RA < bottom_percentile, slot = "data")
CLR_only$Level[Neg_ids] <- "Bottom Quartile"

target_cell_distance <- subset(CLR_only, subset = classified_cluster_anno_l2 == "CLP_Manual")
target_cell_distance$Level = factor(target_cell_distance$Level, levels = c("Bottom Quartile", "Intermediate", "Top Quartile"))
table(target_cell_distance$Level)
saveRDS(target_cell_distance, file = "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Plots/CLP_IL7RA_to_THY1MSC.rds")

#Plot 1 
target_cell_distance <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Plots/CLP_IL7RA_to_THY1MSC.rds")

# Calculate the 99th percentile of dist_to_THY1MSC for each Level
percentile_99 <- target_cell_distance@meta.data %>%
  group_by(Level) %>%
  summarise(dist_to_THY1MSC_99 = quantile(dist_to_THY1MSC*0.5069, 0.99))

# Filter out data above the 99th percentile for each Level
filtered_data <- target_cell_distance@meta.data %>%
  left_join(percentile_99, by = "Level") %>%
  filter(dist_to_THY1MSC*0.5069 <= dist_to_THY1MSC_99)

# Plotting the filtered data
p4 <- ggplot(filtered_data, aes(x = Level, y = dist_to_THY1MSC*0.5069, fill = Level)) +
  geom_violin(scale = "width", color = NA) + 
  geom_boxplot(width=0.2, fill = "NA", outlier.shape = NA) +
  theme_bw() + 
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means(comparisons = list(c("Bottom Quartile", "Top Quartile")),
                     method.args = list(alternative = "greater")) +
  ggtitle("SEC to THY1+ MSC") +
  labs(x = "IL7RA Level", y = "Min Distance to THY1+ MSC (um)") + scale_fill_manual(values = c("steelblue", "#F3D4FE", "#DC143C")) + coord_cartesian(ylim = c(0,max(filtered_data$dist_to_THY1MSC_99-20)))       

ggsave(p4, filename = "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/CLP_to_THY1MSC_Final.pdf", device = "pdf", width = 3, height = 3)

