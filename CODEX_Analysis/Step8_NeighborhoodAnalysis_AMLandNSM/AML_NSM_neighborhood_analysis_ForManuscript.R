library(ComplexHeatmap)
library(readr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)

# read in objects
AML1_183_mapped <- readRDS("objects/AML1_183_mapped.RDS")
AML1_382_mapped <- readRDS("objects/AML1_382_mapped.RDS")
AML2_191_mapped <- readRDS("objects/AML2_191_mapped.RDS")
AML3_1329_mapped <- readRDS("objects/AML3_1329_mapped.RDS")
AML3_1443_mapped <- readRDS("objects/AML3_1443_mapped.RDS")


# Integrate AML and NSM samples
AML1_183_mapped@meta.data$Sample_Group <- "AML"
AML1_382_mapped@meta.data$Sample_Group<- "AML"
AML2_191_mapped@meta.data$Sample_Group <- "AML"
AML3_1329_mapped@meta.data$Sample_Group <- "AML"
AML3_1443_mapped@meta.data$Sample_Group <- "AML"
NSM.combined@meta.data$Sample_Group <- "NSM"

AML1_183_mapped@meta.data$Sample_Timepoint <- "Dx"
AML1_382_mapped@meta.data$Sample_Timepoint<- "Post-Tx"
AML2_191_mapped@meta.data$Sample_Timepoint <- "Dx"
AML3_1329_mapped@meta.data$Sample_Timepoint <- "Dx"
AML3_1443_mapped@meta.data$Sample_Timepoint <- "Post-Tx"
NSM.combined@meta.data$Sample_Timepoint <- "NSM"

All.combined <- merge(x=NSM.combined, , y = c(AML1_183_mapped, AML1_382_mapped, AML2_191_mapped, AML3_1329_mapped, AML3_1443_mapped), add.cell.ids = c("NSM.combined", "AML1_183", "AML1_382", "AML2_191", "AML3_1329", "AML_1443"))

# change cluster anno l2 to blasts
Blast_ids <- colnames(subset(All.combined, subset = Adjusted_Class == "NPM1_Mutant")) # Blast Rename
All.combined$classified_cluster_anno_l2[Blast_ids] <- "NPM1 Mutant Blast"

# Test for markers enriched pre and post therapy, GATA1 text reference ----
All.combined <- SetIdent(All.combined, value = "Adjusted_Class")
FindMarkers(All.combined, subset.ident = "NPM1_Mutant", group.by = "Sample_Timepoint", ident.2 = "Dx", ident.1 = "Post-Tx")
# GATA1 increase - p.adj = 0.000123, log2FC = 0.378

# remove artifacts
All.combined <- subset(All.combined, subset = classified_cluster_anno_l2 != "Artifact" & classified_cluster_anno_l2 != "CD44+ Undetermined" & classified_cluster_anno_l2 != "Undetermined" & classified_cluster_anno_l2 != "Autofluorescent")

# export for nb analysis
# Add metadata columns and write annotated CSV for neighborhood analysis
All.combined@meta.data$x.coord_rev <- -(All.combined$x.coord)
All.combined@meta.data$y.coord_rev <- -(All.combined$y.coord)
All.combined@meta.data <- mutate(All.combined@meta.data, Region = "reg001")
write.csv(x = All.combined@meta.data, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/For_Neighborhoods/NSM_AML_combined_annotated_ForNeighborhoods.csv")

All.combined@meta.data <- All.combined@meta.data %>% mutate(Timepoint = ifelse(All.combined@meta.data$orig.ident %in% c("SB67_NBM46_AML1_382_CODEX_Mesmer", "SB67_NBM54_AML3_1443_CODEX_Mesmer"), yes = "Post-Therapy", no = "Diagnosis"))

### Use Jupyter notebook to perform neighborhood analysis and return here ###


# import neighborhood csvs
neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_NSM_combined.neighborhood.csv")
fc <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_NSM.combined_neighborhood_fc.csv")

neighborhood_mat <- as.matrix(fc[2:34]) # not all cell types were mapped to

# analyze neighborhood frequency per condition
nbs_aml <- neighborhoods %>% mutate(Sample = ifelse(Sample_Group == "NSM",yes="NSM", no="Diagnosis"))
nbs_aml <- neighborhoods %>% mutate(Sample = case_when(
  (Sample_Group == "NSM") ~ "NSM",
  (orig.ident %in% c("SB67_NBM46_AML1_382_CODEX_Mesmer", "SB67_NBM54_AML3_1443_CODEX_Mesmer")) ~ "Post-Therapy",
  (Sample_Name %in% c("AML1_Dx", "AML2_Dx", "AML3_Dx")) ~ "Diagnosis"
))
# calculate frequency
df_frequency <- xtabs(~ neighborhood10  + Sample, data = nbs_aml)
df_frequency_df <- as.data.frame.matrix(df_frequency)
df_frequency_df$neighborhood10 <- rownames(df_frequency_df)
df_frequency_df <- df_frequency_df %>% pivot_longer(cols = c(1:3), names_to = "Sample", values_to = "Frequency")

df_frequency_df <- df_frequency_df %>%
  group_by(Sample) %>%
  mutate(total = sum(Frequency))
  
df_frequency_df$rel_freq <-  df_frequency_df$Frequency / df_frequency_df$total
df_frequency_df %>% group_by(Sample) %>% summarise(sum(rel_freq)) # just confirm it adds up to one





# Figure 7G - Neighborhood Membership Frequency by Sample Type and Timepoint ----
neighborhood_order <- c(0,2,3,4,5,6,7,8,9,10,11,13, 1,12,14)
df_frequency_df$neighborhood10 <- factor(df_frequency_df$neighborhood10, levels = neighborhood_order)
levels(df_frequency_df$neighborhood10) <- c("CN1", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8", "CN9", "CN10", "CN11", "CN12", "CN14", "CN2", "CN13", "CN15")
df_frequency_df  %>%  ggplot(mapping = aes(x = neighborhood10, y = rel_freq, fill = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  xlab("Neighborhood") +
  ylab("Frequency") + scale_fill_manual(values = c("#FA8072", "#89CFF0", "#734F96")) + theme_minimal() + RotatedAxis()+
  ggtitle("Frequency of Neighborhood Membership by Mutation Status and Timepoint")

# Figure 7F - Plot neighborhoods heatmap -----
neighborhood_mat <- as.matrix(fc[2:34])
rownames(neighborhood_mat) <- rownames(fc)
neighborhood_order <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)

neighborhood_counts <- as.numeric(table(factor(neighborhoods$neighborhood10, levels = neighborhood_order)))

ha = HeatmapAnnotation(cell_counts = anno_barplot(x = neighborhood_counts,
                                                  bar_width = 1, 
                                                  gp = gpar(col = "white", fill = "grey"), 
                                                  border = TRUE,
                                                  axis_param = list(at = c(0, 2e4, 4e4),
                                                                    labels = c("0", "20k", "40k")),
                                                  height = unit(4, "cm"), width = unit(1,"cm"), show_annotation_name = FALSE), annotation_label = c("Counts"),annotation_name_side = 'top', annotation_name_rot = 360, annotation_name_align = TRUE,
                       border = TRUE, which = 'row')

# Heatmap(neighborhood_mat, right_annotation = ha, border = TRUE)

htmp <- Heatmap(neighborhood_mat, name = "mat", rect_gp = gpar(col = "black", lwd = 2), 
                column_title = "Bone Marrow Neighborhood Enrichment",right_annotation =  ha, column_names_gp = grid::gpar(fontsize = 14), row_title_gp = grid::gpar(fontsize = 13))
draw(htmp, heatmap_legend_side="left")

# Supplemental Figure S7E AML Bone Proximity Boxplot----
ranks <- read_csv("~/Documents/NBM_Microenvironment/NBM_Atlas_CODEX_Samples/Non-Cell Microenvironment Analysis/Most_Uptodate/combined_v3_aml_neighbourhood.csv")
ranks %>% dplyr::filter(ranks$structure == "bone") %>% ggplot(aes(x = as.factor(neighbourhood+1), y = normalized_rank, fill = structure))  + geom_boxplot() + scale_fill_manual(values = c("#FAA59E","#CADCEB","#C2E7B9","#E8DBEC","#FED194","#ffe2db")) + scale_linetype_manual(values=c("solid", "longdash")) + theme_minimal()


