library(ComplexHeatmap)
library(readr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
setwd("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6")

# utility functions
require(dplyr)
require(stringr)

coalesce_join <- function(x, 
                          y, 
                          by = NULL, 
                          keep = c("left", "right"), # "left" means keep value from left table if values exist in both tables.
                          suffix = c(".x", ".y"), # Same as the suffix argument in dplyr joins. 
                          join = c("full_join","left_join", "right_join", "inner_join") # Choose a join type from the list. The default is full_join.
) { 
  keep = match.arg(keep) 
  join = match.arg(join) 
  join = match.fun(join) # Confirm the join argument is in the list and match the string to the function
  
  # Depends on the keep argument, overwrite the duplicate value
  # If keep = "left", the value from the left table will be kept, vice versa.
  if (keep == "left") suffix_ = suffix else suffix_ = rev(suffix)
  
  join(x, y, by = by, suffix = suffix) %>% 
    mutate(
      across( # Apply the coalesce function to all overlapped columns
        ends_with(suffix_[1]), # Select columns ended with .x if keep = "left"; or .y if keep = "right"
        ~coalesce(.,
                  get(str_replace(cur_column(), suffix_[1], suffix_[2])) # Replace .x in var.x with .y to generate var.y, if keep = "left"; or vice versa.
        ),
        .names = "{str_remove(.col, suffix_[1])}" # Remove the suffix from the combined columns
      ),
      .keep = "unused") # Remove the temporary columns ended with suffix
}



# read in objects
AML1_183_mapped <- readRDS("objects/AML1_183_mapped.RDS")
AML1_382_mapped <- readRDS("objects/AML1_382_mapped.RDS")
AML2_191_mapped <- readRDS("objects/AML2_191_mapped.RDS")
AML3_1329_mapped <- readRDS("objects/AML3_1329_mapped.RDS")
AML3_1443_mapped <- readRDS("objects/AML3_1443_mapped.RDS")
NSM.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/NSM.combined.mapped_unfiltered.RDS")

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

# export for neighborhood analysis ----

## combined objects 
All.combined_mapped <- merge(NSM.combined,c(AML1_183_mapped, AML1_382_mapped, AML2_191_mapped, AML3_1329_mapped, AML3_1443_mapped))

## change cluster anno l2 to blasts
Blast_ids <- colnames(subset(All.combined_mapped, subset = Adjusted_Class == "NPM1_Mutant")) # Blast Rename
All.combined_mapped$classified_cluster_anno_l2[Blast_ids] <- "NPM1 Mutant Blast"

## remove artifacts
All.combined_mapped <- subset(All.combined_mapped, subset = classified_cluster_anno_l2 != "Artifact" & classified_cluster_anno_l2 != "CD44+ Undetermined" & classified_cluster_anno_l2 != "Undetermined" & classified_cluster_anno_l2 != "Autofluorescent")


## Add metadata columns and write annotated CSV for neighborhood analysis
All.combined_mapped@meta.data$x.coord_rev <- -(All.combined_mapped$x.coord)
All.combined_mapped@meta.data$y.coord_rev <- -(All.combined_mapped$y.coord)
All.combined_mapped@meta.data <- mutate(All.combined_mapped@meta.data, Region = "reg001")
write.csv(x = All.combined_mapped@meta.data, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/For_Neighborhoods/NSM_AML_combined_annotated_ForNeighborhoods_062723_test.csv")

# Separate blasts
## change cluster anno l2 to blasts
All.combined_mapped$classified_cluster_anno_l2_separate <- All.combined_mapped$classified_cluster_anno_l2
Dx_Blast_ids <- colnames(subset(All.combined_mapped, subset = Adjusted_Class == "NPM1_Mutant" & Sample_Timepoint == "Dx")) # Blast Rename
Post_Tx_Blast_ids <- colnames(subset(All.combined_mapped, subset = Adjusted_Class == "NPM1_Mutant" & Sample_Timepoint == "Post-Tx")) # Blast Rename

All.combined_mapped$classified_cluster_anno_l2_separate[Dx_Blast_ids] <- "Dx NPM1 Mutant Blast"
All.combined_mapped$classified_cluster_anno_l2_separate[Post_Tx_Blast_ids] <- "Post_Tx NPM1 Mutant Blast"

saveRDS(All.combined_mapped, file = "objects/All.combined_mapped_blasts_separated.RDS")

## Add metadata columns and write annotated CSV for neighborhood analysis for blasts separated
All.combined_mapped <- SetIdent(All.combined_mapped, value = "classified_cluster_anno_l2")
small_clusters <- names(which(table(Idents(All.combined_mapped)) < 10))
All.combined_mapped_filtered <- subset(All.combined_mapped, idents = small_clusters, invert = TRUE)
write.csv(x = All.combined_mapped_filtered@meta.data, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/For_Neighborhoods/NSM_AML_combined_annotated_ForNeighborhoods_112023_blast_separated.csv")

AMLonly.combined_mapped <- subset(All.combined_mapped, Sample_Group != "NSM")
# remove cell types with less than 10 cells
AMLonly.combined_mapped <- SetIdent(AMLonly.combined_mapped, value = "classified_cluster_anno_l2")
small_clusters <- names(which(table(Idents(AMLonly.combined_mapped)) < 10))

AMLonly.combined_mapped <- subset(AMLonly.combined_mapped, idents = small_clusters, invert = TRUE)

write.csv(x = AMLonly.combined_mapped@meta.data, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/For_Neighborhoods/AML_only_combined_annotated_ForNeighborhoods_112023_blast_separated.csv")

# also do dx and post-tx separately
AMLonly.dx_mapped <- subset(AMLonly.combined_mapped, Sample_Timepoint == "Dx")
AMLonly.posttx_mapped <- subset(AMLonly.combined_mapped, Sample_Timepoint == "Post-Tx")
write.csv(x = AMLonly.dx_mapped@meta.data, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/For_Neighborhoods/AML_only_combined_annotated_ForNeighborhoods_112023_blast_separated_DxOnly.csv")
write.csv(x = AMLonly.posttx_mapped@meta.data, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/For_Neighborhoods/AML_only_combined_annotated_ForNeighborhoods_112023_blast_separated_PostTxOnly.csv")



# repeat for NSM (need to have mapped per Step 5)
NSM.combined <- subset(NSM.combined, subset = classified_cluster_anno_l2 != "Artifact" & classified_cluster_anno_l2 != "CD44+ Undetermined" & classified_cluster_anno_l2 != "Undetermined" & classified_cluster_anno_l2 != "Autofluorescent")

# remove cell types with less than 10 cells
NSM.combined <- SetIdent(NSM.combined, value = "classified_cluster_anno_l2")
small_clusters <- names(which(table(Idents(NSM.combined)) < 10))

NSM.combined <- subset(NSM.combined, idents = small_clusters, invert = TRUE)

# Add metadata columns and write annotated CSV for neighborhood analysis
NSM.combined@meta.data$x.coord_rev <- -(NSM.combined$x.coord)
NSM.combined@meta.data$y.coord_rev <- -(NSM.combined$y.coord)
NSM.combined@meta.data <- mutate(NSM.combined@meta.data, Region = "reg001")

write.csv(x = NSM.combined@meta.data, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/For_Neighborhoods/NSM_filtered_ForNeighborhoods.csv")





### Use Jupyter notebook to perform neighborhood analysis and return here ###

# Generate Figures from AML Mapping ----
All.combined_mapped <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/objects/All.combined_mapped_blasts_separated.RDS") # minor note - the export for neighborhood analysis was done with filtering clusters <10 cells, but the rest of the analysis was done unfiltered since these analyses were largely agnostic to cell type (e.g.. pre and post therapy marker staining patterns in NPM1 blasts)


# import neighborhood csvs for Dx ----
neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only_combined.neighborhood_blasts_separated_DxOnly.csv")
fc <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only.combined_neighborhood_fc_blasts_separated_DxOnly.csv")
neighborhoods$unique_CellID <- paste0(neighborhoods$Sample_Name,"_",neighborhoods$CellID)
aml_nbs_tojoin <- data.frame(neighborhoods$unique_CellID, neighborhoods$neighborhood10)
colnames(aml_nbs_tojoin) <- c("unique_CellID", "neighborhood10")
All.combined_mapped@meta.data$unique_CellID <- paste0(All.combined_mapped@meta.data$Sample_Name,"_",All.combined_mapped@meta.data$CellID)

All.combined_mapped@meta.data <- All.combined_mapped@meta.data %>% left_join(aml_nbs_tojoin, by = "unique_CellID")

# write out for visualization

neighborhood_mat <- as.matrix(fc[2:31]) # not all cell types were mapped to

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
df_frequency_df <- df_frequency_df %>% pivot_longer(cols = c(1:2), names_to = "Sample", values_to = "Frequency") # should be c(1:3) if using NSM and AML combined heatmap

df_frequency_df <- df_frequency_df %>%
  group_by(Sample) %>%
  mutate(total = sum(Frequency))
  
df_frequency_df$rel_freq <-  df_frequency_df$Frequency / df_frequency_df$total
df_frequency_df %>% group_by(Sample) %>% summarise(sum(rel_freq)) # just confirm it adds up to one

# Figure 7C GATA1 figure ----
VlnPlot(subset(All.combined_mapped, classified_cluster_anno_l2 == "NPM1 Mutant Blast"), features = c("codex_GATA1"), group.by = "Sample_Timepoint", pt.size = 0) 
## Test for markers enriched pre and post therapy, GATA1 text reference ----
data_dx <- as.data.frame(t(as.data.frame(GetAssayData(subset(All.combined_mapped, subset = Sample_Timepoint == "Dx" & Adjusted_Class == "NPM1_Mutant"), slot = "data", assay = "CODEX"))))
data_mrd <- as.data.frame(t(as.data.frame(GetAssayData(subset(All.combined_mapped, subset = Sample_Timepoint == "Post-Tx" & Adjusted_Class == "NPM1_Mutant"), slot = "data", assay = "CODEX"))))

t.test(data_mrd$GATA1,data_dx$GATA1) # t=36.553, p < 2.2e-16

# Figure 7D MSC Frequency in AML vs. NSM ----

# Rename orig.idents 
levels(as.factor(All.combined_mapped$orig.ident))
#orig.ident_annotations <- c("AML2_Dx","AML1_PostTx","NSM3_1996", "NSM1_1720", "NSM2_1086", "AML1_Dx", "AML3_Dx", "AML3_PostTx")
orig.ident_annotations <- c("Dx","PostTx","NSM", "NSM", "NSM", "Dx", "Dx", "PostTx")
All.combined_mapped <- SetIdent(All.combined_mapped, value = "orig.ident")
names(orig.ident_annotations) <- levels(as.factor(All.combined_mapped$orig.ident))
All.combined_mapped <- RenameIdents(All.combined_mapped, orig.ident_annotations)
All.combined_mapped@meta.data["Sample_Timepoint"] <- All.combined_mapped@active.ident 

levels(as.factor(All.combined_mapped$orig.ident))
orig.ident_annotations <- c("AML2_Dx","AML1_PostTx","NSM3_1996", "NSM1_1720", "NSM2_1086", "AML1_Dx", "AML3_Dx", "AML3_PostTx")
All.combined_mapped <- SetIdent(All.combined_mapped, value = "orig.ident")
names(orig.ident_annotations) <- levels(as.factor(All.combined_mapped$orig.ident))
All.combined_mapped <- RenameIdents(All.combined_mapped, orig.ident_annotations)
All.combined_mapped@meta.data["Sample_Name"] <- All.combined_mapped@active.ident 


## make a bar plot just showing the sample groups (Dx vs. PostTx vs. NSM) MSC frequencies and perform t.test -----
All_ct_freqs <- All.combined_mapped@meta.data %>%
  group_by(Sample_Name, classified_cluster_anno_l2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

All_ct_freqs <- All_ct_freqs %>% 
  mutate(Sample_Timepoint = case_when(
    grepl("Dx$", Sample_Name) ~ "Dx",          # Ends with Dx
    grepl("PostTx$", Sample_Name) ~ "PostTx",  # Ends with PostTx
    grepl("^NSM", Sample_Name) ~ "NSM",        # Starts with NSM
    TRUE ~ as.character(Sample_Name)           # Default case
  ))

All_ct_freqs_toplot <- All_ct_freqs %>% group_by(Sample_Timepoint, classified_cluster_anno_l2) %>% summarise(freq=freq, mean_freq = mean(freq), sd = sd(freq)) 
All_ct_freqs_toplot <- All_ct_freqs_toplot %>% dplyr::filter(classified_cluster_anno_l2 == "MSC" | classified_cluster_anno_l2 == "THY1+ MSC")

All_ct_freqs_toplot$Sample_Timepoint <- factor(All_ct_freqs_toplot$Sample_Timepoint, levels = c("Dx", "PostTx", "NSM"))
p1 <- All_ct_freqs_toplot %>% ggplot(aes(x = Sample_Timepoint, y=mean_freq, fill = classified_cluster_anno_l2)) + 
  geom_bar(stat="identity", position = "dodge", color = 'black') + 
  geom_errorbar(aes(ymin=mean_freq-sd, ymax=mean_freq+sd), width=.2,
                position=position_dodge(.9)) + scale_fill_manual(values = c("#89CFF0", "#FA8072")) + theme_minimal() + RotatedAxis() # Figure 7D graph

ggsave(p1, device = "pdf",filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/Figures/Figure7D_MSC_SubTypeFrequency_By_SampleTimepoint_123123.pdf", height = 5, width = 7.5)

## repeat by individual sample for revision -----
# make a bar plot just showing the MSC frequencies and perform t.test
All_ct_freqs <- All.combined_mapped@meta.data %>%
  group_by(Sample_Name, classified_cluster_anno_l2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

All_ct_freqs <- All_ct_freqs %>% 
  mutate(Sample_Timepoint = case_when(
    grepl("Dx$", Sample_Name) ~ "Dx",          # Ends with Dx
    grepl("MRD$", Sample_Name) ~ "PostTx",  # Ends with PostTx
    grepl("^NSM", Sample_Name) ~ "NSM",        # Starts with NSM
    TRUE ~ as.character(Sample_Name)           # Default case
  ))

All_ct_freqs_toplot <- All_ct_freqs %>% group_by(Sample_Name, classified_cluster_anno_l2) %>% summarise(freq=freq, mean_freq = mean(freq), sd = sd(freq)) 
All_ct_freqs_toplot <- All_ct_freqs_toplot %>% dplyr::filter(classified_cluster_anno_l2 == "MSC" | classified_cluster_anno_l2 == "THY1+ MSC")

All_ct_freqs_toplot$Sample_Name <- factor(All_ct_freqs_toplot$Sample_Name, levels = c("NSM1_1720","NSM2_1086", "NSM3_1996","AML1_Dx", "AML1_PostTx", "AML2_Dx", "AML3_Dx", "AML3_PostTx"))
p2 <- All_ct_freqs_toplot %>% ggplot(aes(x = Sample_Name, y=mean_freq, fill = classified_cluster_anno_l2)) + 
  geom_bar(stat="identity", position = "dodge", color = 'black') + 
  geom_errorbar(aes(ymin=mean_freq-sd, ymax=mean_freq+sd), width=.2,
                position=position_dodge(.9)) + scale_fill_manual(values = c("#89CFF0", "#FA8072")) + theme_minimal() + RotatedAxis() 



aml_sample_names <- c("AML1_Dx", "AML1_MRD", "AML2_Dx", "AML3_Dx", "AML3_MRD")
NSM_sample_names <- c("NSM_1086", "NSM_1720", "NSM_1996")

aml_sample_names <- c("AML1_Dx", "AML1_MRD", "AML2_Dx", "AML3_Dx", "AML3_MRD")
NSM_sample_names <- c("NSM_1086", "NSM_1720", "NSM_1996")

## Welch Two Sample t-test for seeing whether the frequency distribution is different between AML and NSM for both Adipo and THY1+ MSCs ----

t.test(subset(All_ct_freqs, subset = Sample_Timepoint != "NSM" & classified_cluster_anno_l2 == "MSC")$freq, subset(All_ct_freqs, subset = Sample_Timepoint == "NSM" & (classified_cluster_anno_l2 == "MSC"))$freq)
# p = 0.00488, t=4.6201, 0.0134 vs. 0.0043
t.test(subset(All_ct_freqs, subset = Sample_Timepoint != "NSM" & (classified_cluster_anno_l2 == "THY1+ MSC"))$freq, subset(All_ct_freqs, subset = Sample_Timepoint == "NSM" & (classified_cluster_anno_l2 == "THY1+ MSC"))$freq)
# p = 0.03651, t=2.7096,  0.0043 vs. 0.00138




# Figure 7G - Neighborhood Membership Frequency by Sample Type and Timepoint ----
neighborhood_order <- c(0,2,3,4,5,6,7,8,9,10,11,13, 1,12,14)
df_frequency_df$neighborhood10 <- factor(df_frequency_df$neighborhood10, levels = neighborhood_order)
levels(df_frequency_df$neighborhood10) <- c("CN1", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8", "CN9", "CN10", "CN11", "CN12", "CN14", "CN2", "CN13", "CN15")
df_frequency_df  %>%  ggplot(mapping = aes(x = neighborhood10, y = rel_freq, fill = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  xlab("Neighborhood") +
  ylab("Frequency") + scale_fill_manual(values = c("#FA8072", "#89CFF0", "#734F96")) + theme_minimal() + RotatedAxis()+
  ggtitle("Frequency of Neighborhood Membership by Mutation Status and Timepoint")

# Figure 7F - Plot neighborhoods heatmap for Dx -----
neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only_combined.neighborhood_blasts_separated_DxOnly.csv")
fc <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only.combined_neighborhood_fc_blasts_separated_DxOnly.csv")

neighborhood_mat <- as.matrix(fc[2:31])
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


htmp <- Heatmap(neighborhood_mat, name = "mat", rect_gp = gpar(col = "black", lwd = 2), column_names_rot = 45,
                column_title = "Dx AML Neighborhood Enrichment",right_annotation =  ha, column_names_gp = grid::gpar(fontsize = 14), row_title_gp = grid::gpar(fontsize = 13))
draw(htmp, heatmap_legend_side="left")

neighborhoods$neighborhood_named <- neighborhood_names_zeroindex[as.character(neighborhoods$neighborhood10)]

nbs_names <- c("Mature Myeloid", "Erythroid/Myeloid", "PC/Arteriolar", "Erythroid", "Arteriolar", "Erythroid", "Lymphoid", "Erythroid/Myeloid/Lymphoid", "Early Myeloid", "Myeloid/Lymphoid", "Intermediate Myeloid", "Erythroid/Myeloid/Lymphoid", "Erythroid/Myeloid", "Early Myeloid", "Peri-Arterolar Lymphoid")
nbs_cols <- c("#A8D37F", "#51C8EB", "#FF00FF", "#FBED24", "#FF0000", "#FBED24", "#89919C", "#BD7CB5", "#4B409A", "#FFA500", "#A15A26", "#BD7CB5", "#51C8EB", "#4B409A", "#FF007F")
names(nbs_cols) <- nbs_names



Dx_nb_names <- c("0" = "Erythroid/Myeloid/Lymphoid", "1" = "Monocyte", "2" = "MSC/Blast Only", "3" = "Erythroid", "4" = "Blast/Mixed", "5" = "Lymphoid 1", "6" = "Erythroid/Myeloid", "7" = "PC/Arteriolar", "8" = "Arterio-Endosteal MSC/Blast", "9" = "Early MP/Arterio-Endosteal", "10" = "Mature Myeloid", "11" = "Intermediate Myeloid 2", "12" = "Intermediate Myeloid 1", "13" = "Lymphoid 2", "14" = "Lymphoid 3") 
neighborhoods$neighborhood_named <- Dx_nb_names[as.character(neighborhoods$neighborhood10)]


## repeat for Post-Therapy AML ----
neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only_combined.neighborhood_blasts_separated_PostTx.csv")
fc <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only.combined_neighborhood_fc_blasts_separated_PostTx.csv")
neighborhood_mat <- as.matrix(fc[2:31])
rownames(neighborhood_mat) <- rownames(fc)
neighborhood_order <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)

neighborhood_counts <- as.numeric(table(factor(neighborhoods$neighborhood10, levels = neighborhood_order)))

ha = HeatmapAnnotation(cell_counts = anno_barplot(x = neighborhood_counts,
                                                  bar_width = 1, 
                                                  gp = gpar(col = "white", fill = "grey"), 
                                                  border = TRUE,
                                                  axis_param = list(at = c(0, 1e4),
                                                                    labels = c("0", "10k")),
                                                  height = unit(4, "cm"), width = unit(1,"cm"), show_annotation_name = FALSE), annotation_label = c("Counts"),annotation_name_side = 'top', annotation_name_rot = 360, annotation_name_align = TRUE,
                       border = TRUE, which = 'row')


htmp <- Heatmap(neighborhood_mat, name = "mat", rect_gp = gpar(col = "black", lwd = 2), column_names_rot = 45, 
                column_title = "Post-Tx AML Neighborhood Enrichment",right_annotation =  ha, column_names_gp = grid::gpar(fontsize = 14), row_title_gp = grid::gpar(fontsize = 13))
draw(htmp, heatmap_legend_side="left")

PostTx_nb_names <- c("0" = "Recovering Erythroid EB", "1" = "Recovering Erythroid Mature 1", "2" = "Recovering Erythroid HSPC/EB", "3" = "Arterio-Endosteal MSC/Blast 1", "4" = "PC/Arteriolar 1", "5" = "Erythroid", "6" = "Arterial/Sinusoidal/Mixed", "7" = "Recovering Myeloid GMP/Early MP", "8" = "Lymphoid/Myeloid", "9" = "Recovering Myeloid Int/Mature 1", "10" = "Recovering Erythroid Mature 2", "11" = "Recovering Myeloid Int/Mature 2", "12" = "PC/Arteriolar 2", "13" = "Arterio-Endosteal MSC/Blast 2", "14" = "Mature Myeloid") 
neighborhoods$neighborhood_named <- PostTx_nb_names[as.character(neighborhoods$neighborhood10)]
p2 <- ggplot(neighborhoods, aes(x = neighborhood_named, fill=orig.ident)) + 
  geom_bar(position="fill", stat="count") + theme_minimal() + RotatedAxis() +
  ylab("Percentage of Total MSCs") + xlab("Sample") # limits=force drops unused levels from the legend

# Supplemental Figure S7E -----
library(cowplot)
plot_grid(p1,p2)

# AML NPM1 Staining Per Cell Type For Revision -----
All.combined_mapped_blasts_separated <- SetIdent(All.combined_mapped_blasts_separated, value = "classified_cluster_anno_l2")
small_clusters <- names(which(table(Idents(All.combined_mapped_blasts_separated)) < 10))
All.combined_mapped__blasts_separated_filtered <- subset(All.combined_mapped_blasts_separated, idents = small_clusters, invert = TRUE)
VlnPlot(subset(All.combined_mapped__blasts_separated_filtered, Sample_Group != "NSM"), features = c("NPM1C"), group.by = "classified_cluster_anno_l2", pt.size = 0, sort = TRUE) + geom_boxplot(outlier.shape = NA) + NoLegend()

VlnPlot(subset(All.combined_mapped__blasts_separated_filtered, Sample_Group != "NSM" & Class %in% c("NPM1_Mutant", "WT")), features = c("NPM1C"), group.by = "Class", pt.size = 0, sort = TRUE) + geom_boxplot(outlier.shape = NA) + NoLegend()

data_mut <- as.data.frame(t(as.data.frame(GetAssayData(subset(All.combined_mapped__blasts_separated_filtered, subset = Sample_Group != "NSM" & Adjusted_Class == "NPM1_Mutant"), slot = "data", assay = "CODEX"))))
data_wt <- as.data.frame(t(as.data.frame(GetAssayData(subset(All.combined_mapped__blasts_separated_filtered, subset = Sample_Group != "NSM" & Adjusted_Class == "WT"), slot = "data", assay = "CODEX"))))

t.test(data_mut$NPM1C,data_wt$NPM1C) # t=165.32, p < 2.2e-16

data_mut <- as.data.frame(t(as.data.frame(GetAssayData(subset(All.combined_mapped__blasts_separated_filtered, subset = Sample_Group != "NSM" & Adjusted_Class == "NPM1_Mutant"), slot = "counts", assay = "CODEX"))))
data_wt <- as.data.frame(t(as.data.frame(GetAssayData(subset(All.combined_mapped__blasts_separated_filtered, subset = Sample_Group != "NSM" & Adjusted_Class == "WT"), slot = "counts", assay = "CODEX"))))

t.test(data_mut$NPM1C,data_wt$NPM1C) # t=219.59, p < 2.2e-16


# Analyze HIF1a levels Supplemental Figure S7H ----
cal2_cols_npm1 <- c(cal2_cols, "#FA8072")
names(cal2_cols_npm1) <- c(cell_order, "NPM1 Mutant Blast")

p1 <- VlnPlot(subset(All.combined_mapped,classified_cluster_anno_l2 %in% c("Mature Myeloid", "Intermediate Myeloid", "Early Myeloid Progenitor","GMP/Myeloblast","GMP", "NPM1 Mutant Blast")), features = c("codex_HIF1A"), slot = "data", pt.size = 0 , cols = cal2_cols, group.by = "classified_cluster_anno_l2", sort = TRUE) + NoLegend()
#ggsave(p1, device = "pdf",filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/Figures/SuppS6C_NPM1_HIF1A_VlnPlot.pdf", height = 5, width = 7.5)

# BCL2 levels Supplemental Figure S7I ----
p1 <- VlnPlot(subset(All.combined_mapped,classified_cluster_anno_l2 %in% c("NPM1 Mutant Blast")), split.by = "Sample_Name", cols = c("#87CEEB", "#20639B", "#3CAEA3","#F6D55C", "#ED553B"), features = c("codex_BCL2"), slot = "data", pt.size = 0, group.by = "classified_cluster_anno_l2", sort = TRUE) 
ggsave(p1, device = "pdf",filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/Figures/SuppS6C_NPM1_BCL2_VlnPlot.pdf", height = 5, width = 7.5)

# Complex IV levels Supplemental Figure S7J ----
p1 <- VlnPlot(subset(All.combined_mapped,classified_cluster_anno_l2 %in% c("NPM1 Mutant Blast")), group.by = "Sample_Name", features = c("codex_OXPHOS"), slot = "data", pt.size = 0, cols = c("#87CEEB", "#20639B", "#3CAEA3","#F6D55C", "#ED553B"),  sort = FALSE) 
ggsave(p1, device = "pdf",filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/Figures/SuppS6C_NPM1_ComplexIV_VlnPlot.pdf", height = 5, width = 7.5)

# t-tests for supplemental figure S7 ----
data_dx <- as.data.frame(t(as.data.frame(GetAssayData(subset(All.combined_mapped, subset = Sample_Timepoint == "Dx" & Adjusted_Class == "NPM1_Mutant"), slot = "data", assay = "CODEX"))))
data_mrd <- as.data.frame(t(as.data.frame(GetAssayData(subset(All.combined_mapped, subset = Sample_Timepoint == "Post-Tx" & Adjusted_Class == "NPM1_Mutant"), slot = "data", assay = "CODEX"))))

t.test(data_mrd$GATA1,data_dx$GATA1) # t=36.553, p < 2.2e-16
t.test(data_mrd$BCL2,data_dx$BCL2) # t=-35.371, p < 2.2e-16
t.test(data_mrd$OXPHOS,data_dx$OXPHOS) # t=24.722, p < 2.2e-16


# Supplemental Figure S7G AML Bone Proximity Boxplot----
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




