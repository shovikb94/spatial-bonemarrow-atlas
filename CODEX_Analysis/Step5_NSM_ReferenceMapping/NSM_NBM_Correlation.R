# Figure S8A - Plot NSM neighborhoods heatmap -----
#df3 <- anti_join(All.combined_mapped@meta.data,neighborhoods, by = "unique_CellID") %>% bind_rows(neighborhoods)
fc <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/NSMonly_neighborhood_fc.csv")
neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/NSMonly.neighborhood.csv")

#neighborhoods$unique_CellID <- paste0(neighborhoods$Sample_Name,"_",neighborhoods$CellID)
#nsm_nbs_tojoin <- data.frame(neighborhoods$unique_CellID, neighborhoods$neighborhood10)
#colnames(nsm_nbs_tojoin) <- c("unique_CellID", "neighborhood10")

#All.combined_mapped@meta.data <- coalesce_join(join = "left", x = All.combined_mapped@meta.data, y=nsm_nbs_tojoin, by = "unique_CellID")


neighborhood_mat <- as.matrix(fc[2:29])
rownames(neighborhood_mat) <- rownames(fc)
neighborhood_order <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
NSM_nb_names <- c("0" = "Erythroid 1", "1" = "Erythroid/Myeloid", "2" = "Early Myeloid/Arteriolar", "3" = "Mature Myeloid 1", "4" = "Myeloid Multiple Stages", "5" = "Mature Myeloid 2", "6" = "PC/Arteriolar", "7" = "Myeloid/Lymphoid 1","8" = "Erythroid 2", "9" = "Erythroid 3", "10" = "Vascular/MSC 1",
                  "11" = "Early Myeloid/Endosteal", "12" = "Intermediate Myeloid", "13" = "Erythroid/Myeloid/Lymphoid", "14" = "Vascular/MSC 2")

neighborhood_counts <- as.numeric(table(factor(neighborhoods$neighborhood10, levels = neighborhood_order)))
rownames(neighborhood_mat) <- NSM_nb_names



ha = HeatmapAnnotation(cell_counts = anno_barplot(x = neighborhood_counts,
                                                  bar_width = 1, 
                                                  gp = gpar(col = "white", fill = "grey"), 
                                                  border = TRUE,
                                                  axis_param = list(at = c(0, 1e4, 2e4),
                                                                    labels = c("0","10k","20k")),
                                                  height = unit(4, "cm"), width = unit(1,"cm"), show_annotation_name = FALSE), annotation_label = c("Counts"),annotation_name_side = 'top', annotation_name_rot = 360, annotation_name_align = TRUE,
                       border = TRUE, which = 'row')


# Heatmap(neighborhood_mat, right_annotation = ha, border = TRUE)

htmp <- Heatmap(neighborhood_mat, name = "mat", rect_gp = gpar(col = "black", lwd = 2), 
                column_title = "Bone Marrow Neighborhood Enrichment",right_annotation =  ha, column_names_gp = grid::gpar(fontsize = 14), row_title_gp = grid::gpar(fontsize = 13))
draw(htmp, heatmap_legend_side="left")




# correlation of NSM neighborhoods and NBM ----

# make general names for nsm nbs
NSM_nb_names_nonunique <- c("0" = "Erythroid", "1" = "Erythroid/Myeloid", "2" = "Early Myeloid/Arteriolar", "3" = "Mature Myeloid", "4" = "Myeloid Multiple Stages", "5" = "Mature Myeloid", "6" = "PC/Arteriolar", "7" = "Myeloid/Lymphoid","8" = "Erythroid", "9" = "Erythroid", "10" = "Vascular/MSC",
                            "11" = "Early Myeloid/Endosteal", "12" = "Intermediate Myeloid", "13" = "Erythroid/Myeloid/Lymphoid", "14" = "Vascular/MSC")




neighborhoods_nbm <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/output/neighborhood.csv")
fc_nbm <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/output/neighborhood_fc.csv")
neighborhood_nbm_mat <- as.matrix(fc_nbm[2:33])

neighborhood_names_zeroindex_unique <- c('0' = "Mature Myeloid",
                                         "1" = "Erythroid/Myeloid_1",
                                         "2" = "PC/Arteriolar",
                                         "3" = "Erythroid_1",
                                         "4" = "Vascular/Endosteal/MSC",
                                         "5" = "Erythroid_2",
                                         "6" = "Lymphoid",
                                         "7" = "Erythroid/Myeloid/Lymphoid_1",
                                         "8" = "Early Myeloid/Endosteal",
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
                                  "4" = "Vascular/Endosteal/MSC",
                                  "5" = "Erythroid",
                                  "6" = "Lymphoid",
                                  "7" = "Erythroid/Myeloid/Lymphoid",
                                  "8" = "Early Myeloid/Endosteal",
                                  "9" = "Myeloid/Lymphoid",
                                  "10" = "Intermediate Myeloid",
                                  "11" = "Erythroid/Myeloid/Lymphoid",
                                  "12" = "Erythroid/Myeloid",
                                  "13" = "Early Myeloid/Arteriolar",
                                  "14" = "Peri-Arterolar Lymphoid")

rownames(neighborhood_nbm_mat) <- neighborhood_names_zeroindex
rownames(neighborhood_mat) <- NSM_nb_names_nonunique
# identify common nbs
common_neighborhoods <- intersect(rownames(neighborhood_mat), rownames(neighborhood_nbm_mat))

# Identify common cell types
common_cell_types <- intersect(colnames(neighborhood_mat), colnames(neighborhood_nbm_mat))


# Initialize a vector to store correlation results
correlation_results <- numeric(length(common_neighborhoods))
names(correlation_results) <- common_neighborhoods

# Iterate through each common neighborhood
# Check if there is only one common cell type
if (length(common_cell_types) == 1) {
  warning("Only one common cell type found. Correlation cannot be computed.")
} else {
  # Iterate through each common neighborhood
  for (neighborhood in common_neighborhoods) {
    # Extract rows corresponding to the current neighborhood
    row1 <- neighborhood_mat[neighborhood, common_cell_types, drop = FALSE]
    row2 <- neighborhood_nbm_mat[neighborhood, common_cell_types, drop = FALSE]
    
    # Ensure that the rows are vectors
    row1 <- as.numeric(row1)
    row2 <- as.numeric(row2)
    
    # Calculate correlation
    # Note: Using use = "complete.obs" to handle missing values, if any
    correlation_results[neighborhood] <- cor(row1, row2, use = "complete.obs")
  }
}

# Print the result
print(correlation_results)
median(correlation_results) # R = 0.71

# Convert the correlation results to a data frame
correlation_df <- data.frame(
  Neighborhood = names(correlation_results),
  Correlation = correlation_results
)

# Calculate the median correlation coefficient
median_correlation <- median(correlation_results, na.rm = TRUE)

# Supplemental Figure S8B NSM vs NBM Correlation ----
p <- ggplot(correlation_df, aes(x = Neighborhood, y = Correlation)) +
  geom_point() + # Add points
  geom_hline(yintercept = median_correlation, linetype = "dashed", color = "blue") + # Add a line for the median
  geom_text(aes(label = paste("Median:", round(median_correlation, 3))),
            x = Inf, y = median_correlation, hjust = 1.1, vjust = 2, 
            color = "blue", size = 3) + # Add text for the median
  theme_minimal() +
  labs(title = "Correlation Coefficients by Neighborhood",
       x = "Neighborhood", y = "Correlation Coefficient") + RotatedAxis()

# Print the plot
print(p)

ggsave(p, device = "pdf", filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/AML_NSM_NeighborhoodAnalysis_Step7/Figures/NSM_NBM_neighborhood_correlation.pdf", height = 4, width = 6)
