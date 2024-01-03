
neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/output/neighborhood.csv")
annotations_col <- "neighborhood10"
integrated_seurat <- neighborhoods
cluster_factor <- as.numeric(as.factor(integrated_seurat$neighborhood10))
integrated_seurat$cluster_levels <- cluster_factor
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
conversion_table <- data.frame(
  Level = levels(as.factor(integrated_seurat$neighborhood10)),
  Numeric = levels(factor(as.numeric(factor(cluster_factor)))), nbs_name = neighborhood_names_zeroindex
)
conversion_table
sample_names <- levels(as.factor(integrated_seurat$orig.ident))
for(name in sample_names){
  print(name)
  print("Subseting")
  sample_subset <- subset(integrated_seurat, subset=orig.ident==paste0(name))
  annotations_df <- as.data.frame(sample_subset$cluster_levels)
  table(unique(colnames(sample_subset)))
  rownames(annotations_df) <- sample_subset$CellID
  annotations_df$CellID <- sample_subset$CellID
  colnames(annotations_df) <- c("Annotation", "CellID")
  print("Writing Annotations")
  write.csv(annotations_df, file = paste0("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas_NeighborhoodMasks/Masks_ForManuscript/NBM_Neighborhood_ClusterCSVs/Filtered_",name,"_",annotations_col,"_","_12312023.csv"))
}

# troubleshoot
sub(".*_", "", colnames(subset(integrated_seurat, subset = orig.ident == paste0(name[1]))))


