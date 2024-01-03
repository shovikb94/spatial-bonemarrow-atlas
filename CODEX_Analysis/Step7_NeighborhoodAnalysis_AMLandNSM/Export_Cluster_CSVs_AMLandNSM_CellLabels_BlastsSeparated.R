neighborhoods <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/NeighborhoodsOutput/AML_only_combined.neighborhood_blasts_separated.csv")
annotations_col <- "classified_cluster_anno_l2_separate"
integrated_seurat <- All.combined_mapped_blastslabeled
cluster_factor <- as.numeric(as.factor(integrated_seurat$classified_cluster_anno_l2_separate))
integrated_seurat$cluster_levels <- cluster_factor
conversion_table <- data.frame(
  Level = levels(as.factor(integrated_seurat$classified_cluster_anno_l2_separate)),
  Numeric = levels(factor(as.numeric(factor(cluster_factor))))
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
  write.csv(annotations_df, file = paste0("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas_NeighborhoodMasks/Masks_ForManuscript/AMLandNSM_ClusterMask_CSVs/Unfiltered_",name,"_",annotations_col,"_","blasts_separated_10252023.csv"))
}

# troubleshoot
sub(".*_", "", colnames(subset(integrated_seurat, subset = orig.ident == paste0(name[1]))))


