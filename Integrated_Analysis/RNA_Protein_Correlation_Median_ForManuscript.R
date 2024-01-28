# script to measure scRNA and CODEX correlation in the NBM manuscript

library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)


# load and preprocess objects ----


## read in scRNA-Seq object 
combined <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Seurat/combined_corrected_NKTremoved.RDS")
combined <- subset(combined, subset = cluster_anno_l2 != "RNAlo MSC") # remove RNAlo MSC
combined$cluster_anno_l2 <- droplevels(combined$cluster_anno_l2) # drop empty factor levels for cla2
table(combined$cluster_anno_l2)

## read in CODEX object 
CODEX <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.filtered_FINAL.RDS")


## read integration table

integration_table <- read_csv("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Integration_BatchCorrection/CODEX_scRNA_labels_integration_harmonized.csv")



# read protein rna conversion table
conversion_table <- read_csv("/mnt/isilon/tan_lab/sussmanj/Temp/CODEX_RNA_Integration/Protein_Conversions.csv")

# try to perform the conversion ----

library(Seurat)
library(dplyr)
library(ggplot2)

# Assuming 'combined' and 'codex_seurat' are your Seurat objects for scRNA-Seq data and CODEX data, respectively
# Also assuming 'conversion_table' and 'integration_table' are loaded as data frames

# Create a named vector for marker name conversion between scRNA-seq and CODEX assays
names_conversion <- setNames(conversion_table$`Protein name`, conversion_table$`RNA name`)


# Aggregate RNA expression data by scRNA-Seq cell type
rna_agg <- AverageExpression(combined, slot = 'data', group.by = 'cluster_anno_l2', assays = "RNA", features = c(names(names_conversion)), return.seurat = FALSE)$RNA
# Extracting the data slot which contains normalized data
data_matrix <- combined@assays$RNA@data
gene_list <- names(names_conversion)

# Subset the data matrix to only include the specified genes
subset_data <- data_matrix[gene_list, ]

# Extract cluster information from metadata
combined <- SetIdent(combined, value = "cluster_anno_l2")
clusters <- Idents(combined)

# Create a data frame for calculation
df <- data.frame(expression = as.vector(subset_data),
                 gene = rep(rownames(subset_data), ncol(subset_data)),
                 cluster = rep(clusters, each = nrow(subset_data)))

# Calculate median expression by gene and cluster
median_expression <- df %>%
  group_by(gene, cluster) %>%
  summarize(median_expr = median(expression, na.rm = TRUE))

# View the results
print(median_expression)

wide_format <- median_expression %>%
  pivot_wider(names_from = cluster, values_from = median_expr) 
rna_rownames <- wide_format$gene
wide_format <- wide_format[,-1]
rownames(wide_format) <- rna_rownames

# View the results
print(wide_format)
rna_agg <- as.matrix(wide_format)
# Use conversion table to align RNA marker names with CODEX protein names
rownames(rna_agg) <- ifelse(rownames(rna_agg) %in% names(names_conversion),
                            names_conversion[rownames(rna_agg)], 
                            rownames(rna_agg))

# Aggregate CODEX protein expression data by CODEX cell type
protein_agg <- AverageExpression(CODEX, slot = 'data', assays = "CODEX", group.by = "cluster_anno_l2", return.seurat = FALSE)$CODEX

# Extracting the data slot which contains normalized data
data_matrix <- CODEX@assays$CODEX@data

# Extract cluster information from metadata
CODEX <- SetIdent(CODEX, value = "cluster_anno_l2")
clusters <- Idents(CODEX)

subset_data <- data_matrix

# Create a data frame for calculation
df <- data.frame(expression = as.vector(subset_data),
                 protein = rep(rownames(subset_data), ncol(subset_data)),
                 cluster = rep(clusters, each = nrow(subset_data)))

# Calculate median expression by gene and cluster
median_expression <- df %>%
  group_by(protein, cluster) %>%
  summarize(median_expr = median(expression, na.rm = TRUE))

# View the results
print(median_expression)

wide_format <- median_expression %>%
  pivot_wider(names_from = cluster, values_from = median_expr) 
codex_rownames <- wide_format$protein
wide_format <- wide_format[,-1]
rownames(wide_format) <- codex_rownames

# View the results
print(wide_format)
protein_agg <- as.matrix(wide_format)

# correlations 

# match the RNA agg and protein agg tables
rna_agg <- as.data.frame(rna_agg)
rna_agg$Marker_Name <- rownames(rna_agg)
protein_agg <- as.data.frame(protein_agg)
protein_agg$Marker_Name <- rownames(protein_agg)

# match genes (note HLA-DR and CD45RA lost due to multiple genes and isoform issue)
rna_agg <- subset(rna_agg, invert = TRUE, rna_agg$Marker_Name %in% protein_agg$Marker_Name)
protein_agg <- subset(protein_agg, invert = TRUE, protein_agg$Marker_Name %in% rna_agg$Marker_Name)

# rename labels
scrna_to_harmonized <- setNames(integration_table$Harmonized_Labels_scRNA, integration_table$scRNA_cla2)
codex_to_harmonized <- setNames(integration_table$Harmonized_Labels_CODEX, integration_table$CODEX_cla2)

scrna_label_translation <- setNames(integration_table$Harmonized_Labels_scRNA, integration_table$scRNA_cla2)
codex_label_translation <- setNames(integration_table$Harmonized_Labels_CODEX, integration_table$CODEX_cla2)


# Rename the columns in the rna_agg matrix
rna_agg_renamed <- rna_agg
colnames(rna_agg_renamed) <- scrna_to_harmonized[colnames(rna_agg_renamed)]
# Check for any NAs and handle them
#na_cols_rna <- is.na(colnames(rna_agg_renamed))
#if (any(na_cols_rna)) {
#  warning("NA values found in rna_agg column names, replacing with original names.")
#  colnames(rna_agg)[na_cols_rna] <- names(na_cols_rna)[na_cols_rna]
#}

# Rename the columns in the protein_agg matrix
protein_agg_renamed <- protein_agg
colnames(protein_agg_renamed) <- codex_to_harmonized[colnames(protein_agg_renamed)]
# Check for any NAs and handle them
#na_cols_protein <- is.na(colnames(protein_agg_renamed))
#if (any(na_cols_protein)) {
#  warning("NA values found in protein_agg column names, replacing with original names.")
#  colnames(protein_agg)[na_cols_protein] <- names(na_cols_protein)[na_cols_protein]
#}


# This function will average columns in a dataframe with the same name
average_duplicate_columns <- function(df) {
  # Find duplicate column names and average them
  duplicated_cols <- which(duplicated(colnames(df)))
  unique_duplicated_cols <- unique(colnames(df)[duplicated_cols])
  
  for (col in unique_duplicated_cols) {
    # Get all columns that match the duplicated name
    cols_to_average <- which(colnames(df) == col)
    
    # Calculate the mean of these columns
    df[[col]] <- rowMeans(df[, cols_to_average, drop = FALSE], na.rm = TRUE)
    
    # Now select the first occurrence and drop the rest
    df <- df[, !duplicated(colnames(df)) | (1:ncol(df)) %in% cols_to_average[1], drop = FALSE]
  }
  
  return(df)
}

# Apply the function to average columns with duplicate names in both dataframes
rna_agg_renamed <- average_duplicate_columns(rna_agg_renamed)
protein_agg_renamed <- average_duplicate_columns(protein_agg_renamed)

# Add "Marker_Names" column as the first column in both dataframes
rna_agg_renamed <- cbind(Marker_Names = rownames(rna_agg_renamed), rna_agg_renamed)
protein_agg_renamed <- cbind(Marker_Names = rownames(protein_agg_renamed), protein_agg_renamed)

# Verify the structure and column names
print(head(rna_agg_renamed))
print(head(protein_agg_renamed))

# Continue with the analysis using the renamed and averaged matrices

# Add "Marker_Names" column as the last column in both dataframes
rna_agg_renamed$Marker_Names <- rownames(rna_agg_renamed)
protein_agg_renamed$Marker_Names <- rownames(protein_agg_renamed)

# Move "Marker_Names" to the first column
rna_agg_renamed <- rna_agg_renamed[, c(ncol(rna_agg_renamed), 1:(ncol(rna_agg_renamed)-1))]
protein_agg_renamed <- protein_agg_renamed[, c(ncol(protein_agg_renamed), 1:(ncol(protein_agg_renamed)-1))]

library(ggplot2)
library(dplyr)

# Assuming rna_agg_renamed and protein_agg_renamed have been created and have 'Marker_Names' as the first column


library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming rna_agg_renamed and protein_agg_renamed have been created with 'Marker_Names' as the first column

# Extract common cell types
common_colnames <- intersect(colnames(rna_agg_renamed)[-1], colnames(protein_agg_renamed)[-1])

rna_agg_renamed <- rna_agg_renamed[,-1]
protein_agg_renamed <- protein_agg_renamed[,-1]


# Pivot the data frames longer
rna_long <- pivot_longer(
  rna_agg_renamed,
  cols = -Marker_Names,  # Exclude the marker names column
  names_to = "Cell_Type",
  values_to = "RNA_Expression"
)

codex_long <- pivot_longer(
  protein_agg_renamed,
  cols = -Marker_Names,  # Exclude the marker names column
  names_to = "Cell_Type",
  values_to = "CODEX_Expression"
)

# Join the two long data frames on Marker_Name and Cell_Type
combined_long <- left_join(rna_long, codex_long, by = c("Marker_Names", "Cell_Type"))

filtered_combined_long <- combined_long %>%
  filter(
    !is.na(RNA_Expression) & RNA_Expression != "No Match" &
      !is.na(CODEX_Expression) & CODEX_Expression != "No Match"
  )

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Assuming filtered_combined_long is already defined and contains the columns:
# 'Marker_Names', 'Cell_Type', 'RNA_Expression', and 'CODEX_Expression'

# Create a PDF file to save the plots
pdf("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Integration_BatchCorrection/median_celltype_correlation_plots.pdf", width = 11, height = 8.5)

# Loop through each unique cell type and create a plot
for(cell_type in unique(filtered_combined_long$Cell_Type)) {
  
  # Subset the data for the current cell type
  current_data <- filtered_combined_long %>% 
    filter(Cell_Type == cell_type) %>%
    na.omit() # Remove rows with NA values
  
  # Calculate correlation coefficient and linear model
  correlation_test <- cor.test(current_data$RNA_Expression, current_data$CODEX_Expression)
  lm_model <- lm(CODEX_Expression ~ RNA_Expression, data = current_data)
  
  # Plot the data and add a regression line
  p <- ggplot(current_data, aes(x = RNA_Expression, y = CODEX_Expression)) +
    geom_point() +
    geom_smooth(method = "lm", color = "blue") +
    geom_text(aes(label = Marker_Names), hjust = 1.5, vjust = 1.5, check_overlap = TRUE) +
    ggtitle(paste("Cell Type:", cell_type, "\nCorrelation R:", round(correlation_test$estimate, 3))) +
    xlab("RNA Expression (normalized)") +
    ylab("Protein Expression (CLR)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Print the plot to the PDF
  print(p)
}

# Close the PDF device
dev.off()


# repeat for scaled

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)

# Assuming filtered_combined_long is already defined and contains the columns:
# 'Marker_Names', 'Cell_Type', 'RNA_Expression', and 'CODEX_Expression'

# Function to scale data to a -1 to 1 range
scale_to_range <- function(x) {
  (2 * (x - min(x)) / (max(x) - min(x))) - 1
}

# Apply scaling within each cell type for each assay
filtered_combined_long <- filtered_combined_long %>%
  group_by(Cell_Type) %>%
  mutate(
    RNA_Expression_Scaled = scale_to_range(RNA_Expression),
    CODEX_Expression_Scaled = scale_to_range(CODEX_Expression)
  ) %>%
  ungroup()


filtered_combined_long %>% group_by(Cell_Type) %>% summarise((range(RNA_Expression_Scaled)))
filtered_combined_long %>% group_by(Cell_Type) %>% summarise(range(CODEX_Expression_Scaled))

# Open PDF to save plots
pdf("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Integration_BatchCorrection/scaled_median_celltype_correlation_plots.pdf")

# Loop through each unique cell type and create a plot
for(cell_type in unique(filtered_combined_long$Cell_Type)) {
  
  # Subset the data for the current cell type
  current_data <- filtered_combined_long %>% 
    filter(Cell_Type == cell_type) %>%
    na.omit() # Remove rows with NA values
  
  # Calculate correlation coefficient and linear model
  correlation_test <- cor.test(current_data$RNA_Expression_Scaled, current_data$CODEX_Expression_Scaled)
  lm_model <- lm(CODEX_Expression_Scaled ~ RNA_Expression_Scaled, data = current_data)
  
  # Plot the data and add a regression line
  p <- ggplot(current_data, aes(x = RNA_Expression_Scaled, y = CODEX_Expression_Scaled)) +
    geom_point() +
    geom_smooth(method = "lm", color = "blue") +
    geom_text(aes(label = Marker_Names), hjust = 1.5, vjust = 1.5, check_overlap = TRUE) +
    ggtitle(paste("Cell Type:", cell_type, "\nCorrelation R:", round(correlation_test$estimate, 3))) +
    xlab("Scaled RNA Expression (normalized)") +
    ylab("Scaled Protein Expression (CLR)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Print the plot to the PDF
  print(p)
}

# Close the PDF device
dev.off()




# correlation approach ----

library(dplyr)
library(tidyr)

library(dplyr)
library(tidyr)

# Assuming filtered_combined_long is already defined and contains columns for 
# 'Marker_Names', 'Cell_Type', 'RNA_Expression', and 'CODEX_Expression'

library(dplyr)
library(tidyr)

# Assuming filtered_combined_long is already defined and contains the columns:
# 'Marker_Names', 'Cell_Type', 'RNA_Expression', and 'CODEX_Expression'

# Reshape data into separate wide formats for RNA and CODEX
rna_wide <- filtered_combined_long %>%
  select(Marker_Names, Cell_Type, RNA_Expression) %>%
  pivot_wider(names_from = Cell_Type, values_from = RNA_Expression) %>%
  column_to_rownames(var = "Marker_Names")

codex_wide <- filtered_combined_long %>%
  select(Marker_Names, Cell_Type, CODEX_Expression) %>%
  pivot_wider(names_from = Cell_Type, values_from = CODEX_Expression) %>%
  column_to_rownames(var = "Marker_Names")

# Ensure the marker names are in the same order for RNA and CODEX
common_markers <- intersect(rownames(rna_wide), rownames(codex_wide))
rna_wide <- rna_wide[common_markers, , drop = FALSE]
codex_wide <- codex_wide[common_markers, , drop = FALSE]

# Calculate the correlation matrix
correlation_matrix <- cor(rna_wide, codex_wide, use = "pairwise.complete.obs",method = "pearson")

# Determine the closest scRNA-Seq cell type for each CODEX cell type
closest_matches <- apply(correlation_matrix, 2, function(x) {
  names(which.max(x))
})

# Create a dataframe for the closest matches
closest_matches_df <- data.frame(
  CODEX_Cell_Type = colnames(correlation_matrix),
  Closest_scRNA_Seq_Cell_Type = closest_matches,
  stringsAsFactors = FALSE
)

# Print the closest matches
print(closest_matches_df)

library(gt)
library(espnscrapeR)
gt_table <- gt(closest_matches_df)
closest_matches_df %>% dplyr::filter(CODEX_Cell_Type != "No Match") %>% 
gt() %>% 
  gt_theme_538() %>% 
  tab_header(title = md("**CODEX/scRNASeq Cell Type Correlation**")) %>% gtsave(filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Integration_BatchCorrection/closest_cell_type_table.pdf")
```{r}
closest_matches_df %>% dplyr::filter(CODEX_Cell_Type != "No Match") %>% 
gt() %>% 
  gt_theme_538() %>% 
  tab_header(title = md("**CODEX/scRNASeq Cell Type Correlation**"))
```



# make a plot of the correlation approach ----

library(ggplot2)
library(dplyr)

# Original color scheme
scrna_cal2_cols <- c("#FFD580", "#FFBF00", "#B0DFE6", "#7DB954", "#64864A", "#8FCACA", "#4682B4", "#CAA7DD", "#B6D0E2", 
                     "#a15891", "#FED7C3", "#A8A2D2", "#CF9FFF", "#9C58A1", "#2874A6", "#96C5D7", "#63BA97", "#BF40BF", 
                     "#953553", "#6495ED", "#E7C7DC", "#5599C8", "#FA8072", "#F3B0C3", "#F89880", "#40B5AD", "#019477", 
                     "#97C1A9", "#C6DBDA", "#CCE2CB", "#79127F", "#ff9d5c", "#FFC8A2", "#DD3F4E")
cal2_col_names <- c("Adipo-MSC", "AEC", "Ba/Eo/Ma", "CD4+ T-Cell", "CD8+ T-Cell", "CLP", "Cycling DCs", "Cycling HSPC", 
                    "Early Myeloid Progenitor", "Erythroblast", "Fibro-MSC", "GMP", "HSC", "Late Erythroid", "Late Myeloid", 
                    "Macrophages", "Mature B", "Megakaryocyte", "MEP", "Monocyte", "MPP", "Neutrophil", "Osteo-MSC", 
                    "Osteoblast", "OsteoFibro-MSC", "pDC", "Plasma Cell", "Pre-B", "Pre-Pro B", "Pro-B", "RBC", "SEC", 
                    "THY1+ MSC", "VSMC")

# New unique labels
new_labels <- c("Adipo-MSC", "AEC", "Autofluorescent", "CD4+ T-Cell", "CD8+ T-Cell", "CLP", "pDC", "HSPC", 
                "Early Myeloid Progenitor", "Erythroblast", "No Match", "GMP", "Early HSPC", "Erythroid", 
                "Intermediate Myeloid", "Macrophages", "B_Cells", "Megakaryocyte", "MEP/Early Erythroblast", "Monocytes", 
                "Mature Myeloid", "Endosteal", "Plasma Cells", "Immature_B_Cell", "SEC", "THY1+ MSC", "VSMC")

# Initialize a vector to store colors
new_label_colors_vector <- rep(NA, length(new_labels))
names(new_label_colors_vector) <- new_labels

# Function to match colors or assign a new one
get_color <- function(label, original_names, original_colors) {
  if(label %in% original_names) {
    return(original_colors[original_names == label])
  } else {
    # Define custom colors for labels not in the original set
    custom_colors <- c(
      "Autofluorescent" = "#D3D3D3",  # Light grey for "Autofluorescent"
      "GMP" = "#CF9FFF",   # Purple for "GMP/Myeloblast"
      "Early HSPC" = "#CAA7DD",       # Light purple for "Early HSPC"
      "Erythroid" = "#ff9d5c",        # Coral for "Erythroid"
      "Intermediate Myeloid" = "#FFC8A2", # Light orange for "Intermediate Myeloid"
      "MEP/Early Erythroblast" = "#a15891", # Dark pink for "MEP/Early Erythroblast"
      "Mature Myeloid" = "#FFBF00",   # Orange for "Mature Myeloid"
      "Endosteal" = "#FA8072",        # Salmon for "Endosteal"
      "Immature_B_Cell" = "#7DB954",  # Green for "Immature_B_Cell"
      "No Match" = "#000000"          # Black for "No Match"
    )
    return(ifelse(label %in% names(custom_colors), custom_colors[label], "#FFFFFF")) # White as default
  }
}

# Assign colors to new labels
for(i in seq_along(new_labels)) {
  label <- new_labels[i]
  new_label_colors_vector[i] <- get_color(label, cal2_col_names, scrna_cal2_cols)
}

# Print the named vector to confirm
new_label_colors_vector["GMP"] <- "#CF9FFF"

# fix the messed up ones
new_label_colors_vector["Plasma Cells"] <- '#228B22'
new_label_colors_vector["B_Cells"] <- '#4F9153'
new_label_colors_vector["Monocytes"] <- '#6495ED'
  
  



# Assuming you have already calculated the correlation matrix and have the closest_matches_df from previous steps
# and 'scrna_cal2_cols' is your color scheme for scRNA-Seq cell types

# Create a long dataframe for plotting, with a row for each correlation and the corresponding RNA cell type
plot_data <- as.data.frame(correlation_matrix) %>%
  tibble::rownames_to_column("RNA_Cell_Type") %>%
  pivot_longer(cols = -RNA_Cell_Type, names_to = "CODEX_Cell_Type", values_to = "Correlation") %>%
  mutate(RNA_Cell_Type = factor(RNA_Cell_Type, levels = names(new_label_colors_vector)))

library(ggplot2)
library(dplyr)

# Assuming you have a data frame 'data' from the previous steps with columns 'CODEX_Cell_Type', 
# 'RNA_Cell_Type', 'Correlation', and that 'highest_correlation_df' is a data frame with columns 
# 'CODEX_Cell_Type' and 'Highest_RNA_Cell_Type' indicating the scRNA-Seq cell type with the highest correlation for each CODEX cell type.

library(ggplot2)
library(dplyr)
library(ggrepel)

# Assuming 'data' is your dataframe and 'highest_correlation_df' is defined as mentioned before

library(ggplot2)
library(dplyr)
library(ggrepel)

# Assuming plot_data has been created with columns: 'CODEX_Cell_Type', 'RNA_Cell_Type', 'Correlation'
# and that 'new_label_colors_vector' is your named vector of colors for the RNA cell types

# Add a flag for the points where the RNA label matches the CODEX label
plot_data <- plot_data %>%
  mutate(Is_Matching = RNA_Cell_Type == CODEX_Cell_Type)

# Now create the plot
p <- ggplot(plot_data, aes(x = CODEX_Cell_Type, y = Correlation, color = RNA_Cell_Type)) +
  geom_point(aes(shape = Is_Matching), size = 3, stroke = 1.5) +
  scale_shape_manual(values = c(19, 21)) +  # Different shapes for matching/non-matching points
  scale_color_manual(values = new_label_colors_vector) +
  geom_text_repel(
    data = filter(plot_data, Is_Matching),
    aes(label = RNA_Cell_Type),
    nudge_y = 0.05,
    nudge_x = 0.25,  # Adjust this value as needed to prevent text overlap
    segment.color = 'grey50'
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") +
  labs(title = "Correlation between CODEX and scRNA-Seq Cell Types",
       x = "CODEX Cell Type", y = "Correlation") +
  guides(shape = FALSE) & RotatedAxis() # Hide the shape legend 

# Print the plot
print(p)




# Save the plot as a file
ggsave("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Integration_BatchCorrection/CODEX_scRNA_Seq_Correlation_Plot_Adjusted_medians.pdf", p, width = 16, height = 8)
