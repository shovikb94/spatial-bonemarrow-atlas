
library(ggrepel)
library(ggplot2)
library(CellChat)
library(dplyr)
library(readr)
library(viridis)

setwd("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/cellchat_integration")

cal2_cols <- c("#FFD580", "#FFBF00", "#B0DFE6", "#7DB954", "#64864A", "#8FCACA", "#4682B4", "#CAA7DD", "#B6D0E2", "#a15891", "#FED7C3", "#A8A2D2", "#CF9FFF", "#9C58A1", "#2874A6", "#96C5D7", "#63BA97", "#BF40BF", "#953553", "#6495ED", "#E7C7DC", "#5599C8", "#FA8072", "#F3B0C3", "#F89880", "#40B5AD", "#019477", "#97C1A9", "#C6DBDA", "#CCE2CB", "#79127F", "#FFC5BF", "#ff9d5c", "#FFC8A2", "#DD3F4E")
cal2_col_names <- c("Adipo-MSC", "AEC", "Ba/Eo/Ma", "CD4+ T-Cell", "CD8+ T-Cell", "CLP", "Cycling DCs", "Cycling HSPC", "Early Myeloid Progenitor", "Erythroblast", "Fibro-MSC", "GMP", "HSC", "Late Erythroid", "Late Myeloid", "Macrophages", "Mature B", "Megakaryocyte", "MEP", "Monocyte", "MPP", "Neutrophil", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "pDC", "Plasma Cell", "Pre-B", "Pre-Pro B", "Pro-B", "RBC", "RNAlo MSC", "SEC", "THY1+ MSC", "VSMC")
names(cal2_cols) <- cal2_col_names

# Example data - replace this with your actual data
set.seed(0)
neighborhood_fc <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/output/neighborhood_fc.csv")
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

# Set row and column names for the neighborhood_fc matrix
neighborhood_fc <- neighborhood_fc[,2:33]
rownames(neighborhood_fc) <- make.unique(neighborhood_names_zeroindex, sep = "_")

neighborhoods <- make.unique(neighborhood_names_zeroindex, sep = "_")
cell_types <- colnames(neighborhood_fc)


# Initialize a list to store neighborhood-specific binary columns
binary_columns <- list()

# Specify the fold change threshold
fold_change_threshold <- 1.0

# Create a function to check if a cell type pair is enriched
is_enriched <- function(fc_values) {
  return(all(fc_values >= fold_change_threshold))
}

# Loop through each neighborhood
for (neighborhood_to_analyze in neighborhoods) {
  
  # Create a logical vector to store binary values for cell type pairs
  is_enriched_vector <- logical(0)
  cp_sum_vector <- numeric(0)
  cell_type_pairs <- character(0)
  
  # Loop through each cell type pair
  for (i in 1:(length(cell_types) - 1)) {
    for (j in (i + 1):length(cell_types)) {
      cell_type_i <- cell_types[i]
      cell_type_j <- cell_types[j]
      
      # Get the fold change values for the specific cell type pair in the current neighborhood
      fc_values <- neighborhood_fc[neighborhood_to_analyze, c(cell_type_i, cell_type_j)]
      
      # Check if the cell type pair is enriched
      is_enriched_value <- is_enriched(fc_values)
      
      # calculate the sum
      cp_sum_value <- sum(neighborhood_fc[neighborhood_to_analyze, c(cell_type_i, cell_type_j)])
      
      # Append the result to the vectors
      is_enriched_vector <- c(is_enriched_vector, is_enriched_value)
      cp_sum_vector <- c(cp_sum_vector, cp_sum_value)
      cell_type_pairs <- c(cell_type_pairs, paste(cell_type_i, cell_type_j, sep = "-"))
    }
  }
  
  # Store the binary vector for the current neighborhood
  binary_columns[[neighborhood_to_analyze]] <- cp_sum_vector
}

# Create the result dataframe with binary columns for each neighborhood
result_df <- data.frame(cell_type_pair = cell_type_pairs)
for (neighborhood_to_analyze in neighborhoods) {
  binary_column_name <- paste0(neighborhood_to_analyze)
  result_df[binary_column_name] <- binary_columns[[neighborhood_to_analyze]]
}

# scale result by columns
result_df_scaled <- as.data.frame(scale(result_df[2:16]))
rownames(result_df_scaled) <- result_df$cell_type_pair

result_df_scaled$CellTypePair <- rownames(result_df_scaled)

result_df_scaled$CellTypePair <- gsub(x = result_df_scaled$CellTypePair, "B-Cell", replacement = "B_Cell")
result_df_scaled$CellTypePair <- gsub(x = result_df_scaled$CellTypePair, "Non-Classical Monocyte", replacement = "Non_Classical Monocyte")
result_df_scaled$CellTypePair <- gsub(x = result_df_scaled$CellTypePair, "T-Cell", replacement = "T_Cell")



# Split the "CellTypePair" column into two new columns
split_names <- strsplit(result_df_scaled$CellTypePair, "-")

# Create "Cell Type 1" and "Cell Type 2" columns
result_df_scaled$CellType1 <- sapply(split_names, function(x) x[1])
result_df_scaled$CellType2 <- sapply(split_names, function(x) x[2])

#save
write.csv(result_df_scaled, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/output/scaled_neighborhood_colocalization_scores.csv")

# read integration table
integration <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/cellchat_integration/CODEX_scRNA_labels_integration.csv")

cellchat <- read_csv("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/CellChat/Figures/cellchat_LRpairs_updated.csv")


library(dplyr)

# Join the dataframes
# Join the dataframes
result1 <- cellchat %>%
  # Join with df3 on source to get CellType1
  left_join(integration, by = c("source" = "scRNA_cla2")) %>%
  rename(CellType1 = 'CODEX_closest_cla2' ) %>%
  # Join again with df3 on target to get CellType2
  left_join(integration, by = c("target" = "scRNA_cla2"), suffix = c("", ".target")) %>%
  rename(CellType2 = "CODEX_closest_cla2") 

library(dplyr)
result1$NEWCOL <- paste(pmin(result1$CellType1, result1$CellType2),
                        pmax(result1$CellType1, result1$CellType2),
                        sep = '')
result_df_scaled$NEWCOL <- paste(pmin(result_df_scaled$CellType1, result_df_scaled$CellType2),
                                 pmax(result_df_scaled$CellType1, result_df_scaled$CellType2),
                                 sep = '')
result_df_scaled$maxscore <- apply(result_df_scaled[,1:15], MARGIN = 1, FUN = max)
cellchat_integrated <- left_join(result1, result_df_scaled, by = 'NEWCOL')
cellchat_integrated <- cellchat_integrated %>% select(-c("CellType1.y", "CellType2.y"))
colnames(cellchat_integrated)[13] <- "CODEX_CellType1"
colnames(cellchat_integrated)[14] <- "CODEX_CellType2"

# Write the data to a new CSV file 
# write.csv(data, "/path/to/your/cellchat_codex_integrated_with_counts.csv", row.names = FALSE)

# Repeat the cellchat / neighborhood correlation for the neighborhood matrix ----
cellchat <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/CellChat/cellchat_032923.RDS")
weights <- as.data.frame(cellchat@net$weight)

data <- cellchat_integrated

# Melt the matrix into a long format
weights_long <- as.data.frame(as.table(as.matrix(weights)))
names(weights_long) <- c("source", "target", "weight")

# Now merge this weights_long dataframe with your existing dataframe
# Assuming 'data' is your existing dataframe
# Average maxscore for each source-target pair
data_avg <- data %>%
  group_by(source, target) %>%
  summarise(average_maxscore = mean(maxscore, na.rm = TRUE)) %>%
  ungroup()

# Merge averaged maxscores back to the original data
data <- data %>%
  left_join(data_avg, by = c("source", "target"))


data <- left_join(data, weights_long, by = c("source", "target"))


# Identify points of interest and create a new column for coloring
data$point_of_interest <- data$average_maxscore > 2 & data$weight > 1

# extract dataframe for just columns of interest
data <- unique(data.frame(data$source, data$target, data$weight, data$average_maxscore, data$point_of_interest))
vec_names <- c("source", "target", "weight", "average_maxscore", "point_of_interest")
colnames(data) <- vec_names
data <- data %>% filter(!is.nan(average_maxscore))
model <- lm(average_maxscore ~ weight, data = data)


# Calculate correlation statistics with averaged maxscore
correlation_test <- cor.test(data$weight, data$average_maxscore, method = "pearson")

# Extract the R and p-value from the test
R_value <- correlation_test$estimate
p_value <- correlation_test$p.value


# Check how many points of interest we have
points_of_interest <- filter(data, average_maxscore > 2 & data$weight > 1)
print(paste("Number of points of interest:", nrow(points_of_interest)))

# If we have points of interest, proceed with the plot
if (nrow(points_of_interest) > 0) {
  
  # Create the scatter plot with a correlation line and labeled points of interest
  p <- data %>% ggplot(aes(x = weight, y = average_maxscore)) +
    geom_point(aes(color = point_of_interest)) + # Color points based on the point of interest
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) + # Assign colors
    #vgeom_smooth(method = "lm", se = FALSE, color = "blue") +
    geom_text_repel(
      data = points_of_interest,
      aes(label = paste0(source,"-",target)), 
      size = 4, 
      color = "red") +
    labs(x = "CellChat Weight", y = "Average Max Score", 
         title = "Correlation Plot",
         subtitle = paste("Pearson R:", round(R_value, 3), "p-value:", format.pval(p_value))) +
    theme_minimal() +
    theme(legend.position = "none")
}

# Print the plot
print(p)

# create combined heatmap

# Calculate the unscaled distance from each point to the line y = x
data$distance_to_line = abs(data$weight - data$average_maxscore) / sqrt(2)

# To scale this distance to a correlation-like metric between -1 and 1, we need to normalize it
# For example, you could simply inverse the distance, so that a smaller distance gives a higher score
data$correlation_metric = 1 - data$distance_to_line / max(data$distance_to_line)

# If you want to differentiate between positive and negative correlation (assuming anticorrelation for weight < average_maxscore)
data$correlation_metric = ifelse(data$weight < data$average_maxscore, -data$correlation_metric, data$correlation_metric)

# Now you can rank the cell type pairs based on this correlation metric
data <- data %>% arrange(desc(correlation_metric))

# Ensure that 'weight' and 'average_maxscore' are min/max scaled
data <- data %>%
  mutate(
    scaled_weight = (weight - min(weight)) / (max(weight) - min(weight)),
    scaled_maxscore = (average_maxscore - min(average_maxscore)) / (max(average_maxscore) - min(average_maxscore))
  )

# Calculate the geometric mean of the scaled values
data <- data %>%
  rowwise() %>%
  mutate(effect_size = sqrt(scaled_weight * scaled_maxscore))

# Now you can rank the cell type pairs based on this effect size
data <- data %>% arrange(desc(effect_size))

# Map source and target to their respective colors
data$source_color <- cal2_cols[data$source]
data$target_color <- cal2_cols[data$target]

# Figure 5F - Create the bubble plot ----

p1 <- ggplot(data, aes(x = source, y = target, size = correlation_metric)) +
  geom_point(shape = 21, aes(fill = effect_size)) + # Use shape 21 which has a border and fill
  scale_fill_gradient(name = "Effect Size", ,
                      low = "white", high = "purple") + # Set fill scale to viridis
  geom_point(shape = 21, fill = NA) + # Overlay for source color outline
  geom_point(shape = 21, fill = NA) + # Overlay for target color outline
  theme_minimal() +# Hide legend, optional
  labs(x = "Source Cell Type", y = "Target Cell Type", title = "Bubble Plot of Cell Type Pairs") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(p1, device = "pdf", height = 7, width = 7.5, filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Neighborhood_Analysis_Step4/cellchat_integration/CODEX_scRNA_Integrated_BubblePlot_purples.pdf")


