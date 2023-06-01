### NBM38_H33 Seurat Object Generation Script
### This script takes as input a channel_names text file and the cell X protein output of our Mesmer pipeline and outputs a Seurat Object for NBM38

library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)

# read in channelnames
channel_names <- read_csv(file = "/mnt/isilon/tan_lab_imaging/FUSION/NBM27_H10_CITRATE_REIMAGE/H10/Scan1/.temp/MarkerList.txt", col_names = FALSE)

#Specify the directory of your data
file_path = '/mnt/isilon/tan_lab_imaging/FUSION/NBM38_H33_citrate/H33/Scan1/mesmercombined_markers.csv'
#Read cell segmentation data using the specified file name, e.g., 'tonsil_cell_seg_data.txt'
#gc_csd_raw<-read_cell_seg_data(file.path(base_path, 'arlene_cell_seg_data.txt'))
gc_csd_raw <- read_csv(file.path(file_path))

# gc_csd_raw<- dplyr::filter(gc_csd_raw, gc_csd_raw$`Ch1Cy1 vs Frequency 4 (reg 1)` == 1)
#Remove any Phenotypes from the data. This does not apply to CODEX data.
# gc_csd_raw<-gc_csd_raw%>%select(-contains("Phenotype"), -contains("UMAP"))
colnames(gc_csd_raw)[5:58] <- channel_names$X1
colnames(gc_csd_raw)[2] <- "Size"
colnames(gc_csd_raw)[1] <- "CellID"
colnames(gc_csd_raw)[3:4] <- c("x.coord", "y.coord")
gc_csd_raw[,3:58] <- gc_csd_raw[,3:58] / gc_csd_raw$Size 

# Check QC for DAPI 
gc_csd_raw %>% ggplot(aes(x=`DAPI`)) + geom_histogram(bins = 100) + xlim(c(0,255))
gc_csd_raw %>% ggplot(aes(x=`x.coord`, y = -`y.coord`)) + geom_point(alpha = 0.1)
gc_csd_raw <- gc_csd_raw %>% dplyr::filter(`DAPI`>10, `DAPI`<250)

# Crop Image if Needed
gc_csd_raw <- gc_csd_raw %>% dplyr::filter(x.coord < 12500, y.coord >3000, x.coord > 7458, y.coord < 11106)
cells_to_remove <- subset(gc_csd_raw, subset = x.coord > 10074 & y.coord < 9440 & x.coord < 10870 & y.coord > 8064) # Use this line to further crop out piece which had artifact
gc_csd_raw <- anti_join(gc_csd_raw, cells_to_remove, by = "CellID")

gc_csd_raw %>% ggplot(aes(x=`x.coord`, y = -`y.coord`)) + geom_point(alpha = 0.1) # Visually confirm crop

# Specify markers to be excluded from clustering, NOTE that if the name of one marker is contained in another this will remove both!
failed <- c("DAPI", "NAKATPASE", "Ki67", "ADIPOQ", "CD56")

#Specify the cell compartment and markers for clustering. Remove DAPI channel
# remove those cells which are positive for everything and negative for everything
gc_csd_raw <- tibble(gc_csd_raw) %>% rowwise() %>%
  dplyr::mutate(RowSum = sum(c_across(cols = ("CD19":"CXCR4"))))
gc_csd_raw %>% ggplot(aes(x=RowSum)) + geom_histogram(bins = 200) + geom_vline(xintercept = quantile(gc_csd_raw$RowSum, probs = (0.99)))
gc_csd_raw <- gc_csd_raw %>% dplyr::filter(RowSum < quantile(gc_csd_raw$RowSum, probs = (0.99)), RowSum > quantile(gc_csd_raw$RowSum, probs = (0.01))) # remove outlier cells with 99th percentile of row sums


# Remove the size and x/y coords
gc_csd2 <- gc_csd_raw %>% select(-starts_with('x.coord'), -contains("Empty"), -starts_with('y.coord'),-contains('size'),-contains("Std.Dev"), -contains("indexes"),-contains('DAPI'), -contains('RowSum'), -contains('Blank'),-starts_with("Ch"), -starts_with("faile"), -contains(failed)) # Note removed failed markers
gc_csd2 <- gc_csd2 %>% select(-contains("Index"),-contains("CellID"), -matches("Python_Index"))
#Rename objects
gc_csd2 = gc_csd2 %>%
  rename_with(~str_remove_all(.x, 'Nucleus Intensity')) %>%
  rename_with(~str_remove_all(.x, '\\(.*\\) ')) %>%
  rename_with(~str_remove_all(.x, ' Mean'))
gc_csd2<-na.omit(gc_csd2)

# Create a Seurat object
t_gc_csd = t(as.matrix(gc_csd2))
colnames(t_gc_csd) <- rownames(gc_csd2)
rownames(t_gc_csd) <- colnames(gc_csd2)
codex_H33 <- CreateSeuratObject(t_gc_csd, project = "SB67_NBM38_H33_CODEX_Mesmer", assay="CODEX", meta.data = gc_csd_raw)
saveRDS(codex_H33, "codex_H33_SeuratObj.RDS")
