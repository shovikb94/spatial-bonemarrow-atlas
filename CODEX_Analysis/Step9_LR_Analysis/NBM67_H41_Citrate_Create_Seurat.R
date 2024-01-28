library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
options(Seurat.object.assay.version = "v4")

# read in channelnames
channel_names <- read_csv(file = "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/NBM67/MarkerList.txt", col_names = FALSE)

#Specify the directory of your data
file_path = '/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/NBM67/mesmer/combined_markers.csv'
gc_csd_raw <- read_csv(file.path(file_path))
z = dim(gc_csd_raw)[2]

colnames(gc_csd_raw)[5:z] <- channel_names$X1
colnames(gc_csd_raw)[2] <- "Size"
colnames(gc_csd_raw)[1] <- "CellID"
colnames(gc_csd_raw)[3:4] <- c("x.coord", "y.coord")
gc_csd_raw[,3:z] <- gc_csd_raw[,3:z] / gc_csd_raw$Size 

# Check QC for DAPI 
gc_csd_raw %>% ggplot(aes(x=`DAPI`)) + geom_histogram(bins = 100) + xlim(c(0,255))
gc_csd_raw %>% ggplot(aes(x=`x.coord`, y = -`y.coord`)) + geom_point(alpha = 0.1)
gc_csd_raw <- gc_csd_raw %>% dplyr::filter(`DAPI`>10, `DAPI`<250)

# Crop Image if Needed
gc_csd_raw <- gc_csd_raw %>% dplyr::filter(y.coord < 39000, x.coord < 15000, x.coord > 2500, y.coord > 20000)
gc_csd_raw %>% ggplot(aes(x=`x.coord`, y = -`y.coord`)) + geom_point(alpha = 0.1) # Visually confirm crop

# Specify markers to be excluded from clustering, NOTE that if the name of one marker is contained in another this will remove both!
failed <- c("DAPI", "NAKATPASE")
first = colnames(gc_csd_raw)[6]
last = colnames(gc_csd_raw)[z]
print(paste0(first, ":", last))

#Specify the cell compartment and markers for clustering. Remove DAPI channel
# remove those cells which are positive for everything and negative for everything
gc_csd_raw <- tibble(gc_csd_raw) %>% rowwise() %>%
  dplyr::mutate(RowSum = sum(c_across(cols = (first:last))))
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
codex_assay <- CreateAssayObject(counts = t_gc_csd, project = "NBM67_H41_Citrate_CODEX", assay="CODEX", meta.data = gc_csd_raw)
codex_seurat <- CreateSeuratObject(t_gc_csd, project = "NBM67_H41_Citrate_CODEX", assay="CODEX", meta.data = gc_csd_raw)
codex_seurat[["CODEX"]] = codex_assay
saveRDS(codex_seurat, "/mnt/isilon/tan_lab/sussmanj/CODEX/NBM_Receptor_Ligand/Seurat_Objects/NBM67_H41_Citrate_CODEX_SeuratObj.RDS")
