library(spatstat)
library(png)
library(stringr)
library(foreach)
library(doParallel)
library(parallel)
library(Seurat)
library(ggpubr)

setwd("Seurat_Analysis/Bone_Marrow")
immune.combined <- readRDS("Immune_Combined_Added_MSC_CD56.RDS")
pattern <- "H\\d{2}"
sample_fixed <- str_extract(immune.combined$orig.ident, pattern)
immune.combined$sample_fixed <- sample_fixed

table(immune.combined$cluster_anno_l2)
DefaultAssay(immune.combined) <- "CODEX2"
sample_names <- c('H10', 'H14', 'H26', 'H27', 'H32', 'H33', 'H35', 'H36', 'H37', 'H38', 'H39', 'H41')
filenames <- c('H10.csv','H14.csv','H26.csv', 'H27.csv', 'H32.csv', 'H33.csv', 'H35.csv', 
               'H36.csv', 'H37.csv', 'H38.csv', 'H39.csv', 'H41.csv')

unique_cell_types <- c("Osteo_MSC", "Fibro_MSC", "THY1+ MSC", "MSC")
osteomsc_dist <- c()
fibromsc_dist <- c()
thy1msc_dist <- c()
adipomsc_dist <- c()

structure = 'bone'
for (i in 1:length(sample_names)) {
  sample_name <- sample_names[i]
  print(sample_name)
  filename <- filenames[i]
  sample_data <- subset(immune.combined, subset = sample_fixed == sample_name)
  
  cts <- read.csv(paste0('Bone_Contours/', filename))
  cells_on_bone <- read.table("Osteo_Fibro_MSC_on_Bone.txt", header = TRUE, sep = '\t')
  cells_on_bone <- cells_on_bone[cells_on_bone$sample == sample_name,]
  cells_on_bone$sample <- NULL
  cts <- cbind("bone", cts)
  colnames(cts)[1] = "celltype"
  cts <- rbind(cts, cells_on_bone)
  
  xmin <- min(cts$x0)
  xmax <- max(cts$x1)
  ymin <- min(cts$y0)
  ymax <- max(cts$y1)
  
  sample_data_crop_xy <- sample_data@meta.data[, c('x.coord', 'y.coord')]
  
  for (j in 1:length(unique_cell_types)) {
    cell_type <- unique_cell_types[j]
    print(cell_type)

    #Add manual annotations of osteo and fibro MSC
    coords_crop <- subset(sample_data, subset = cluster_anno_l2 == cell_type)
    x <- coords_crop$x.coord
    y <- coords_crop$y.coord

    population <- length(x)
    print(population)

    cells <- ppp(x, y, c(xmin-10, xmax+10), c(ymin-10, ymax+10))
    bone_cts <- psp(cts$x0, cts$y0, cts$x1, cts$y1, window = cells$window)
    
    bm <-list(cell=cells,bone=bone_cts)
    bm$dbone <- distfun(bm$bone)
       
    distance <- bm$dbone(x,y)
    # remove NA values
    distance <- distance[!is.na(distance)]
    # remove infinite values
    distance <- distance[!is.infinite(distance)]
    
    #Append to lists 
    unique_cell_types <- c("Osteo_MSC", "Fibro_MSC", "THY1+ MSC", "MSC")
    
    if(cell_type == "Osteo_MSC"){
      osteomsc_dist <- c(osteomsc_dist, distance)
    }
    if(cell_type == "Fibro_MSC"){
      fibromsc_dist <- c(fibromsc_dist, distance)
    }
    if(cell_type == "THY1+ MSC"){
      thy1msc_dist <- c(thy1msc_dist, distance)
    }
    if(cell_type == "MSC"){
      adipomsc_dist <- c(adipomsc_dist, distance)
    }
  }
}

#Figure 4E Make violin plot ----
msc_data <- data.frame(
  Group = rep(c("Osteo_MSC", "Fibro_MSC", "THY1+ MSC", "Adipo-MSC"), 
              times = c(length(osteomsc_dist), length(fibromsc_dist), length(thy1msc_dist), length(adipomsc_dist))),
  Value = c(osteomsc_dist, fibromsc_dist, thy1msc_dist, adipomsc_dist)
)
msc_data$Value = msc_data$Value * 0.5069
msc_data$Group <- factor(msc_data$Group, levels = c("Fibro_MSC", "Osteo_MSC", "THY1+ MSC", "Adipo-MSC"))
saveRDS(msc_data, "MSC_Distance_to_Bone.RDS")

ggplot(msc_data, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(width = 1.2) +
  geom_boxplot(width = 0.18, outlier.shape = NA, position = position_dodge(0.75), fill = "white") +
  theme_bw() + RotatedAxis() +
  labs(y = "Distance to Bone (um)") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 300)) +
  stat_compare_means(comparisons = list(c("Fibro_MSC", "Osteo_MSC"), c("Fibro_MSC", "THY1+ MSC"), c("Fibro_MSC", "Adipo-MSC"),
                                        c("Osteo_MSC", "THY1+ MSC"), c("Osteo_MSC", "Adipo-MSC"), c("THY1+ MSC", "Adipo-MSC")),
                     method = "wilcox.test")