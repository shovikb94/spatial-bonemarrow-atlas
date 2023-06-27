library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(rstatix)


immune.filtered <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.filtered_FINAL.RDS")
immune.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.combined_final_040523.RDS")

DefaultAssay(immune.filtered) <- "integrated"
DefaultAssay(immune.combined) <- "integrated"

# Give coarse definition for combined object
levels(as.factor(immune.combined$cluster_anno_l2))
coarse_annotations <- c("Adipocyte", "AEC","Artifact","Artifact","Lymphoid", "HSPC", "Lymphoid","Artifact","Lymphoid", "HSPC", "Myeloid", "Endosteal", "Erythroid", "Erythroid", "Megakaryocyte", "Megakaryocyte", 
                        "HSPC", "Myeloid", "HSPC", "HSPC", "Lymphoid", "Myeloid", "Macrophages","Myeloid", "HSPC", "Monocytes", "MSC", "Monocytes", "pDC", "Lymphoid", "Schwann Cells", "SEC", "HSPC", "MSC","Artifact","VSMC")
immune.combined <- SetIdent(immune.combined, value = "cluster_anno_l2")
names(coarse_annotations) <- levels(as.factor(immune.combined$cluster_anno_l2))
immune.combined <- RenameIdents(immune.combined, coarse_annotations)
immune.combined@meta.data["cluster_anno_coarse"] <- immune.combined@active.ident 

# Rename MSC to Adipo-MSC based on scRNA-Seq
# AdipoMSC_ids <- colnames(subset(immune.combined, subset = cluster_anno_l2 == "MSC")) 
# immune.combined$cluster_anno_l2[AdipoMSC_ids] <- "Adipo-MSC"

immune.filtered <- RunUMAP(immune.filtered, dims = 1:30, reduction = "pca", reduction.name = "UMAP_dim30", return.model = TRUE)
immune.combined <- RunUMAP(immune.combined, dims = 1:30, reduction = "pca", reduction.name = "UMAP_dim30", return.model = TRUE)

# saveRDS(immune.combined, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.combined_final_040523.RDS")

immune.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/Annotate_Cells_Step3/objects/immune.combined_final_040523.RDS")

# relevel to be consistent with color scheme 
#immune.combined$cluster_anno_l2 <- factor(immune.combined$cluster_anno_l2, levels = c("Adipocyte", "AEC", "Artifact", "Autofluorescent", "B-Cells", "CD34+ CD61+", "CD4+ T-Cell", "CD44+ Undetermined", "CD8+ T-Cell", "CLP", "Early Myeloid Progenitor", "Endosteal", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "GMP", "GMP/Myeloblast", "HSC", "HSPC", "Immature_B_Cell", "Intermediate Myeloid", "Macrophages", "Mature Myeloid", "MEP/Early Erythroblast", "Monocytes", "Adipo-MSC", "Non-Classical Monocytes", "pDC", "Plasma Cells", "Schwann Cells", "SEC", "SPINK2+ HSPC", "THY1+ MSC", "Undetermined", "VSMC"))

# read in AML objects
AML1_183 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/objects/codex_AML1_183_SeuratObj.RDS") 
AML1_382 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/objects/codex_AML1_382_SeuratObj.RDS")
AML2_191 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/objects/codex_AML2_191_SeuratObj.RDS")
AML2_380 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/objects/codex_AML2_380_SeuratObj.RDS")
AML3_1329 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/objects/codex_AML3_1329_SeuratObj.RDS")
AML3_1443 <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/objects/codex_AML3_1443_SeuratObj.RDS")

# correct orig.ident patient numbering system to be consistent (Penn/CHOP note format is AML(tan_lab_numbering)_(last digits of penn pathology numbering))
AML1_183$orig.ident <- "SB67_NBM51_AML1_183_CODEX_Mesmer"
AML2_191$orig.ident <- "SB67_NBM44_AML2_191_CODEX_Mesmer"

# export metadata
AML1_183_md <- AML1_183@meta.data
AML1_382_md <- AML1_382@meta.data
AML2_191_md <- AML2_191@meta.data
AML2_380_md <- AML2_380@meta.data
AML3_1329_md <- AML3_1329@meta.data
AML3_1443_md <- AML3_1443@meta.data

write.csv(AML1_183_md, "ReferenceMap_AML_Step6/seurat_metadata/AML1_183_md.csv")
write.csv(AML1_382_md, "ReferenceMap_AML_Step6/seurat_metadata/AML1_382_md.csv")
write.csv(AML2_191_md, "ReferenceMap_AML_Step6/seurat_metadata/AML2_191_md.csv")
write.csv(AML2_380_md, "ReferenceMap_AML_Step6/seurat_metadata/AML2_380_md.csv")
write.csv(AML3_1329_md, "ReferenceMap_AML_Step6/seurat_metadata/AML3_1329_md.csv")
write.csv(AML3_1443_md, "ReferenceMap_AML_Step6/seurat_metadata/AML3_1443_md.csv")

ob.list <- list(AML1_183,AML1_382,AML2_191,AML2_380, AML3_1329, AML3_1443)
for (i in 1:length(ob.list)) {
  # ob.list[[i]] <- subset(ob.list[[i]], cells = cells.use)
  ob.list[[i]] <- NormalizeData(object = ob.list[[i]], normalization.method = "CLR", margin = 1)
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  all_features = rownames(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], features = all_features, approx = FALSE)  
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:30, reduction.name = "UMAP_dim30", reduction.key = "UMAP_dim30_")
  ob.list[[i]] <- FindNeighbors(ob.list[[i]], features = all_features, k.param=30) #features = all_features
  ob.list[[i]] <- FindClusters(ob.list[[i]], algorithm = 2, resolution = 1.5) # Louvain
  print('finished processing {i}th samples')
} 

AML1_183 <- ob.list[[1]]
AML1_382 <- ob.list[[2]]
AML2_191 <- ob.list[[3]]
AML2_380 <- ob.list[[4]]
AML3_1329 <- ob.list[[5]]
AML3_1443 <- ob.list[[6]]

# merge objects - omitting 380 because of challenges in determining true MRD cells
AML.combined <- merge(x=AML1_183, y = c(AML1_382,AML2_191, AML3_1329, AML3_1443), add.cell.ids = c("AML1_183", "AML1_382", "AML2_191","AML3_1329", "AML3_1443"))

# Perform reference mapping
anchors <- FindTransferAnchors(
  reference = immune.combined, k.filter = NA,
  query = AML.combined,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:30
)
# saveRDS(anchors, "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/objects/AML_anchors_final_380removed_FINAL_040423.RDS")

AML.combined <- MapQuery(
  anchorset = anchors,
  query = AML.combined,
  reference = immune.combined,
  refdata = list(
    refmap = "cluster_anno_coarse"
  ),
  reference.reduction = "pca", 
  reduction.model = "UMAP_dim30"
)

AML.combined$classified_cluster_anno_coarse <- AML.combined$predicted.refmap
AML.combined$classified_cluster_anno_coarse_score <- AML.combined$predicted.refmap.score

AML.combined <- MapQuery(
  anchorset = anchors,
  query = AML.combined,
  reference = immune.combined,
  refdata = list(
    refmap = "cluster_anno_l2"
  ),
  reference.reduction = "pca", 
  reduction.model = "UMAP_dim30"
)

AML.combined$classified_cluster_anno_l2 <- AML.combined$predicted.refmap
AML.combined$classified_cluster_anno_l2_score <- AML.combined$predicted.refmap.score

# Rename orig.idents 
AML.combined$orig.ident <- factor(AML.combined$orig.ident, levels = c("SB67_NBM51_AML1_183_CODEX_Mesmer", "SB67_NBM46_AML1_382_CODEX_Mesmer", "SB67_NBM44_AML2_191_CODEX_Mesmer", "SB67_NBM52_AML3_1329_CODEX_Mesmer","SB67_NBM54_AML3_1443_CODEX_Mesmer"))
levels(AML.combined$orig.ident)
orig.ident_annotations <- c("AML1_Dx","AML1_MRD", "AML2_Dx","AML3_Dx", "AML3_MRD")
AML.combined <- SetIdent(AML.combined, value = "orig.ident")
names(orig.ident_annotations) <- levels(as.factor(AML.combined$orig.ident))
AML.combined <- RenameIdents(AML.combined, orig.ident_annotations)
AML.combined@meta.data["Sample_Name"] <- AML.combined@active.ident 

# Read in negative staging marrow object
NSM.combined <- readRDS("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_NSM_Step5/NSM.combined.mapped_unfiltered.RDS")

# Rename MSC to Adipo-MSC based on scRNA-Seq
# AdipoMSC_ids <- colnames(subset(All.combined_mapped, subset = classified_cluster_anno_l2 == "MSC")) 
# All.combined_mapped$classified_cluster_anno_l2[AdipoMSC_ids] <- "Adipo-MSC"

# Rename orig.idents 
levels(as.factor(NSM.combined$orig.ident))
orig.ident_annotations <- c("NSM_1996","NSM_1720", "NSM_1086")
NSM.combined <- SetIdent(NSM.combined, value = "orig.ident")
names(orig.ident_annotations) <- levels(as.factor(NSM.combined$orig.ident))
NSM.combined <- RenameIdents(NSM.combined, orig.ident_annotations)
NSM.combined@meta.data["Sample_Name"] <- NSM.combined@active.ident 

# Figure 7B Cell Composition of AML and NSM ----
## set coarse colors
names_coarse <- c("Adipocyte", "AEC", "Artifact", "Lymphoid", "HSPC", "Myeloid", "Endosteal", "Erythroid", "Megakaryocyte", "Macrophages", "Monocytes", "MSC", "pDC", "Schwann Cells", "SEC", "VSMC")
cols_coarse <- c("#FDDA0D", "#FFBF00", "#5A5A5A", "#AFE1AF", "#E0B0FF", "#A7C7E7", "#c3cede", "#BDB5D5", "#79127F", "#96C5D7", "#6495ED", "#FFB6C1", "#40B5AD", "#FF69B4", "#ff9d5c", "#DD3F4E")
names(cols_coarse) <- names_coarse

## set cal2 colors
cal2_cols <- c("#CF9FFF", "#E7C7DC", "#CAA7DD", "#A8A2D2", "#B6D0E2", "#2874A6", "#5599C8", "#AEC6CF", "#6495ED", "#64b8ed", "#96C5D7", "#40B5AD", "#8FCACA", "#CCE2CB", "#63BA97", "#7DB954", "#64864A", "#019477", "#953553", "#A1045A", "#a15891", "#9C58A1", "#79127F", "#BF40BF", "#FFD580", "#FFC8A2", "#FDDA0D", "#F3B0C3", "#FFBF00", "#ff9d5c", "#DD3F4E", "#FF69B4")
names(cal2_cols) <- c("HSC", "SPINK2+ HSPC", "HSPC", "GMP", "GMP/Myeloblast", "Early Myeloid Progenitor", "Intermediate Myeloid", "Mature Myeloid", "Monocytes", "Non-Classical Monocyte", "Macrophages", "pDC", "CLP", "Immature_B_Cell", "B-Cells", "CD4+ T-Cell", "CD8+ T-Cell", "Plasma Cells", "MEP/Early Erythroblast", "CD34+ CD61+", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "Adipo-MSC", "THY1+ MSC", "Adipocyte", "Endosteal", "AEC", "SEC", "VSMC", "Schwann Cells")
cell_order <- c("HSC", "SPINK2+ HSPC", "HSPC", "GMP", "GMP/Myeloblast", "Early Myeloid Progenitor", "Intermediate Myeloid", "Mature Myeloid", "Monocytes", "Non-Classical Monocyte", "Macrophages", "pDC", "CLP", "Immature_B_Cell", "B-Cells", "CD4+ T-Cell", "CD8+ T-Cell", "Plasma Cells", "MEP/Early Erythroblast", "CD34+ CD61+", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "Adipo-MSC", "THY1+ MSC", "Adipocyte", "Endosteal", "AEC", "SEC", "VSMC", "Schwann Cells")



## Create dataframe for plotting
NSM.combined$Sample_Name <- factor(NSM.combined$Sample_Name,  levels = c("NSM_1086", "NSM_1720", "NSM_1996"))
cell_counts_l2_AML <- as.data.frame(table(AML.combined$classified_cluster_anno_l2,AML.combined$Sample_Name))
cell_counts_l2_NSM <- as.data.frame(table(NSM.combined$classified_cluster_anno_l2,NSM.combined$Sample_Name))
cell_counts_l2 <- rbind(cell_counts_l2_NSM, cell_counts_l2_AML)
cell_counts_l2 <- dplyr::filter(cell_counts_l2, Var1 != "Artifact" & Var1 != "CD44+ Undetermined" & Var1 != "Undetermined" & Var1 != "Autofluorescent")
cell_counts_l2$Var1 <- factor(cell_counts_l2$Var1, levels = cell_order)
## Plot cluster_anno_l2 cell composition
cell_counts_l2 %>% arrange(Var1) %>% dplyr::filter(Var1 != "Artifact" & Var1 != "CD44+ Undetermined" & Var1 != "Undetermined" & Var1 != "Autofluorescent") %>% ggplot(aes(x=Var2,y=Freq, fill=droplevels(Var1))) +
  geom_bar(position="fill",stat="identity", width = 0.9) + 
  theme_bw() + 
  labs(x='Sample') + 
  labs(y='Cell Type') + scale_fill_manual(values=na.omit(cal2_cols)) + theme(axis.text = element_text(size=12)) + ggtitle("NSM/AML Predicted Cell Frequency") + RotatedAxis()

## Create combined object containing NSM and AML samples
AML1_183_mapped <- subset(AML.combined, subset = orig.ident == "SB67_NBM51_AML1_183_CODEX_Mesmer")
AML1_382_mapped <- subset(AML.combined, subset = orig.ident == "SB67_NBM46_AML1_382_CODEX_Mesmer")
AML2_191_mapped <- subset(AML.combined, subset = orig.ident == "SB67_NBM44_AML2_191_CODEX_Mesmer")
AML3_1329_mapped <- subset(AML.combined, subset = orig.ident == "SB67_NBM52_AML3_1329_CODEX_Mesmer")
AML3_1443_mapped <- subset(AML.combined, subset = orig.ident == "SB67_NBM54_AML3_1443_CODEX_Mesmer")
All.combined <- merge(AML1_183_mapped, c(AML1_382_mapped, AML2_191_mapped, AML3_1329_mapped, AML3_1443_mapped, NSM.combined))

## Perform proportion test to confirm significance of myeloid populations enriched in AML samples vs. NSM ----
myeloid_cell_types <- c("GMP", "GMP/Myeloblast", "Early Myeloid Progenitor", "Intermediate Myeloid", "Monocyte","Non-Classical Monocyte", "pDC", "Macrophage") # make sure there are non npm1 mutant blast annotation since we are just checking the unsupervised labeling here, also removed mature myeloid since AML should be immature myeloid or monocytic/dendritic lineage
All.combined@meta.data <- mutate(All.combined@meta.data, IsMyeloid = classified_cluster_anno_l2 %in% myeloid_cell_types)
cell_counts_l2 <- mutate(cell_counts_l2, IsAML = Var2 %in% orig.ident_annotations)
table(subset(All.combined, subset = Sample_Name %in% c("AML1_Dx","AML2_Dx", "AML3_Dx"))$IsMyeloid)
table(subset(All.combined, subset = Sample_Name %in% c("NSM_1086", "NSM_1720", "NSM_1996"))$IsMyeloid)
prop.test(x= c(38108, 32945), n = c((71701+38108), (132355+32945)))
# results - p <2.2e-16 , prop1 = 0.321 , prop2 = 0.199

# Map mutant cells ----
## Load classified spatially joined dataframes to map mutant cell IDs
df_classified_183 <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/qupath_results/df_classified_bonedistsSB67_NBM51_AML1_183_CODEX_Mesmer.csv")
df_classified_382 <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/qupath_results/df_classified_bonedistsSB67_NBM46_AML1_382_CODEX_Mesmer.csv")
df_classified_191 <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/qupath_results/df_classified_bonedistsSB67_NBM44_AML2_191_CODEX_Mesmer.csv")
df_classified_1329 <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/qupath_results/df_classified_bonedistsSB67_NBM52_AML3_1329_CODEX_Mesmer.csv")
df_classified_1443 <- read_csv("/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/qupath_results/df_classified_bonedistsSB67_NBM54_AML3_1443_CODEX_Mesmer.csv")

## Pull out columns of interest from spatial join results
df_classified_382 <- df_classified_382 %>% distinct(CellID, .keep_all = TRUE) # 60/68876 cells with ties in classified_382, removed by taking unique Cell IDs
df_classified_183 <- tibble(df_classified_183 %>% dplyr::select(CellID, Class, `Distance to annotation with Bone µm`))
df_classified_382 <- tibble(df_classified_382 %>% dplyr::select(CellID, Class, `Distance to annotation with Bone µm`))
df_classified_191 <- tibble(df_classified_191 %>% dplyr::select(CellID, Class, `Distance to annotation with Bone µm`))
df_classified_1329 <- tibble(df_classified_1329 %>% dplyr::select(CellID, Class, `Distance to annotation with Bone µm`))
df_classified_1443 <- tibble(df_classified_1443 %>% dplyr::select(CellID, Class, `Distance to annotation with Bone µm`))
## Pull out cell IDs for each patient
AML1_183_cellids <- rownames(AML1_183_mapped@meta.data)
AML1_382_cellids <- rownames(AML1_382_mapped@meta.data)
AML2_191_cellids <- rownames(AML2_191_mapped@meta.data)
AML3_1329_cellids <- rownames(AML3_1329_mapped@meta.data)
AML3_1443_cellids <- rownames(AML3_1443_mapped@meta.data)
## join the matrices
AML1_183_mapped@meta.data <- left_join(AML1_183_mapped@meta.data, df_classified_183, by = "CellID")
AML1_382_mapped@meta.data <- left_join(AML1_382_mapped@meta.data, df_classified_382, by = "CellID")
AML2_191_mapped@meta.data <- left_join(AML2_191_mapped@meta.data, df_classified_191, by = "CellID")
AML3_1329_mapped@meta.data <- left_join(AML3_1329_mapped@meta.data, df_classified_1329, by = "CellID")
AML3_1443_mapped@meta.data <- left_join(AML3_1443_mapped@meta.data, df_classified_1443, by = "CellID")
## rename the rows appropriately
rownames(AML1_183_mapped@meta.data) <- AML1_183_cellids
rownames(AML1_382_mapped@meta.data) <- AML1_382_cellids
rownames(AML2_191_mapped@meta.data) <- AML2_191_cellids
rownames(AML3_1329_mapped@meta.data) <- AML3_1329_cellids
rownames(AML3_1443_mapped@meta.data) <- AML3_1443_cellids

# UMAP and see if the cells cluster
DimPlot(AML1_183_mapped, group.by = "Class") + coord_fixed()

# filter cells of HSPC, myeloid, erythroid, and meg lineage

AML1_183_mapped@meta.data <- AML1_183_mapped@meta.data %>% mutate(Adjusted_Class = if_else(Class == "NPM1_Mutant" & classified_cluster_anno_l2 %in% c("Early Myeloid Progenitor", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "GMP", "GMP/Myeloblast", "HSPC", "Intermediate Myeloid", "Macrophages", "MEP/Early Erythroblast", "Monocytes", "Non-Classical Monocyte", "pDC", "SPINK2+ HSPC"), true = "NPM1_Mutant",false =  "WT", missing = "missing"))
table(AML1_183_mapped@meta.data$Adjusted_Class, AML1_183_mapped@meta.data$classified_cluster_anno_l2)
AML1_382_mapped@meta.data <- AML1_382_mapped@meta.data %>% mutate(Adjusted_Class = if_else(Class == "NPM1_Mutant" & classified_cluster_anno_l2 %in% c("Early Myeloid Progenitor", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "GMP", "GMP/Myeloblast", "HSPC", "Intermediate Myeloid", "Macrophages", "MEP/Early Erythroblast", "Monocytes", "Non-Classical Monocyte", "pDC", "SPINK2+ HSPC"), true = "NPM1_Mutant",false =  "WT", missing = "missing"))
table(AML1_382_mapped@meta.data$Adjusted_Class, AML1_382_mapped@meta.data$classified_cluster_anno_l2)
AML2_191_mapped@meta.data <- AML2_191_mapped@meta.data %>% mutate(Adjusted_Class = if_else(Class == "NPM1_Mutant" & classified_cluster_anno_l2 %in% c("Early Myeloid Progenitor", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "GMP", "GMP/Myeloblast", "HSPC", "Intermediate Myeloid", "Macrophages", "MEP/Early Erythroblast", "Monocytes", "Non-Classical Monocyte", "pDC", "SPINK2+ HSPC"), true = "NPM1_Mutant",false =  "WT", missing = "missing"))
table(AML2_191_mapped@meta.data$Adjusted_Class, AML2_191_mapped@meta.data$classified_cluster_anno_l2)
AML3_1329_mapped@meta.data <- AML3_1329_mapped@meta.data %>% mutate(Adjusted_Class = if_else(Class == "NPM1_Mutant" & classified_cluster_anno_l2 %in% c("Early Myeloid Progenitor", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "GMP", "GMP/Myeloblast", "HSPC", "Intermediate Myeloid", "Macrophages", "MEP/Early Erythroblast", "Monocytes", "Non-Classical Monocyte", "pDC", "SPINK2+ HSPC"), true = "NPM1_Mutant",false =  "WT", missing = "missing"))
table(AML3_1329_mapped@meta.data$Adjusted_Class, AML3_1329_mapped@meta.data$classified_cluster_anno_l2)
AML3_1443_mapped@meta.data <- AML3_1443_mapped@meta.data %>% mutate(Adjusted_Class = if_else(Class == "NPM1_Mutant" & classified_cluster_anno_l2 %in% c("Early Myeloid Progenitor", "Erythroblast", "Erythroid", "GATA1neg_Mks", "GATA1pos_Mks", "GMP", "GMP/Myeloblast", "HSPC", "Intermediate Myeloid", "Macrophages", "MEP/Early Erythroblast", "Monocytes", "Non-Classical Monocyte", "pDC", "SPINK2+ HSPC"), true = "NPM1_Mutant",false =  "WT", missing = "missing"))
table(AML3_1443_mapped@meta.data$Adjusted_Class, AML3_1443_mapped@meta.data$classified_cluster_anno_l2)

# save objects
#saveRDS(AML1_183_mapped, "objects/AML1_183_mapped.RDS")
#saveRDS(AML1_382_mapped, "objects/AML1_382_mapped.RDS")
#saveRDS(AML2_191_mapped, "objects/AML2_191_mapped.RDS")
#saveRDS(AML3_1329_mapped, "objects/AML3_1329_mapped.RDS")
#saveRDS(AML3_1443_mapped, "objects/AML3_1443_mapped.RDS")

# read objects
AML1_183_mapped <- readRDS("objects/AML1_183_mapped.RDS")
AML1_382_mapped <- readRDS("objects/AML1_382_mapped.RDS")
AML2_191_mapped <- readRDS("objects/AML2_191_mapped.RDS")
AML3_1329_mapped<- readRDS("objects/AML3_1329_mapped.RDS")
AML3_1443_mapped <- readRDS("objects/AML3_1443_mapped.RDS")



# combined objects 
All.combined_mapped <- merge(AML1_183_mapped, c(AML1_382_mapped, AML2_191_mapped, AML3_1329_mapped, AML3_1443_mapped, NSM.combined))
AML.combined_mapped <- merge(AML1_183_mapped, c(AML1_382_mapped, AML2_191_mapped, AML3_1329_mapped, AML3_1443_mapped))

# change cluster anno l2 to blasts
Blast_ids <- colnames(subset(All.combined_mapped, subset = Adjusted_Class == "NPM1_Mutant")) # Blast Rename
All.combined_mapped$classified_cluster_anno_l2[Blast_ids] <- "NPM1 Mutant Blast"
Blast_ids <- colnames(subset(AML.combined_mapped, subset = Adjusted_Class == "NPM1_Mutant")) # Blast Rename
AML.combined_mapped$classified_cluster_anno_l2[Blast_ids] <- "NPM1 Mutant Blast"

# save objects
## saveRDS(All.combined_mapped, file = "objects/All.combined_mapped_blastslabeled_040423.RDS")
## saveRDS(AML.combined_mapped, file = "objects/AML.combined_mapped_blastslabeled_040423.RDS")
## saveRDS(All.combined, file = "objects/All.combined_mapped_unlabeled_040423.RDS")

# Related to Figure 7C
FeatureScatter(AML3_1329_mapped, feature1 = "x.coord", feature2 = "y.coord", group.by = "classified_cluster_anno_l2", cols = cal2_cols, pt.size = 0.5)
FeatureScatter(subset(All.combined, Sample_Name == "NSM_1720") , feature1 = "x.coord", feature2 = "y.coord", group.by = "classified_cluster_anno_l2", cols = cal2_cols, pt.size = 0.1)

# Figure 7D ----
# make a bar plot just showing the MSC frequencies and perform t.test
All_ct_freqs <- All.combined_mapped@meta.data %>%
  group_by(Sample_Name, classified_cluster_anno_l2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

All_ct_freqs <- All_ct_freqs %>% 
  mutate(Sample_Group = substr(Sample_Name, start = 1, stop = 3))

All_ct_freqs_toplot <- All_ct_freqs %>% group_by(Sample_Group, classified_cluster_anno_l2) %>% summarise(freq=freq, mean_freq = mean(freq), sd = sd(freq)) 
All_ct_freqs_toplot <- All_ct_freqs_toplot %>% dplyr::filter(classified_cluster_anno_l2 == "MSC" | classified_cluster_anno_l2 == "THY1+ MSC")
All_ct_freqs_toplot %>% ggplot(aes(x = classified_cluster_anno_l2, y = freq, fill = Sample_Group)) + 
  geom_boxplot()

All_ct_freqs_toplot$Sample_Group <- factor(All_ct_freqs_toplot$Sample_Group, levels = c("NSM", "AML"))
  
p1 <- All_ct_freqs_toplot %>% ggplot(aes(x = classified_cluster_anno_l2, y = mean_freq, fill = Sample_Group)) + 
  geom_bar(stat="identity", position = "dodge", color = 'black') +
  geom_errorbar(aes(ymin=mean_freq-sd, ymax=mean_freq+sd), width=.2,
                position=position_dodge(.9)) + scale_fill_manual(values = c("#89CFF0", "#FA8072")) + theme_minimal()

# ggsave(p1, device = "pdf",filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/Figures/Figure7D_MSC_SubTypeFrequency.pdf", height = 5, width = 7.5)


aml_sample_names <- c("AML1_Dx", "AML1_MRD", "AML2_Dx", "AML3_Dx", "AML3_MRD")
NSM_sample_names <- c("NSM_1086", "NSM_1720", "NSM_1996")

aml_sample_names <- c("AML1_Dx", "AML1_MRD", "AML2_Dx", "AML3_Dx", "AML3_MRD")
NSM_sample_names <- c("NSM_1086", "NSM_1720", "NSM_1996")

## t-test for seeing whether the frequency distribution is different between AML and NSM for both Adipo and THY1+ MSCs ----
t.test(subset(All_ct_freqs, subset = Sample_Group == "AML" & (classified_cluster_anno_l2 == "MSC"))$freq, subset(All_ct_freqs, subset = Sample_Group == "NSM" & (classified_cluster_anno_l2 == "MSC"))$freq)
# p=0.004654
t.test(subset(All_ct_freqs, subset = Sample_Group == "AML" & (classified_cluster_anno_l2 == "THY1+ MSC"))$freq, subset(All_ct_freqs, subset = Sample_Group == "NSM" & (classified_cluster_anno_l2 == "THY1+ MSC"))$freq)
# p = 0.03737


# Analyze HIF1a levels Supplemental Figure S7F ----
cal2_cols_npm1 <- c(cal2_cols, "#FA8072")
names(cal2_cols_npm1) <- c(cell_order, "NPM1 Mutant Blast")

p1 <- VlnPlot(subset(All.combined_mapped,classified_cluster_anno_l2 %in% c("Mature Myeloid", "Intermediate Myeloid", "Early Myeloid Progenitor","GMP/Myeloblast","GMP", "NPM1 Mutant Blast")), features = c("codex_HIF1A"), slot = "data", pt.size = 0 , cols = cal2_cols, group.by = "classified_cluster_anno_l2", sort = TRUE) + NoLegend()
#ggsave(p1, device = "pdf",filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/Figures/SuppS6C_NPM1_HIF1A_VlnPlot.pdf", height = 5, width = 7.5)

# BCL2 levels Supplemental Figure S7G ----
p1 <- VlnPlot(subset(All.combined_mapped,classified_cluster_anno_l2 %in% c("NPM1 Mutant Blast")), split.by = "Sample_Name", cols = c("#87CEEB", "#20639B", "#3CAEA3","#F6D55C", "#ED553B"), features = c("codex_BCL2"), slot = "data", pt.size = 0, group.by = "classified_cluster_anno_l2", sort = TRUE) 
ggsave(p1, device = "pdf",filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/Figures/SuppS6C_NPM1_BCL2_VlnPlot.pdf", height = 5, width = 7.5)

# Complex IV levels Supplemental Figure S7H ----
p1 <- VlnPlot(subset(All.combined_mapped,classified_cluster_anno_l2 %in% c("NPM1 Mutant Blast")), group.by = "Sample_Name", features = c("codex_OXPHOS"), slot = "data", pt.size = 0, cols = c("#87CEEB", "#20639B", "#3CAEA3","#F6D55C", "#ED553B"),  sort = FALSE) 
ggsave(p1, device = "pdf",filename = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/Figures/SuppS6C_NPM1_ComplexIV_VlnPlot.pdf", height = 5, width = 7.5)

# export for nb analysis
# Add metadata columns and write annotated CSV for neighborhood analysis
All.combined_mapped@meta.data$x.coord_rev <- -(All.combined_mapped$x.coord)
All.combined_mapped@meta.data$y.coord_rev <- -(All.combined_mapped$y.coord)
All.combined_mapped@meta.data <- mutate(All.combined_mapped@meta.data, Region = "reg001")
# write.csv(x = All.combined_mapped@meta.data, file = "/mnt/isilon/tan_lab_imaging/Analysis/bandyopads/NBM_CODEX_Atlas/Combined_Analysis/Seurat/ReferenceMap_AML_Step6/For_Neighborhoods/NSM_AML_combined_annotated_ForNeighborhoods_040523.csv")


