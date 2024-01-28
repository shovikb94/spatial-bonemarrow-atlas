library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

# load combined scRNA-Seq object and prepare CellChat Obj
scrna_cal2_cols <- c("#FFD580", "#FFBF00", "#B0DFE6", "#7DB954", "#64864A", "#8FCACA", "#4682B4", "#CAA7DD", "#B6D0E2", "#a15891", "#FED7C3", "#A8A2D2", "#CF9FFF", "#9C58A1", "#2874A6", "#96C5D7", "#63BA97", "#BF40BF", "#953553", "#6495ED", "#E7C7DC", "#5599C8", "#FA8072", "#F3B0C3", "#F89880", "#40B5AD", "#019477", "#97C1A9", "#C6DBDA", "#CCE2CB", "#79127F", "#ff9d5c", "#FFC8A2", "#DD3F4E")
cal2_col_names <- c("Adipo-MSC", "AEC", "Ba/Eo/Ma", "CD4+ T-Cell", "CD8+ T-Cell", "CLP", "Cycling DCs", "Cycling HSPC", "Early Myeloid Progenitor", "Erythroblast", "Fibro-MSC", "GMP", "HSC", "Late Erythroid", "Late Myeloid", "Macrophages", "Mature B", "Megakaryocyte", "MEP", "Monocyte", "MPP", "Neutrophil", "Osteo-MSC", "Osteoblast", "OsteoFibro-MSC", "pDC", "Plasma Cell", "Pre-B", "Pre-Pro B", "Pro-B", "RBC", "SEC", "THY1+ MSC", "VSMC")
names(scrna_cal2_cols) <- cal2_col_names
cell_order <- c("HSC", "MPP", "Cycling HSPC", "MEP", "Erythroblast", "Late Erythroid", "RBC", "Megakaryocyte", "GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil", "Monocyte", "Macrophages", "Ba/Eo/Ma", "Cycling DCs", "pDC", "CLP", "Pre-Pro B", "Pro-B", "Pre-B", "Mature B", "Plasma Cell", "CD4+ T-Cell", "CD8+ T-Cell", "AEC", "SEC", "VSMC", "THY1+ MSC", "Adipo-MSC", "OsteoFibro-MSC", "Osteo-MSC", "Osteoblast", "Fibro-MSC")


combined <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Seurat/SB66_combined_corrected.RDS")
combined <- subset(combined, subset = cluster_anno_l2 != "NKT Cell" & cluster_anno_l2 != "RNAlo MSC") # remove RNAlo MSCs and likely artifactual NKT cell population
combined$cluster_anno_l2 <- droplevels(combined$cluster_anno_l2)
combined$cluster_anno_l2 <- factor(combined$cluster_anno_l2, levels = c("HSC", "MPP", "Cycling HSPC", "MEP", "Erythroblast", "Late Erythroid", "RBC", "Megakaryocyte", "GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil", "Monocyte", "Macrophages", "Ba/Eo/Ma", "Cycling DCs", "pDC", "CLP", "Pre-Pro B", "Pro-B", "Pre-B", "Mature B", "Plasma Cell", "CD4+ T-Cell", "CD8+ T-Cell", "AEC", "SEC", "VSMC", "THY1+ MSC", "Adipo-MSC", "OsteoFibro-MSC", "Osteo-MSC", "Osteoblast", "Fibro-MSC"))
cell_order <- c("HSC", "MPP", "Cycling HSPC", "MEP", "Erythroblast", "Late Erythroid", "RBC", "Megakaryocyte", "GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil", "Monocyte", "Macrophages", "Ba/Eo/Ma", "Cycling DCs", "pDC", "CLP", "Pre-Pro B", "Pro-B", "Pre-B", "Mature B", "Plasma Cell", "CD4+ T-Cell", "CD8+ T-Cell", "AEC", "SEC", "VSMC", "THY1+ MSC", "Adipo-MSC", "OsteoFibro-MSC", "Osteo-MSC", "Osteoblast", "Fibro-MSC")
DimPlot(combined, group.by = "cluster_anno_l2", cols = scrna_cal2_cols, label = TRUE, repel= TRUE)+ NoLegend()

cellchat <- createCellChat(object = combined, group.by = "cluster_anno_l2")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4) # do parallel
options(future.globals.maxSize= 891289600, future.seed=TRUE)



# Run CellChat ----
## Identify Genes for CellChat ----
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

## Project expression to PPI network to impute expression ----
## project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)

## Compute Interaction Network  ----
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

## convert cellchat to DF 
df.cellchat <- subsetCommunication(cellchat)

## Supplemental Table S4 ----
write.csv(df.cellchat, "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/CellChat/Figures/cellchat_LRpairs_updated.csv")

## Calculate Network Metrics ----
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways


saveRDS(cellchat, "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/CellChat/cellchat_032923.RDS")
cellchat <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/CellChat/cellchat_032923.RDS")


# Plot relevant cellchat plots ----

# Show general signaling patterns for each lineage
sources.use <- c("Adipo-MSC", "THY1+ MSC", "Fibro-MSC","Osteo-Fibro MSC","Osteo-MSC","Osteoblast","VSMC", "AEC", "SEC")

targets_hspc <- c("HSC", "GMP", "CLP", "MEP")
targets_hspc <- c("HSC","MPP", "Cycling HSPC","GMP", "Early Myeloid Progenitor", "CLP", "MEP")

targets_myeloid <- c("GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil")
targets_mege <- c("MEP", "Megakaryocyte", "Erythroblast", "Late Erythroid", "RBC")
targets_lymphoid <- c("CLP", "CD4+ T-Cell", "CD8+ T-Cell", "Pre-Pro B", "Pre-B", "Pro-B", "Mature B")



## Plot Chord Diagrams for Specific Pairs Figure 3B ----
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = "CXCL", geneLR.return = FALSE)
pathways.show <- "CXCL"
LR.show <- pairLR.CXCL[4,] # show one ligand-receptor pair
netVisual_individual(cellchat, signaling = pathways.show,  remove.isolate = TRUE, pairLR.use = LR.show, sources.use = sources.use, targets.use = targets_hspc, vertex.receiver = vertex.receiver)
pathways.show <- "SELE"
netVisual_individual(cellchat, signaling = pathways.show,  remove.isolate = TRUE, pairLR.use = "SELE_CD44", sources.use = sources.use, targets.use = targets_hspc, vertex.receiver = vertex.receiver)

# Add choird diagrams for specific interactions to validate later Supplemental Figure S8B ---- 
pathways.show <- "COLLAGEN"
sources.use <- c("Osteoblast")
targets_hspc <- c("Monocyte")
netVisual_individual(cellchat, signaling = pathways.show,  remove.isolate = TRUE, pairLR.use = "COL1A1_CD44", sources.use = sources.use, targets.use = targets_hspc, vertex.receiver = vertex.receiver)

pathways.show <- "COLLAGEN"
sources.use <- c("Osteoblast")
targets_hspc <- c("Monocyte")
netVisual_individual(cellchat, signaling = pathways.show,  remove.isolate = TRUE, pairLR.use = "COL1A1_CD44", sources.use = sources.use, targets.use = targets_hspc, vertex.receiver = vertex.receiver)

pathways.show <- "ANGPT"
sources.use <- c("MEP", "Erythroblast", "Megakaryocyte")
targets_hspc <- c("AEC", "SEC")
netVisual_individual(cellchat, signaling = pathways.show,  remove.isolate = TRUE, pairLR.use = "ANGPT1_TEK", sources.use = sources.use, targets.use = targets_hspc, vertex.receiver = vertex.receiver)

pathways.show <- "IL2"
sources.use <- c("Adipo-MSC", "THY1+ MSC", "Fibro-MSC","Osteo-Fibro MSC","Osteo-MSC","Osteoblast","VSMC", "AEC", "SEC")
targets_hspc <- c("CLP")
netVisual_individual(cellchat, signaling = pathways.show,  remove.isolate = TRUE, pairLR.use = "IL7_IL7R_IL2RG", sources.use = sources.use, targets.use = targets_hspc, vertex.receiver = vertex.receiver)

pathways.show <- "NOTCH"
sources.use <- c("Adipo-MSC", "THY1+ MSC", "Fibro-MSC","Osteo-Fibro MSC","Osteo-MSC","Osteoblast","VSMC", "AEC", "SEC")
targets_hspc <- c("AEC","VSMC")
netVisual_individual(cellchat, signaling = pathways.show,  remove.isolate = TRUE, pairLR.use = "DLL4_NOTCH2", sources.use = sources.use, targets.use = targets_hspc, vertex.receiver = vertex.receiver)



## Plot Signaling Biaxial Plot Figure 3D ----
netAnalysis_signalingRole_scatter(cellchat, color.use = scrna_cal2_cols,  font.size = 14, do.label = TRUE)

## Plot Pathway Level Heatmaps Figure 3E-F ----
pathways.show <- c("CXCL","PTN", "KIT", "ANGPT", "CSF","NOTCH", "FGF","PDGF","VWF", "NGF", "NCAM","RANKL","CSPG4","SELE", "SELPLG","IL1","IL2","IL6","TNF","FLT3","TGFb", "MHC-I", "MHC-II", "WNT", "VCAM", "COLLAGEN")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling =  pathways.show, color.heatmap = "Purples", color.use = scrna_cal2_cols[cell_order])
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling =  pathways.show, color.heatmap = "Purples", color.use = scrna_cal2_cols[cell_order])
ht1 + ht2

## Plot NMF Results Figure 3G, S4E ----
selectK(cellchat, pattern = "outgoing")
nPatterns = 7
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, width = 14, height = 24, color.use = scrna_cal2_cols[cell_order])
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, width = 14, height = 24, color.use = scrna_cal2_cols[cell_order])



# extra 
#plot heatmap of edge weights
library(ComplexHeatmap)
weights <- cellchat@net$weight
Heatmap(scale(t(weights)))

#plot heatmap of interaction counts
library(ComplexHeatmap)
count <- cellchat@net$count
Heatmap(scale(t(count)))




