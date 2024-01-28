library(Seurat)
library(AUCell)
library(GSEABase)

setwd("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Inflammation_Comparison")

combined <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Seurat/combined_corrected_NKTremoved.RDS")
combined <- RunUMAP(combined, dims = 1:30, reduction = "pca", reduction.name = "UMAP_dim30", return.model = TRUE) # rerun find umap to return the umap model
hca <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Inflammation_Comparison/hcabm40k.hca.RDS")
peds <- readRDS("/mnt/isilon/tan_lab/chenc6/MLLr_Project/HTAN_Data_Submission/scRNASeq_Level4/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC_Oct25_2021.rds")
peds = UpdateSeuratObject(object = peds)
peds$orig.ident <- peds$sample
combined_aspirate <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/combined_aspirate.RDS")
Li_MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/ExternalDatasets_Revisions/mapped_mscs/Li_MSCs.RDS")

# merge objs
combined$Dataset <- "SB"
peds$Dataset <- "CC_Peds"
hca$Dataset <- "HCA_Adult"
combined_aspirate$Dataset <- "SB_aspirate"
Li_MSCs$Dataset <- "Li_MSCs" 

merged_atlases <- merge(combined, y = c(hca, peds, combined_aspirate, Li_MSCs), add.cell.ids = c("SB_femoralhead", "HCA", "CC_Peds", "SB_aspirate", "Li_aspirate"))

# Supplemental Figure 2A Inflammation ----- 
#Run AUCell
combined.small <- subset(combined, downsample = 300) # note the results may be subtly different between runs because of downsampling step
fmeta <- data.frame(symbol = rownames(combined.small)) 
rownames(fmeta) <- fmeta$symbol
eset_combined.small <- new("ExpressionSet",
                           assayData = assayDataNew("environment", exprs=combined.small@assays$RNA@counts, 
                                                    norm_exprs = combined.small@assays$RNA@data),
                           phenoData =  new("AnnotatedDataFrame", data = combined.small@meta.data),
                           featureData = new("AnnotatedDataFrame", data = fmeta))
exprMatrix <- exprs(eset_combined.small)
geneSets <- getGmt("~/Downloads/h.all.v2022.1.Hs.symbols.gmt")
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets)))
# add random noise
# Random
set.seed(321)
extraGeneSets <- c(
  GeneSet(sample(rownames(exprMatrix), 50), setName="Random (50g)"),
  GeneSet(sample(rownames(exprMatrix), 500), setName="Random (500g)"))

countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
# Housekeeping-like
extraGeneSets <- c(extraGeneSets,
                   GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK-like (100g)"))

geneSets <- GeneSetCollection(c(geneSets,extraGeneSets))
names(geneSets)
# Start to build the cell rankings
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=FALSE)
#save(cells_rankings, file="AUCell/cells_rankings.RData")

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
#save(cells_AUC, file="AUCell/AUCecells_AUC.RData")

set.seed(123)
par(mfrow=c(3,3)) 
#cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=FALSE, assign=TRUE)
#save(cells_assignment, file="AUCell/cellassignment_AUC.RData")

auc_test <- as.data.frame(t(cells_AUC@assays@data$AUC))
combined_AUC <- AddMetaData(combined.small, auc_test)


FeaturePlot(object = combined_AUC, features = c("HALLMARK_INFLAMMATORY_RESPONSE"), coord.fixed = TRUE) & NoAxes()




## Supplemental Figure S2A Inflammation UMAP----- 
### For details on AUCell run navigate to Figure 5F code, the same combined_AUC object is used
### note due to computationally necessary downsampling results may vary slightly between runs
combined_AUC <- subset(combined_AUC, cells = rownames(auc_test))
p1 <- FeaturePlot(object = combined_AUC, raster = TRUE, raster.dpi = c(1028,1028), features = c("HALLMARK_INFLAMMATORY_RESPONSE"),, coord.fixed = TRUE) & NoAxes()

p1 <- FeaturePlot(object = combined_AUC, raster = TRUE, raster.dpi = c(1028,1028), features = c("HALLMARK_INFLAMMATORY_RESPONSE"),, coord.fixed = TRUE) & NoAxes()
p1_raster <- rasterize(p1)

ggsave(p1_raster,file = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure1/InflammatoryResponse_UMAP.pdf", device = "pdf", width = 15, height = 2)

# Supplemental Figure S2B Inflammation Per Cell Type----
p2 <- VlnPlot(object = combined_AUC, features = c("HALLMARK_INFLAMMATORY_RESPONSE"), group.by = "cluster_anno_l2",  cols = cal2_cols ,pt.size = 0, sort = TRUE) & NoLegend()
ggsave(p2, file = "~/Documents/Manuscripts/NBM_Atlas/Figures/Supplemental_Figures/Related_to_Figure1/InflammationResponse_Violins.pdf", device = "pdf", width = 10, height = 5)


VlnPlot(combined, fill.by = "ident", stack = TRUE, features = c("IFNG", "TNF", "IL6", "CXCL8", "IL1B"), group.by = "cluster_anno_l2", flip = TRUE, cols = cal2_cols)
azimuth <- NormalizeData(azimuth)
VlnPlot(azimuth, assay = "refAssay", fill.by = "ident", stack = TRUE, features = c("IFNG", "TNF", "IL6", "CXCL8", "IL1B"), group.by = "celltype.l2", flip = TRUE, cols = cal2_cols)

# Supplemental Figure S2C Inflammatory Cytokines Across Datasets ----
VlnPlot(merged_atlases, group.by = "Dataset", features = c("IL1B", "IFNG", "TNF", "IL6", "CSF2", "IL18"), pt.size = 0)

# Prepare AUCell

geneset_path <- "h.all.v2022.1.Hs.symbols.gmt"

runAUCell <- function(object, geneset_path, cluster_anno, n_cells){
object <- SetIdent(object, value = "orig.ident")
combined.small <- subset(object, downsample = n_cells) # note the results may be subtly different between runs because of downsampling step
fmeta <- data.frame(symbol = rownames(combined.small)) 
rownames(fmeta) <- fmeta$symbol
eset_combined.small <- new("ExpressionSet",
                           assayData = assayDataNew("environment", exprs=combined.small@assays$RNA@counts, 
                                                    norm_exprs = combined.small@assays$RNA@data),
                           phenoData =  new("AnnotatedDataFrame", data = combined.small@meta.data),
                           featureData = new("AnnotatedDataFrame", data = fmeta))
exprMatrix <- exprs(eset_combined.small)
geneSets <- getGmt(geneset_path)
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets)))
# add random noise
# Random
set.seed(321)
extraGeneSets <- c(
  GeneSet(sample(rownames(exprMatrix), 50), setName="Random (50g)"),
  GeneSet(sample(rownames(exprMatrix), 500), setName="Random (500g)"))

countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
# Housekeeping-like
extraGeneSets <- c(extraGeneSets,
                   GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK-like (100g)"))

geneSets <- GeneSetCollection(c(geneSets,extraGeneSets))
names(geneSets)

# Start to build the cell rankings
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=FALSE)
save(cells_rankings, file="AUCell/peds_cells_rankings.RData")

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file="AUCell/peds_AUCecells_AUC.RData")

set.seed(123)
par(mfrow=c(3,3)) 
#cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=FALSE, assign=TRUE)
# save(cells_assignment, file="AUCell/cellassignment_AUC.RData")

auc_test <- as.data.frame(t(cells_AUC@assays@data$AUC))
combined_AUC <- AddMetaData(combined.small, auc_test)
return(combined_AUC) 
}


# read in inflamed and control MSCs from GSE115149
inflamed_mscs <- Read10X_h5(filename = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Inflammation_Comparison/GSE115149_filtered_gene_bc_matrices_h5.h5")
inflamed_mscs <- CreateSeuratObject(inflamed_mscs, project = "inflamed_MSCs", assay = "RNA", min.cells = 3, min.features = 100)

EM17_lic <- Read10X(data.dir = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Inflammation_Comparison/SeparateSampleMatrices/EM17_lic/filtered_gene_bc_matrices/hg19/")
EM17_ctrl <- Read10X(data.dir = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Inflammation_Comparison/SeparateSampleMatrices/EM17_un/filtered_gene_bc_matrices/hg19/")
EM18_lic <- Read10X(data.dir = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Inflammation_Comparison/SeparateSampleMatrices/EM18_lic/filtered_gene_bc_matrices/hg19/")
EM18_ctrl <- Read10X(data.dir = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Inflammation_Comparison/SeparateSampleMatrices/EM18_un/filtered_gene_bc_matrices/hg19/")
EM01_lic <- Read10X(data.dir = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Inflammation_Comparison/SeparateSampleMatrices/EM1_lic/filtered_gene_bc_matrices/hg19/")
EM01_ctrl <- Read10X(data.dir = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Inflammation_Comparison/SeparateSampleMatrices/EM1_un/filtered_gene_bc_matrices/hg19/")

EM17_lic <- CreateSeuratObject(EM17_lic, project = "EM17_lic", assay = "RNA", min.cells = 3, min.features = 100)
EM17_ctrl <- CreateSeuratObject(EM17_ctrl, project = "EM17_ctrl", assay = "RNA", min.cells = 3, min.features = 100)
EM18_lic <- CreateSeuratObject(EM18_lic, project = "EM18_lic", assay = "RNA", min.cells = 3, min.features = 100)
EM18_ctrl <- CreateSeuratObject(EM18_ctrl, project = "EM18_ctrl", assay = "RNA", min.cells = 3, min.features = 100)
EM01_lic <- CreateSeuratObject(EM01_lic, project = "EM01_lic", assay = "RNA", min.cells = 3, min.features = 100)
EM01_ctrl <- CreateSeuratObject(EM01_ctrl, project = "EM01_ctrl", assay = "RNA", min.cells = 3, min.features = 100)

inflamed_mscs <- merge(EM17_lic, y = c(EM18_lic, EM01_lic))
ctrl_lic_mscs <-  merge(EM17_ctrl, y = c(EM18_ctrl, EM01_ctrl))

ob.list <- list(inflamed_mscs, ctrl_lic_mscs)
for (i in 1:length(ob.list)) {
  ob.list[[i]][["percent.mt"]] <- PercentageFeatureSet(ob.list[[i]], pattern = "^MT-")
  ob.list[[i]] <- subset(ob.list[[i]], subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 10)
  # ob.list[[i]] <- subset(ob.list[[i]], cells = cells.use)
  ob.list[[i]] <- NormalizeData(ob.list[[i]])
  ob.list[[i]] <- FindVariableFeatures(ob.list[[i]])
  ob.list[[i]] <- ScaleData(ob.list[[i]])
  ob.list[[i]] <- RunPCA(ob.list[[i]], features = VariableFeatures(object = ob.list[[i]]))
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:30, reduction.name = "UMAP_dim30", reduction.key = "UMAP_dim30_")
  ob.list[[i]] <- RunUMAP(ob.list[[i]], dims = 1:50, reduction.name = "UMAP_dim50", reduction.key = "UMAP_dim50_")
  ob.list[[i]] <- FindNeighbors(ob.list[[i]], dims = 1:30)
  ob.list[[i]] <- FindClusters(ob.list[[i]], resolution = 1)
}

inflamed_mscs <- ob.list[[1]]
ctrl_lic_mscs <- ob.list[[2]]

# read in SB_MSC atlas
MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Seurat/MSCs_Final_091523.RDS")
MSCs <- RunUMAP(MSCs, dims = 1:50, reduction = "MSC_pca", reduction.name = "MSC_UMAP_dim50", return.model = TRUE) # rerun find umap to return the umap model

# Run AUCell and add Dataset Names

combined_inflamed_AUC <- runAUCell(object = inflamed_mscs, geneset_path = geneset_path, n_cells = 5000)
combined_inflamed_AUC$Dataset <- "Inflamed_MSC"

combined_ctrl_AUC <- runAUCell(object = ctrl_lic_mscs, geneset_path = geneset_path, n_cells = 5000)
combined_ctrl_AUC$Dataset <- "Control_MSC"

combined_SB_AUC <- runAUCell(object = MSCs, geneset_path = geneset_path, n_cells = 5000)
combined_SB_AUC$Dataset <- "SB_MSC"
Li_MSCs_AUC <- runAUCell(object = Li_MSCs, geneset_path = geneset_path, n_cells = 5000)
Li_MSCs_AUC$Dataset <- "Li_MSCs"



merged_inflamed <- merge(combined_SB_AUC, y = c(combined_ctrl_AUC, combined_inflamed_AUC, Li_MSCs_AUC))
VlnPlot(subset(merged_inflamed, predicted.MSC_refmap %in% c("Adipo-MSC", "Osteo-MSC", "Osteoblast", "Fibro-MSC", "OsteoFibro-MSC")) , group.by = "Dataset", features = "HALLMARK_INFLAMMATORY_RESPONSE", pt.size = 0)

# Supplemental Figure S2D - AUCel Inflammatory Response Signatures of Control, Inflamed, Li et al, and our atlas's MSCs (note RNAlo MSCs removed) ----
VlnPlot(merged_inflamed , group.by = "Dataset", features = "HALLMARK_INFLAMMATORY_RESPONSE", pt.size = 0) + geom_boxplot(outlier.shape = NA) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.3)) +
  stat_compare_means(comparisons = list(c("Control_MSC", "Inflamed_MSC") , c("Li_MSCs", "SB")),method = "t.test")

saveRDS(merged_inflamed, file = "/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/Revisions_Inflammation_Comparison/merged_inflamed.RDS")

# Calculate Statistics for S2C ----
t.test(combined_SB_AUC$HALLMARK_INFLAMMATORY_RESPONSE, combined_ctrl_AUC$HALLMARK_INFLAMMATORY_RESPONSE)
t.test(combined_SB_AUC$HALLMARK_INFLAMMATORY_RESPONSE, combined_inflamed_AUC$HALLMARK_INFLAMMATORY_RESPONSE)

t.test(combined_SB_AUC$HALLMARK_INFLAMMATORY_RESPONSE, Li_MSCs_AUC$HALLMARK_INFLAMMATORY_RESPONSE)

