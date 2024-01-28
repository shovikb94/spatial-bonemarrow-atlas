# Run cytotrace from MSC_Final object

library(readr)
library(CytoTRACE)
library(tibble)
library(Seurat)


# Run CytoTrace
MSCs <- readRDS("/mnt/isilon/tan_lab/bandyopads/SB66_scRNASeq_S2/CytoTRACE/input/MSCs_final.RDS")
# Prepare CytoTRACE input
cti <- GetAssayData(MSCs, slot = 'data', assay = 'RNA')
cti <- as.data.frame(cti) # cytotrace error does not mention but need to convert matrix to DF to run
# datasets <- list(marrow_10x_expr, marrow_plate_expr)

#No need to subsample
cytotrace <- CytoTRACE(cti, ncores = 8) 

MSCs <- AddMetaData(MSCs, cytotrace$CytoTRACE, col.name = "CytoTRACE_score")

saveRDS(MSCs, "/mnt/isilon/tan_lab/sussmanj/Temp/CytoTrace_BM/MSCs_CytoTRACE_computed.RDS")
write.csv(cytotrace$CytoTRACE, "/mnt/isilon/tan_lab/sussmanj/Temp/CytoTrace_BM/MSCs_CytoTRACE_scores.csv")

MSCs$original.seurat.clusters <- Idents(MSCs)
ctype.sele = subset(MSCs@meta.data, select = c('original.seurat.clusters'))
phe.vec = as.character(ctype.sele$original.seurat.clusters)
names(phe.vec) = rownames(ctype.sele)

plotCytoTRACE(cytotrace, phe.vec, 
              emb = MSCs@reductions$MSC_UMAP_dim50@cell.embeddings)

plotCytoGenes(cytotrace)


# Figure 2D - CytoTRACE Analysis  ----
## See Cytotrace script, run by Jon Sussman 
MSCs_CT <- readRDS("~/Documents/NBM_Microenvironment/NBM_Atlas_scRNA/Final_scRNA_Analysis/CytoTRACE/MSCs_CytoTRACE_computed.RDS")
p1 <- FeaturePlot(MSCs_CT, features = "CytoTRACE_score", reduction = "MSC_UMAP_dim50", pt.size = 1,coord.fixed = TRUE) + NoAxes() + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p1_raster <- rasterize(p1, dpi = 300)
ggsave(p1_raster, device = "pdf", filename = "~/Documents/Manuscripts/NBM_Atlas/Figures/Figure2/panel/CytoTRACE_Score_FeaturePlot.pdf", units = "in", width = 5, height = 5)




