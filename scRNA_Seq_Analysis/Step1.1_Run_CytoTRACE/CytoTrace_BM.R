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

