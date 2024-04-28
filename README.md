# spatial-bonemarrow-atlas
This repository contains the code used to generate figures and perform analysis for the manuscript titled "Transcriptomic and Spatial Proteomic Profiling Reveals the Cellular Composition and Spatial Organization of the Human Bone Marrow Microenvironment". The data from this manuscript is interactively browsable through our [Vitessce-based platform](https://cscb.research.chop.edu/index.php/bm-data).

**Data Availability**
_______________________________________
Raw and processed scRNA-Seq data can be accessed through [GSE253355](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253355). 
Raw and processed CODEX data can be accessed through our repository on [FigShare](https://doi.org/10.25452/figshare.plus.c.7174914). 

**Figure Manifest**
____________________________________
**General scRNA-Seq Analysis**
Figure 1B-D, S1C-F,Figure 2A-C, Figure 2G-H, Figure S3A-C, Figure S3E , Supplemental Table S2, Figure 3A, Figure 5E, Figure S7B- scRNA_Seq_Analysis/Step1_AtlasGeneration_and_Figures/SB_Manuscript_scRNA_Analysis.R

**scRNA-Seq Inflammation Analysis**
Figure S2A-D
scRNA_Seq_Analysis/Step1.2_MapExternalDatasets/Inflammation_Comparison_ForManuscript.R

**CytoTRACE Analysis**
Figure 2D – 
scRNA_Seq_Analysis/Step1.1_Run_CytoTRACE/CytoTrace_BM.R

**External Dataset Mapping**
Figure 2I, Supplemental Figure S4A scRNA_Seq_Analysis/Step1.2_MapExternalDatasets/Map_MSC_Datasets_Manuscript.R
scRNA_Seq_Analysis/Step1.2_MapExternalDatasets/Integrate_Mapped_Datasets_ForManuscript_v3.R

Figure S4B - 
scRNA_Seq_Analysis/Step1.2_MapExternalDatasets/Bandyopadhyay_Jardine_integration_ForManuscript.R

**CellChat Analysis**
Figure S4A, S4E, S8B - scRNA_Seq_Analysis/Step2_CellChat/CellChat_ForMaunscript.R

**CODEX Cell Type Annotation**
Figure 4C-4D - CODEX_Analysis/Step3_AnnotateCells/NBM_CODEX_Atlas_Combined_Step3_AnnotateCells.R

**CODEX MSC Subset Manual Annotation and Distance Analysis**
Figure S6E - CODEX_Analysis/Step1_CreateSeuratObjs/Fibro_Osteo_MSC_Create_Seurat.R
Figure 4E - CODEX_Analysis/Step8_Distance_Analysis/MSC_Distance_Analysis_VlnPlots.R

**CODEX NBM Neighborhood Analysis**
Figure 5A, 5D, S7D-E - CODEX_Analysis/Step4_NBM_NeighborhoodAnalysis/NeighborhoodAnalysis_Step4_FINAL.R
CODEX_Analysis/Step4_NBM_NeighborhoodAnalysis/Combined_NeighborhoodAnalysis.ipynb
Neighborhood Enrichment Stats - CODEX_Analysis/Step4_NBM_NeighborhoodAnalysis/Hypergeometric_Statistics_NBM.R

**CODEX Structure Distance Analysis**
Figure 6A-E, S9B, S11D - CODEX_Analysis/Step8_Distance_Analysis/Plot_Rank_Data_ForManuscript.R

**CODEX NSM Annotation**
Figure 7A, S10A, S10D-E, also needed for 7B-C
CODEX_Analysis/Step5_NSM_ReferenceMapping/ReferenceMap_NSM_ForManuscript.R

**CODEX AML Blast Mapping**
Figure 7B
CODEX_Analysis/Step6_AML_ReferenceMapping/RefMap_AML_ForManuscript_Final.R
CODEX_Analysis/Step6_AML_ReferenceMapping/SpatialJoin_AML.R

**CODEX AML Analysis**
Figure 7C-D, 7F-G, S11B,S11D, S11E-G - CODEX_Analysis/Step7_NeighborhoodAnalysis_AML/AML_NSM_neighborhood_analysis_ForManuscript.R
Neighborhood Enrichment Stats - CODEX_Analysis/Step7_NeighborhoodAnalysis_AML/Hypergeometric_Statistics_AML.R

**CODEX NSM Neighborhood Analysis**
Figure S10D-E - CODEX_Analysis/Step5_NSM_ReferenceMapping/NSM_NBM_Correlation.R
CODEX_Analysis/Step5_NSM_ReferenceMapping/NSM_Only_NeighborhoodAnalysis-.ipynb
Neighborhood Enrichment Stats - CODEX_Analysis/Step5_NSM_ReferenceMapping/Hypergeometric_Statistics_NSM.R

**Integrated Analysis**
Figure 5F
Integrated_Analysis/Integrate_CellChat_and_NeighborhoodAnalysis_ForManuscript_v2.R

Figure S5F-G
Integrated_Analysis/RNA_Protein_Correlation_Median_ForManuscript.R

**Ligand-Receptor CODEX Analysis**
Figure S8C
CODEX_Analysis/Step9_LR_Analysis/LR_Distance_Analysis_ForManuscript.R

