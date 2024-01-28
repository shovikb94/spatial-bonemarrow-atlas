# spatial-bonemarrow-atlas
This repository contains the code used to generate figures and perform analysis for the manuscript titled "Transcriptomic and Spatial Proteomic Profiling Reveals the Cellular Composition and Spatial Organization of the Human Bone Marrow Microenvironment"

**Figure Manifest**
____________________________________
**General scRNA-Seq Analysis**
Figure 1B-D, S1C-G,Figure 2A-C, Figure 2G-H, Figure S3A-B, Figure S3E-F , Supplemental Table S2, Figure 3A, Figure 5E, Figure S7B- scRNA_Seq_Analysis/Step1_AtlasGeneration_and_Figures/SB_Manuscript_scRNA_Analysis.R

**scRNA-Seq Inflammation Analysis**
Figure S2A-D
scRNA_Seq_Analysis/Step1.2_MapExternalDatasets/Inflammation_Comparison_ForManuscript.R

**CytoTRACE Analysis**
Figure 2D – 
scRNA_Seq_Analysis/Step1.1_Run_CytoTRACE/CytoTrace_BM.R

**External Dataset Mapping**
Figure 2I, Supplemental Figure S3G scRNA_Seq_Analysis/Step1.2_MapExternalDatasets/Map_MSC_Datasets_Manuscript.R
scRNA_Seq_Analysis/Step1.2_MapExternalDatasets/Integrate_Mapped_Datasets_ForManuscript_v3.R

Figure S3H - 
scRNA_Seq_Analysis/Step1.2_MapExternalDatasets/Bandyopadhyay_Jardine_integration_ForManuscript.R

**CellChat Analysis**
Figure S3G, S4E, S8B - scRNA_Seq_Analysis/Step2_CellChat/CellChat_ForMaunscript.R

**CODEX Cell Type Annotation**
Figure 4C-4D - CODEX_Analysis/Step3_AnnotateCells/NBM_CODEX_Atlas_Combined_Step3_AnnotateCells.R

**CODEX MSC Subset Manual Annotation and Distance Analysis**
Figure 4E - CODEX_Analysis/Step8_Distance_Analysis/MSC_Distance_Analysis_VlnPlots.R
CODEX_Analysis/Step8_Distance_Analysis/Osteo_Fibro_MSC_Get_Bone_Overlap.ipynb

**CODEX NBM Neighborhood Analysis**
Figure 5A, 5D, S7D-E - CODEX_Analysis/Step4_NBM_NeighborhoodAnalysis/NeighborhoodAnalysis_Step4_FINAL.R
CODEX_Analysis/Step4_NBM_NeighborhoodAnalysis/Combined_NeighborhoodAnalysis.ipynb

**CODEX Structure Distance Analysis**
Figure 6A-E, S9B, S10G - CODEX_Analysis/Step8_Distance_Analysis/Plot_Rank_Data_ForManuscript.R

**CODEX NSM Annotation**
Figure 7A, S10A, also needed for 7B-C and Supplemental Figure S11
CODEX_Analysis/Step5_NSM_ReferenceMapping/ReferenceMap_NSM_ForManuscript.R

**CODEX AML Blast Mapping**
Figure 7B
CODEX_Analysis/Step6_AML_ReferenceMapping/RefMap_AML_ForManuscript_Final.R
CODEX_Analysis/Step6_AML_ReferenceMapping/SpatialJoin_AML.R

**CODEX AML Analysis**
Figure 7C-D, 7F-G, S10E, S10G-J - CODEX_Analysis/Step7_NeighborhoodAnalysis_AMLandNSM/AML_NSM_neighborhood_analysis_ForManuscript.R

**CODEX NSM Neighborhood Analysis**
Figure S11A-B - CODEX_Analysis/Step5_NSM_ReferenceMapping/NSM_NBM_Correlation.R
CODEX_Analysis/Step5_NSM_ReferenceMapping/NSM_Only_NeighborhoodAnalysis-.ipynb

**Integrated Analysis**
Figure 5F
Integrated_Analysis/Integrate_CellChat_and_NeighborhoodAnalysis_ForManuscript_v2.R

Figure S5F-G
Integrated_Analysis/RNA_Protein_Correlation_Median_ForManuscript.R

