#!/bin/bash
#SBATCH -n 6
#SBATCH --job-name=S11
#SBATCH --mem=50G
#SBATCH -t 30-0
#SBATCH -o slurm_out/slurm-%j.out

/mnt/isilon/tan_lab/xuj5/software/cellranger4-0/cellranger-4.0.0/cellranger count --fastqs=/mnt/isilon/tan_lab/bandyopads/SB15_normal_huMSC_scRNA_deJong/data/S11/ --id=S11 --transcriptome=/mnt/isilon/tan_lab/chenc6/Tools/SingleCellAnnotation/GRCh38

