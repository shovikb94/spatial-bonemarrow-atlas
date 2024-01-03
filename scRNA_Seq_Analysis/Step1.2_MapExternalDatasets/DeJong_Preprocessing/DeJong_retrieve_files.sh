#!/bin/bash
#SBATCH --job-name=SB15_retrievefastq
#SBATCH -n 2
#SBATCH --mem=50G
#SBATCH -o slurm_out/%x-%A_%a.out
#SBATCH -t 30-0
mkdir S9 
cd S9
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5941_N_S9_L001_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5941_N_S9_L001_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5941_N_S9_L001_I1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5941_N_S9_L002_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5941_N_S9_L002_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5941_N_S9_L002_I1_001.fastq.gz
cd ..
mkdir S10
cd S10
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5950_N_S10_L001_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5950_N_S10_L001_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5950_N_S10_L001_I1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5950_N_S10_L002_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5950_N_S10_L002_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/2020_5950_N_S10_L002_I1_001.fastq.gz
cd ..
mkdir S11
cd S11
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D03_N_S11_L001_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D03_N_S11_L001_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D03_N_S11_L001_I1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D03_N_S11_L002_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D03_N_S11_L002_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D03_N_S11_L002_I1_001.fastq.gz
cd ..
mkdir S12
cd S12
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D04_N_S12_L001_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D04_N_S12_L001_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D04_N_S12_L001_I1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D04_N_S12_L002_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D04_N_S12_L002_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D04_N_S12_L002_I1_001.fastq.gz
cd ..
mkdir S13
cd S13
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D05_N_S13_L001_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D05_N_S13_L001_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D05_N_S13_L001_I1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D05_N_S13_L002_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D05_N_S13_L002_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-9139/HBM_D05_N_S13_L002_I1_001.fastq.gz

