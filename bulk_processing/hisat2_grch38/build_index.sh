#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=hisat2_index.log


cd /users/jsundstr/hicks_home/bulk_processing/hisat2_grch38

# load HISAT2 module 
module load hisat2/2.2.1     

# Build HISAT2 index
/jhpce/shared/libd/core/hisat2/2.2.1/hisat2-2.2.1/hisat2-build -p 16 GRCh38.primary_assembly.genome.fa hisat2_GRCh38primary

echo "Done"