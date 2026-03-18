#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH --job-name=Clean_and_Subset
#SBATCH -o logs/Clean_and_Subset.log
#SBATCH -e logs/Clean_and_Subset.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3.x
Rscript Clean_and_Subset_hb_only.R

echo "**** Job ends ****"
date