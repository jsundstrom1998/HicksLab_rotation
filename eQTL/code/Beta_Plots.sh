#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --job-name=Beta_Plots
#SBATCH -o logs/Beta_Plots.log
#SBATCH -e logs/Beta_Plots.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.5
module load liftover/1.0
Rscript Beta_Plots.R

echo "**** Job ends ****"
date