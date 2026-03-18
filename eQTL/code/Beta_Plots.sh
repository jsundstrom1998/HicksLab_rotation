#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --job-name=08_beta_plots
#SBATCH -o logs/08_beta_plots.log
#SBATCH -e logs/08_beta_plots.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load conda_R/4.3.x
module load liftover/1.0
Rscript 08_beta_plots.R

echo "**** Job ends ****"
date