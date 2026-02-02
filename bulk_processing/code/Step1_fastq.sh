#!/bin/bash
#SBATCH --job-name=step1-fastqc-${EXPERIMENT}.${PREFIX}
#SBATCH --chdir=.
#SBATCH --output=logs/fastqc-${EXPERIMENT}.%a.txt
#SBATCH --error=logs/fastqc-${EXPERIMENT}.%a.txt
#SBATCH --array=1-${NUM}%100
#SBATCH --dependency=afterok:pipeline_setup,step00-merge-${EXPERIMENT}.${PREFIX}
#SBATCH --time=4:00:00
#SBATCH --mem=${MEM_SLURM}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${USER}@jh.edu

echo "**** Job starts ****"
date
echo "User: ${USER}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Hostname: ${HOSTNAME}"
echo "****"

