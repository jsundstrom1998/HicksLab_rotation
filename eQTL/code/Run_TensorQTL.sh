#!/bin/bash

#SBATCH -c 1
##SBATCH --gpus=1
#SBATCH --mem=40G
#SBATCH --job-name=Run_TensorQTL
#SBATCH -o /dev/null
#SBATCH -e /dev/null

run_mode=interaction
interaction_cov=tot_Thal

if [[ "$run_mode" == "interaction" ]]; then
    log_path="logs/run_tensorqtl_${run_mode}_${interaction_cov}.log"
else
    log_path="logs/run_tensorqtl_${run_mode}.log"
fi

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

module load tensorqtl/1.0.8
python Run_TensorQTL.py $run_mode $interaction_cov

echo "**** Job ends ****"
date
} > $log_path 2>&1