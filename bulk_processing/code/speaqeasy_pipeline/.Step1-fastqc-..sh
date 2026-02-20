#!/bin/bash
#SBATCH --job-name=Step1-fastqc-.
#SBATCH --output=./logs/fastqc-.%a.txt
#SBATCH --error=./logs/fastqc-.%a.txt
#SBATCH --array=1-1%100
#SBATCH --mem=7G 

##SBATCH --dependency=afterok:pipeline_setup,step00-merge-.

#SBATCH --mail-type=ALL      
#SBATCH --chdir=.
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"
echo "Sample id: $(cat /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk '{print $NF}' | awk "NR==${SLURM_ARRAY_TASK_ID}")"
echo "****"

## Locate file and ids
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk "NR==${SLURM_ARRAY_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk "NR==${SLURM_ARRAY_TASK_ID}")
fi
ID=$(cat /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk '{print $NF}' | awk "NR==${SLURM_ARRAY_TASK_ID}")

mkdir -p /users/jsundstr/hicks_home/bulk_processing//FastQC/Untrimmed/${ID}

module load fastqc/0.11.8

if [ TRUE == "TRUE" ]
then 
    fastqc ${FILE1} ${FILE2} --outdir=/users/jsundstr/hicks_home/bulk_processing//FastQC/Untrimmed/${ID} --extract
else
    fastqc ${FILE1} --outdir=/users/jsundstr/hicks_home/bulk_processing//FastQC/Untrimmed/${ID} --extract
fi

echo "**** Job ends ****"
date
