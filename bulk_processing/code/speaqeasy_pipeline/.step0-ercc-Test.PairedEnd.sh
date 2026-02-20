#!/bin/bash

#Define SLURM parameters
#SBATCH --job-name=step0-ercc-Test.PairedEnd
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --array=1-1%20
#SBATCH --output=logs/ercc-Test.%a.out
#SBATCH --error=logs/ercc-Test.%a.err

echo "**** Job starts ****"
date
echo "Node: $(hostname)"
echo "Task ID: $SLURM_ARRAY_TASK_ID"

#Get file names and IDs from manifest
ID=$(awk '{print $NF}' /dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/samples.manifest | awk "NR==$SLURM_ARRAY_TASK_ID")
FILE1=$(awk '{print $1}' /dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/samples.manifest | awk "NR==$SLURM_ARRAY_TASK_ID")

#Use second file for paired-end data
if [[ "TRUE" == "TRUE" ]]; then
    FILE2=$(awk '{print $3}' /dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/samples.manifest | awk "NR==$SLURM_ARRAY_TASK_ID")
fi

#Create output directory
mkdir -p /dcs04/hicks/data/jsundstr/bulk_processing/code/speaqeasy_pipeline/Ercc/${ID}

#Run Kallisto based on paried-end or single-end options

module load kallisto

if [[ "TRUE" == "TRUE" ]]; then
    kallisto quant \
        -i /dcs04/hicks/data/jsundstr/bulk_processing/ERCC/ERCC92.idx \
        -o /dcs04/hicks/data/jsundstr/bulk_processing/code/speaqeasy_pipeline/Ercc/${ID} \
        -t 8 --rf-stranded \
        ${FILE1} ${FILE2}
else
    kallisto quant \
        -i /dcs04/hicks/data/jsundstr/bulk_processing/ERCC/ERCC92.idx \
        -o /dcs04/hicks/data/jsundstr/bulk_processing/code/speaqeasy_pipeline/Ercc/${ID} \
        -t 8 --single --rf-stranded \
        ${FILE1}
fi
echo "**** Job ends ****"
date
