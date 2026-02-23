#!/bin/bash
#SBATCH --job-name=Step2-trim-Test.PairedEnd
#SBATCH --chdir=.
#SBATCH --output=./logs/trim-Test.%a.txt
#SBATCH --error=./logs/trim-Test.%a.txt
#SBATCH --array=1-1%5
#SBATCH --cpus-per-task=8
#SBATCH --mem=7G     

##SBATCH --dependency=afterok:pipeline_setup,step1-fastqc-Test.PairedEnd

#SBATCH --mail-type=ALL          

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"
echo "Sample id: $(cat /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk '{print $NF}' | awk "NR==${SLURM_ARRAY_TASK_ID}")"
echo "****"

## Locate file and idsf
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk "NR==${SLURM_ARRAY_TASK_ID}")
FILEBASE1=$(basename ${FILE1} | sed 's/.fq.gz//; s/.fq//; s/.fastq.gz//; s/.fastq//')
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk "NR==${SLURM_ARRAY_TASK_ID}")
    FILEBASE2=$(basename ${FILE2} | sed 's/.fq.gz//; s/.fq//; s/.fastq.gz//; s/.fastq//')
fi
ID=$(cat /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk '{print $NF}' | awk "NR==${SLURM_ARRAY_TASK_ID}")

if [ TRUE == "TRUE" ] ; then 
	REPORT1=/users/jsundstr/hicks_home/bulk_processing//FastQC/Untrimmed/${ID}/${FILEBASE1}_fastqc/summary.txt
	REPORT2=/users/jsundstr/hicks_home/bulk_processing//FastQC/Untrimmed/${ID}/${FILEBASE2}_fastqc/summary.txt
	RESULT1=$(grep "Adapter Content" $REPORT1 | cut -c1-4)
	RESULT2=$(grep "Adapter Content" $REPORT2 | cut -c1-4)

	if [[ $RESULT1 == "FAIL" || $RESULT2 == "FAIL" ]] ; then
		## trim, rerun fastQC
		echo "End 1 adapters: $RESULT1"
		echo "End 2 adapters: $RESULT2"
		echo "Trimming will occur."
		
		mkdir -p /users/jsundstr/hicks_home/bulk_processing//FastQC/trimmed_fq
		FP=/users/jsundstr/hicks_home/bulk_processing//FastQC/trimmed_fq/${ID}_trimmed_forward_paired.fastq
		FU=/users/jsundstr/hicks_home/bulk_processing//FastQC/trimmed_fq/${ID}_trimmed_forward_unpaired.fastq
		RP=/users/jsundstr/hicks_home/bulk_processing//FastQC/trimmed_fq/${ID}_trimmed_reverse_paired.fastq
		RU=/users/jsundstr/hicks_home/bulk_processing//FastQC/trimmed_fq/${ID}_trimmed_reverse_unpaired.fastq
		
		## trim adapters
        module load trimmomatic/0.39
        trimmomatic PE -threads 8 -phred33             ${FILE1} ${FILE2} $FP $FU $RP $RU             ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:1             LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

		## rerun fastqc

        module load fastqc/0.11.8

		mkdir -p /users/jsundstr/hicks_home/bulk_processing//FastQC/Trimmed/${ID}
		fastqc 
		$FP $FU $RP $RU 		--outdir=/users/jsundstr/hicks_home/bulk_processing//FastQC/Trimmed/${ID} --extract
	else
		echo "No trimming required!"
	fi

else
	## reads are single-end
	REPORT1=/users/jsundstr/hicks_home/bulk_processing//FastQC/Untrimmed/${ID}/${FILEBASE1}_fastqc/summary.txt
	RESULT1=$(grep "Adapter Content" $REPORT1 | cut -c1-4)

	if [[ $RESULT1 == "FAIL" ]] ; then
		## trim, rerun fastQC
		echo "Adapters: $RESULT1"
		echo "Trimming will occur."
		
		mkdir -p /users/jsundstr/hicks_home/bulk_processing//FastQC/trimmed_fq
		OUT=/users/jsundstr/hicks_home/bulk_processing//FastQC/trimmed_fq/${ID}_trimmed.fastq
		
		## trim adapters
        module load trimmomatic/0.39
        trimmomatic SE -threads 8 -phred33             ${FILE1} $OUT             ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:1             LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		## rerun fastqc
		mkdir -p /users/jsundstr/hicks_home/bulk_processing//FastQC/Trimmed/${ID}
		fastqc 		--outdir=/users/jsundstr/hicks_home/bulk_processing//FastQC/Trimmed/${ID} --extract
	else
		echo "No trimming required!"
	fi
fi

	
echo "**** Job ends ****"
date
