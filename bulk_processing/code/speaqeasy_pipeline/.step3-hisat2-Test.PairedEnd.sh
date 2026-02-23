#!/bin/bash
#SBATCH --job-name=step3-hisat2-Test.PairedEnd
#SBATCH --chdir=.
#SBATCH --output=./logs/hisat2-Test.%a.txt
#SBATCH --error=./logs/hisat2-Test.%a.txt
#SBATCH --array=1-1%15
#SBATCH --cpus-per-task=8
#SBATCH --mem=25G

##SBATCH --dependency=afterok:pipeline_setup,step2-trim-Test.PairedEnd

#SBATCH --mail-type=ALL

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: jsundstr"
echo "Job id: 29137152"
echo "Job name: Step3-histat2.sh"
echo "Hostname: compute-153"
echo "Task id: "
echo "****"
echo "Sample id: "
echo "****"

# Directories
mkdir -p /users/jsundstr/hicks_home/bulk_processing//HISAT2_out/align_summaries
mkdir -p /users/jsundstr/hicks_home/bulk_processing//HISAT2_out/infer_strandness

if [ FALSE == "TRUE" ]
then
    mkdir -p /users/jsundstr/hicks_home/bulk_processing//HISAT2_out/unaligned
fi

## Load HISAT2 module

module load hisat2/2.2.1

## Locate file and ids
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk "NR==${SLURM_ARRAY_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk "NR==${SLURM_ARRAY_TASK_ID}")
fi
ID=$(cat /users/jsundstr/hicks_home/bulk_processing//fastq_test/samples.manifest | awk '{print $NF}' | awk "NR==${SLURM_ARRAY_TASK_ID}")

if [ -f /users/jsundstr/hicks_home/bulk_processing//trimmed_fq/${ID}_trimmed_forward_paired.fastq ] ; then
	## Trimmed, paired-end
	echo "HISAT2 alignment run on trimmed paired-end reads"
	FP=/users/jsundstr/hicks_home/bulk_processing//trimmed_fq/${ID}_trimmed_forward_paired.fastq
	FU=/users/jsundstr/hicks_home/bulk_processing//trimmed_fq/${ID}_trimmed_forward_unpaired.fastq
	RP=/users/jsundstr/hicks_home/bulk_processing//trimmed_fq/${ID}_trimmed_reverse_paired.fastq
	RU=/users/jsundstr/hicks_home/bulk_processing//trimmed_fq/${ID}_trimmed_reverse_unpaired.fastq
	
	hisat2 -p 8 	-x  -1 $FP -2 $RP -U ${FU},${RU} 	-S /users/jsundstr/hicks_home/bulk_processing//HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF --phred33      	2>/users/jsundstr/hicks_home/bulk_processing//HISAT2_out/align_summaries/${ID}_summary.txt
	
elif  [ -f /users/jsundstr/hicks_home/bulk_processing//trimmed_fq/${ID}_trimmed.fastq ] ; then
	## Trimmed, single-end
	echo "HISAT2 alignment run on trimmed single-end reads"
	hisat2 -p 8 	-x  -U /users/jsundstr/hicks_home/bulk_processing//trimmed_fq/${ID}_trimmed.fastq 	-S /users/jsundstr/hicks_home/bulk_processing//HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF --phred33 	2>/users/jsundstr/hicks_home/bulk_processing//HISAT2_out/align_summaries/${ID}_summary.txt

elif [ TRUE == "TRUE" ] ; then
	## Untrimmed, pair-end
	echo "HISAT2 alignment run on original untrimmed paired-end reads"
	hisat2 -p 8 	-x  -1 ${FILE1} -2 ${FILE2} 	-S /users/jsundstr/hicks_home/bulk_processing//HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF --phred33      	2>/users/jsundstr/hicks_home/bulk_processing//HISAT2_out/align_summaries/${ID}_summary.txt

else
	## Untrimmed, single-end
	echo "HISAT2 alignment run on original untrimmed single-end reads"
	hisat2 -p 8 	-x  -U ${FILE1} 	-S /users/jsundstr/hicks_home/bulk_processing//HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF --phred33 	2>/users/jsundstr/hicks_home/bulk_processing//HISAT2_out/align_summaries/${ID}_summary.txt
fi


###sam to bam
SAM=/users/jsundstr/hicks_home/bulk_processing//HISAT2_out/${ID}_hisat_out.sam
ORIGINALBAM=/users/jsundstr/hicks_home/bulk_processing//HISAT2_out/${ID}_accepted_hits.bam
SORTEDBAM=/users/jsundstr/hicks_home/bulk_processing//HISAT2_out/${ID}_accepted_hits.sorted

#filter unmapped segments
echo "**** Filtering unmapped segments ****"
date

module load samtools/1.18

samtools view -bh -F 4  > 
samtools sort -@ 8  
samtools index .bam

## Clean up
rm 
rm 

echo "**** Job ends ****"
date
