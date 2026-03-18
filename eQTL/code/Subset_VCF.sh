#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --job-name=Subset_VCF
#SBATCH -o logs/Subset_VCF.log
#SBATCH -e logs/Subset_VCF.log

big_vcf=/dcs04/hicks/data/data_proc/bulkDNAgenotype_habenula_human_pilot/v1/habenula_genotypes.vcf.gz
out_dir=/users/jsundstr/hicks_home/eQTL/

paired_variants=${out_dir}/DEA_paired_variants.txt
small_vcf=${out_dir}/paired_DEA_eQTL_SNPs.vcf

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

zcat $big_vcf \
    | grep -E "^#|$(paste -sd "|" $paired_variants)" \
    > $small_vcf
echo "**** Job ends ****"
date