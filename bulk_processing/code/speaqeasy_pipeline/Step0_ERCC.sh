#!/bin/bash
#Prevent silent errors
set -euo pipefail

#Define variables 
TEMP=$(getopt -o x:p:c:l:s:h \
  --long experiment:,prefix:,cores:,large:,stranded:,help \
  -n 'step0-ercc' -- "$@")
eval set -- "$TEMP"

#Default parameters
LARGE="FALSE"
CORES=8
STRANDED="FALSE"

#Assign command line arguments to variables
while true; do
    case "$1" in
        -x|--experiment) EXPERIMENT=$2; shift 2 ;;
        -p|--prefix) PREFIX=$2; shift 2 ;;
        -c|--cores) CORES=${2:-8}; shift 2 ;;
        -l|--large) LARGE=${2:-FALSE}; shift 2 ;;
        -s|--stranded) STRANDED=${2:-FALSE}; shift 2 ;;
        -h|--help)
            echo "Usage: step0-ercc-slurm.sh --experiment --prefix [--cores] [--large] [--stranded]"
            exit 0 ;;
        --) shift; break ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

#Defining more variables
SOFTWARE=/dcs04/hicks/data/jsundstr/bulk_processing
MAINDIR=/dcs04/hicks/data/jsundstr/bulk_processing/code
SHORT="ercc-${EXPERIMENT}"
JOBNAME="step0-${SHORT}.${PREFIX}" 

FILELIST=/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/samples.manifest
NUM=$(awk '{print $NF}' ${FILELIST} | uniq | wc -l)

ERCC_JOBID=""

mkdir -p logs

#Define memory
if [[ $LARGE == "TRUE" ]]; then
    MEM="10G"                    
else
    MEM="5G"
fi

#Paired-end flag 
if [[ -f ".paired_end" ]]; then
    PE="TRUE"
else
    PE="FALSE"
fi

#Strandedness flags
STRANDOPTION=""
case "${STRANDED}" in  
    FALSE) ;;
    forward) STRANDOPTION="--fr-stranded" ;;
    reverse) STRANDOPTION="--rf-stranded" ;;
    *) echo "Invalid stranded option: ${STRANDED}"; exit 1 ;;
esac

#Create ERCC array job and submit to SLURM
cat > ${MAINDIR}/.${JOBNAME}.sh <<EOF
#!/bin/bash

#Define SLURM parameters
#SBATCH --job-name=${JOBNAME}
#SBATCH --cpus-per-task=${CORES}
#SBATCH --mem=${MEM}
#SBATCH --array=1-${NUM}%20
#SBATCH --output=logs/${SHORT}.%a.out
#SBATCH --error=logs/${SHORT}.%a.err

echo "**** Job starts ****"
date
echo "Node: \$(hostname)"
echo "Task ID: \$SLURM_ARRAY_TASK_ID"

#Get file names and IDs from manifest
ID=\$(awk '{print \$NF}' ${FILELIST} | awk "NR==\$SLURM_ARRAY_TASK_ID")
FILE1=\$(awk '{print \$1}' ${FILELIST} | awk "NR==\$SLURM_ARRAY_TASK_ID")

#Use second file for paired-end data
if [[ "${PE}" == "TRUE" ]]; then
    FILE2=\$(awk '{print \$3}' ${FILELIST} | awk "NR==\$SLURM_ARRAY_TASK_ID")
fi

#Create output directory
mkdir -p ${MAINDIR}/Ercc/\${ID}

#Run Kallisto based on paried-end or single-end options

module load kallisto

if [[ "${PE}" == "TRUE" ]]; then
    kallisto quant \\
        -i /dcs04/hicks/data/jsundstr/bulk_processing/ERCC/ERCC92.idx \\
        -o ${MAINDIR}/Ercc/\${ID} \\
        -t ${CORES} ${STRANDOPTION} \\
        \${FILE1} \${FILE2}
else
    kallisto quant \\
        -i /dcs04/hicks/data/jsundstr/bulk_processing/ERCC/ERCC92.idx \\
        -o ${MAINDIR}/Ercc/\${ID} \\
        -t ${CORES} --single ${STRANDOPTION} \\
        \${FILE1}
fi
echo "**** Job ends ****"
date
EOF

#Submit ERCC job w/ dependency on merge job
if squeue -h -o "%A" | grep -q "pipeline_setup" && squeue -h -o "%A" | grep -q "step00-merge-${EXPERIMENT}.${PREFIX}"; then
    ERCC_JOBID=$(sbatch --dependency=afterok:pipeline_setup:step00-merge-${EXPERIMENT}.${PREFIX} .${JOBNAME}.sh | awk '{print $4}')
else
    ERCC_JOBID=$(sbatch .${JOBNAME}.sh | awk '{print $4}')
fi

echo "ERCC job submitted with JobID: ${ERCC_JOBID}"

