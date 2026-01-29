#!/bin/bash
set -euo pipefail

## Define variables 
TEMP=$(getopt -o x:p:c:l:s:h \
  --long experiment:,prefix:,cores:,large:,stranded:,help \
  -n 'step0-ercc' -- "$@")
eval set -- "$TEMP"

LARGE="FALSE"
CORES=8
STRANDED="FALSE"

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

## Setup 
SOFTWARE=/dcs04/hicks/data/jsundstr/bulk_processing
MAINDIR=/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test
SHORT="ercc-${EXPERIMENT}"
JOBNAME="step0-${SHORT}.${PREFIX}" 

FILELIST=${MAINDIR}/samples.manifest
NUM=$(awk '{print $NF}' ${FILELIST} | uniq | wc -l)


mkdir -p logs

## Define memory
if [[ $LARGE == "TRUE" ]]; then
    MEM="10G"                    
else
    MEM="5G"
fi

## Paired-end flag 
if [[ -f ".paired_end" ]]; then
    PE="TRUE"
else
    PE="FALSE"
fi

## Strandedness flags
STRANDOPTION=""
if [[ $STRANDED == "forward" ]]; then
    STRANDOPTION="--fr-stranded"
elif [[ $STRANDED == "reverse" ]]; then
    STRANDOPTION="--rf-stranded"
elif [[ $STRANDED != "FALSE" ]]; then
    echo "Invalid stranded option"
    exit 1
fi

## Write ERCC array job
cat > .${JOBNAME}.slurm <<EOF
#!/bin/bash
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

ID=\$(awk '{print \$NF}' ${FILELIST} | awk "NR==\$SLURM_ARRAY_TASK_ID")
FILE1=\$(awk '{print \$1}' ${FILELIST} | awk "NR==\$SLURM_ARRAY_TASK_ID")

if [[ "${PE}" == "TRUE" ]]; then
    FILE2=\$(awk '{print \$3}' ${FILELIST} | awk "NR==\$SLURM_ARRAY_TASK_ID")
fi

mkdir -p ${MAINDIR}/Ercc/\${ID}

if [[ "${PE}" == "TRUE" ]]; then
    ${SOFTWARE}/kallisto quant \\
        -i /dcs04/hicks/data/jsundstr/bulk_processing/ERCC/ERCC92.idx \\
        -o ${MAINDIR}/Ercc/\${ID} \\
        -t ${CORES} ${STRANDOPTION} \\
        \${FILE1} \${FILE2}
else
    ${SOFTWARE}/kallisto quant \\
        -i /dcs04/hicks/data/jsundstr/bulk_processing/ERCC/ERCC92.idx \\
        -o ${MAINDIR}/Ercc/\${ID} \\
        -t ${CORES} --single ${STRANDOPTION} \\
        \${FILE1}
fi
echo "**** Job ends ****"
date
EOF

## Submit ERCC job WITH dependency
ERCC_JOBID=$(sbatch --dependency=afterok:${MERGE_JOBID} .${JOBNAME}.slurm | awk '{print $4}')

echo "ERCC job submitted with JobID: ${ERCC_JOBID}"
echo "Dependency: afterok:${MERGE_JOBID}"
