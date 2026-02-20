#!/bin/bash

## Usage information:
# bash step1-fastqc.sh --help

# Define variables
TEMP=$(getopt -o x:p:l:h --long experiment:,prefix:,large:,help -n 'Step1-fastqc' -- "$@")
eval set -- "$TEMP"

LARGE="FALSE"

while true; do
    case "$1" in
        -x|--experiment)
            case "$2" in
                "") shift 2 ;;
                *) EXPERIMENT=$2 ; shift 2;;
            esac;;
        -p|--prefix)
            case "$2" in
                "") shift 2 ;;
                *) PREFIX=$2 ; shift 2;;
            esac;;
        -l|--large)
            case "$2" in
                "") LARGE="FALSE" ; shift 2;;
                *) LARGE=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step1-fastqc.sh -x -p -l (default:FALSE)\nLong options:\n  bash step1-fastqc.sh --experiment --prefix --large (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

SOFTWARE=/users/jsundstr/hicks_home/bulk_processing/code/speaqeasy_pipeline
MAINDIR=/users/jsundstr/hicks_home/bulk_processing/
SHORT="fastqc-${EXPERIMENT}"
sname="Step1-${SHORT}.${PREFIX}"

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=10G,h_vmem=14G,h_fsize=100G"
else
    MEM="mem_free=5G,h_vmem=7G,h_fsize=100G"
fi

if [ -f ".send_emails" ]
then
    EMAIL="e"
else
    EMAIL="a"
fi

if [ -f ".queue" ]
then
    SGEQUEUE="$(cat .queue),"
else
    SGEQUEUE=""
fi

if [ -f ".paired_end" ]
then
    PE="TRUE"
else
    PE="FALSE"
fi

# Directories
mkdir -p ${MAINDIR}/fastq_test/FastQC/Untrimmed
mkdir -p ${MAINDIR}/logs

# Construct shell files
FILELIST=${MAINDIR}/fastq_test/samples.manifest
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/code/speaqeasy_pipeline/.${sname}.sh <<EOF
#!/bin/bash
#SBATCH --job-name=${sname}
#SBATCH --output=./logs/${SHORT}.%a.txt
#SBATCH --error=./logs/${SHORT}.%a.txt
#SBATCH --array=1-${NUM}%100
#SBATCH --mem=7G 
             
##SBATCH --dependency=afterok:pipeline_setup,step00-merge-${EXPERIMENT}.${PREFIX}
#SBATCH --mail-type=ALL      
#SBATCH --chdir=.
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SLURM_ARRAY_TASK_ID}"
echo "****"
echo "Sample id: \$(cat ${MAINDIR}/fastq_test/samples.manifest | awk '{print \$NF}' | awk "NR==\${SLURM_ARRAY_TASK_ID}")"
echo "****"

## Locate file and ids
FILE1=\$(awk 'BEGIN {FS="\t"} {print \$1}' ${FILELIST} | awk "NR==\${SLURM_ARRAY_TASK_ID}")
if [ $PE == "TRUE" ] 
then
    FILE2=\$(awk 'BEGIN {FS="\t"} {print \$3}' ${FILELIST} | awk "NR==\${SLURM_ARRAY_TASK_ID}")
fi
ID=\$(cat ${FILELIST} | awk '{print \$NF}' | awk "NR==\${SLURM_ARRAY_TASK_ID}")

mkdir -p ${MAINDIR}/FastQC/Untrimmed/\${ID}

module load fastqc/0.11.8

if [ $PE == "TRUE" ]
then 
    fastqc \
\${FILE1} \${FILE2} \
--outdir=${MAINDIR}/FastQC/Untrimmed/\${ID} --extract
else
    fastqc \${FILE1} \
--outdir=${MAINDIR}/FastQC/Untrimmed/\${ID} --extract
fi

echo "**** Job ends ****"
date
EOF

call="sbatch .${sname}.sh"
echo $call
$call

