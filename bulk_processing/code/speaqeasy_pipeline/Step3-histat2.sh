#!/bin/bash

## Usage information:
# bash step3-hisat2.sh --help

# Define variables
TEMP=$(getopt -o x:p:i:b:l:s:u:h --long experiment:,prefix:,index:,bed:,large:,stranded:,unaligned:,help -n 'step3-hisat2' -- "$@")
eval set -- "$TEMP"

LARGE="FALSE"
STRANDED="FALSE"
UNALIGNED="FALSE"

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
        -i|--index)
            case "$2" in
                "") shift 2 ;;
                *) HISATIDX=$2 ; shift 2;;
            esac;;
        -b|--bed)
            case "$2" in
                "") shift 2 ;;
                *) BED=$2 ; shift 2;;
            esac;;
        -l|--large)
            case "$2" in
                "") LARGE="FALSE" ; shift 2;;
                *) LARGE=$2; shift 2;;
            esac ;;
        -s|--stranded)
            case "$2" in
                "") STRANDED="FALSE" ; shift 2;;
                *) STRANDED=$2; shift 2;;
            esac ;;
        -u|--unaligned)
            case "$2" in
                "") UNALIGNED="FALSE" ; shift 2;;
                *) UNALIGNED=$2; shift 2;;
            esac ;;
        -h|--help)
            echo -e "Usage:\nShort options:\n  bash step3-hisat2.sh -x -p -i -b -l (default:FALSE) -s (default:FALSE) -u (default:FALSE)\nLong options:\n  bash step3-hisat2.sh --experiment --prefix --index --bed --large (default:FALSE) --stranded (default:FALSE) --unaligned (default:FALSE)"; exit 0; shift ;;
            --) shift; break ;;
        *) echo "Incorrect options!"; exit 1;;
    esac
done

#SOFTWARE=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software
MAINDIR=/users/jsundstr/hicks_home/bulk_processing/
SHORT="hisat2-${EXPERIMENT}"
sname="step3-${SHORT}.${PREFIX}"
CORES=8

if [[ $LARGE == "TRUE" ]]
then
    MEM="mem_free=12G,h_vmem=14G,h_fsize=200G"
else
    MEM="mem_free=5G,h_vmem=6G,h_fsize=200G"
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
    if [ ${STRANDED} == "FALSE" ]
    then
        STRANDOPTION=""
    elif [ ${STRANDED} == "forward" ]
    then
        STRANDOPTION="--rna-strandness FR"
    elif [ ${STRANDED} == "reverse" ]
    then
        STRANDOPTION="--rna-strandness RF"
    else
        echo "The option --stranded has to either be 'FALSE', 'forward' or 'reverse'."
        exit 1
    fi
else
    PE="FALSE"
    if [ ${STRANDED} == "FALSE" ]
    then
        STRANDOPTION=""
    elif [ ${STRANDED} == "forward" ]
    then
        STRANDOPTION="--rna-strandness F"
    elif [ ${STRANDED} == "reverse" ]
    then
        STRANDOPTION="--rna-strandness R"
    else
        echo "The option --stranded has to either be 'FALSE', 'forward' or 'reverse'."
        exit 1
    fi
fi

if [ ${UNALIGNED} == "FALSE" ]
then
    UNALIGNEDOPT=""
elif [ ${UNALIGNED} == "TRUE" ]
then
    mkdir -p ${MAINDIR}/HISAT2_out/unaligned
    UNALIGNEDOPT="--un-conc ${MAINDIR}/HISAT2_out/unaligned/\${ID}.fastq"
else
    echo "The option --unaligned has to either be 'FALSE' or 'TRUE'"
    exit 1
fi
    

# Construct shell files
FILELIST=${MAINDIR}/fastq_test/samples.manifest
NUM=$(cat $FILELIST | awk '{print $NF}' | uniq | wc -l)
echo "Creating script ${sname}"

cat > ${MAINDIR}/code/speaqeasy_pipeline/.${sname}.sh <<EOF
#!/bin/bash
#SBATCH --job-name=${sname}
#SBATCH --chdir=.
#SBATCH --output=./logs/${SHORT}.%a.txt
#SBATCH --error=./logs/${SHORT}.%a.txt
#SBATCH --array=1-${NUM}%15
#SBATCH --cpus-per-task=${CORES}
#SBATCH --mem=25G

##SBATCH --dependency=afterok:pipeline_setup,step2-trim-${EXPERIMENT}.${PREFIX}

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
echo "Sample id: $(cat ${MAINDIR}/fastq_test/samples.manifest | awk '{print \$NF}' | awk "NR==\${SLURM_ARRAY_TASK_ID}")"
echo "****"

# Directories
mkdir -p ${MAINDIR}/HISAT2_out/align_summaries
mkdir -p ${MAINDIR}/HISAT2_out/infer_strandness

if [ ${UNALIGNED} == "TRUE" ]
then
    mkdir -p ${MAINDIR}/HISAT2_out/unaligned
fi

## Load HISAT2 module

module load hisat2/2.2.1

## Locate file and ids
FILE1=\$(awk 'BEGIN {FS="\t"} {print \$1}' ${FILELIST} | awk "NR==\${SLURM_ARRAY_TASK_ID}")
if [ $PE == "TRUE" ] 
then
    FILE2=\$(awk 'BEGIN {FS="\t"} {print \$3}' ${FILELIST} | awk "NR==\${SLURM_ARRAY_TASK_ID}")
fi
ID=\$(cat ${FILELIST} | awk '{print \$NF}' | awk "NR==\${SLURM_ARRAY_TASK_ID}")

if [ -f ${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_paired.fastq ] ; then
	## Trimmed, paired-end
	echo "HISAT2 alignment run on trimmed paired-end reads"
	FP=${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_paired.fastq
	FU=${MAINDIR}/trimmed_fq/\${ID}_trimmed_forward_unpaired.fastq
	RP=${MAINDIR}/trimmed_fq/\${ID}_trimmed_reverse_paired.fastq
	RU=${MAINDIR}/trimmed_fq/\${ID}_trimmed_reverse_unpaired.fastq
	
	hisat2 -p ${CORES} \
	-x $HISATIDX -1 \$FP -2 \$RP -U \${FU},\${RU} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam ${STRANDOPTION} --phred33 \
    ${UNALIGNEDOPT} \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt
	
elif  [ -f ${MAINDIR}/trimmed_fq/\${ID}_trimmed.fastq ] ; then
	## Trimmed, single-end
	echo "HISAT2 alignment run on trimmed single-end reads"
	hisat2 -p ${CORES} \
	-x $HISATIDX -U ${MAINDIR}/trimmed_fq/\${ID}_trimmed.fastq \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam ${STRANDOPTION} --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt

elif [ $PE == "TRUE" ] ; then
	## Untrimmed, pair-end
	echo "HISAT2 alignment run on original untrimmed paired-end reads"
	hisat2 -p ${CORES} \
	-x $HISATIDX -1 \${FILE1} -2 \${FILE2} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam ${STRANDOPTION} --phred33 \
    ${UNALIGNEDOPT} \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt

else
	## Untrimmed, single-end
	echo "HISAT2 alignment run on original untrimmed single-end reads"
	hisat2 -p ${CORES} \
	-x $HISATIDX -U \${FILE1} \
	-S ${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam ${STRANDOPTION} --phred33 \
	2>${MAINDIR}/HISAT2_out/align_summaries/\${ID}_summary.txt
fi


###sam to bam
SAM=${MAINDIR}/HISAT2_out/\${ID}_hisat_out.sam
ORIGINALBAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.bam
SORTEDBAM=${MAINDIR}/HISAT2_out/\${ID}_accepted_hits.sorted

#filter unmapped segments
echo "**** Filtering unmapped segments ****"
date

module load samtools/1.18

samtools view -bh -F 4 ${SAM} > ${ORIGINALBAM}
samtools sort -@ ${CORES} ${ORIGINALBAM} ${SORTEDBAM}
samtools index ${SORTEDBAM}.bam

## Clean up
rm ${SAM}
rm ${ORIGINALBAM}

echo "**** Job ends ****"
date
EOF

call="sbatch ${MAINDIR}/code/speaqeasy_pipeline/.${sname}.sh"
echo $call
$call

