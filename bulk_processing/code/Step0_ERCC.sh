#!/bin/bash
set -euo pipefail

## =========================
## Define variables
## =========================
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

## =========================
## Setup
## =========================
SOFTWARE=/dcs04/hicks/data/jsundstr/bulk_processing
MAINDIR=$/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test
FILELIST=${MAINDIR}/samples.manifest

SHORT="ercc-${EXPERIMENT}"
JOBNAME="step0-${SHORT}.${PREFIX}"

NUM=$(awk '{print $NF}' ${FILELIST} | uniq | wc -l)

mkdir -p logs

## Memory
if [[ "${LARGE}" == "TRUE" ]]; then
    MEM="10G"
else
    MEM="5G"
fi

## Paired-end?
if [[ -f ".paired_end" ]]; then
    PE="TRUE"
else
    PE="FALSE"
fi

## Strandedness
STRANDOPTION=""
case "${STRANDED}" in
    FALSE) ;;
    forward) STRANDOPTION="--fr-stranded" ;;
    reverse) STRANDOPTION="--rf-stranded" ;;
    *)
        echo "Invalid stranded option: ${STRANDED}"
        exit 1 ;;
esac