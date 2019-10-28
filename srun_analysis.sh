#!/bin/bash

module load conda
conda init bash
conda activate rnaseqmva
cd /shared/projects/rnaseqmva/RNAseqMVA
LOG_DIR=/shared/projects/rnaseqmva/logs
mkdir -p ${LOG_DIR}

cd /shared/projects/rnaseqmva/RNAseqMVA

START_DATE=`date +%Y-%m-%d_%H%M%S`
OUT_FILE=${LOG_DIR}/log_${START_DATE}_out.txt
ERR_FILE=${LOG_DIR}/log_${START_DATE}_err.txt

## Submit the job to slurm job scheduler via srun
srun --mem=32GB \
  --output ${OUT_FILE} \
  --error ${ERR_FILE} \
  Rscript --vanilla misc/main_processes.R
