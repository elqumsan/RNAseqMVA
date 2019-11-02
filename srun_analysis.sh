#!/bin/bash


## Note: the loading of conda module and environment must apparently be done
## in the sbatch that calls this script, and not within the script itself
# module load conda
# conda init bash
# conda activate rnaseqmva
cd /shared/projects/rnaseqmva/RNAseqMVA
WORKSPACE=/shared/projects/rnaseqmva/RNAseqMVA_workspace
LOG_DIR=${WORKSPACE}/logs
mkdir -p ${LOG_DIR}

cd /shared/projects/rnaseqmva/RNAseqMVA
START_DATE=`date +%Y-%m-%d_%H%M%S`


## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Second analysis

## Parameters
RECOUNT_ID=SRP066834
FEATURE=gene
PREFIX=${RECOUNT_ID}_${FEATURE}_${START_DATE}

## Submit the job to slurm job scheduler via srun
srun --mem=32GB \
  --output ${LOG_DIR}/${PREFIX}_out.txt \
  --error ${LOG_DIR}/${PREFIX}_err.txt \
  Rscript --vanilla misc/main_processes.R ${RECOUNT_ID} ${FEATURE}

