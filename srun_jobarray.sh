#!/bin/bash

## Script: srun_jobarray.sh
##
## Authors: Julien Seiler and Jacques van Helden
##
## DESCRIPTION
## Run the analysis of the different study cases (recount IDs) and feature types
## on IFB core cluster via a job array handled by slurm job scheduler.
## All jobs are sent in one shot, and are then handled in parallel according to
## the availability of resources on the cluster).
##
## USAGE
## Change directory to RNAseqMVA package
## cd /shared/projects/rnaseqmva/RNAseqMVA
##
## ## Initiate the environment
## ## This is necessary for the sbatch command !
## module load conda
## conda init bash
## conda activate rnaseqmva
##
## ## Send the script to the job scheduler
## sbatch srun_jobarray.sh

#SBATCH --array=0-13  # Define the IDs for the job array
#SBATCH --mem=32GB # Request 32Gb per job
#SBATCH --cpus=50  # Request 50 CPUs per job
#SBATCH --partition=long  # the long partition (>1 day) is required for some study cases
#SBATCH -o slurm_logs/test_%a_%A_%j_out.txt  # file to store standard output
#SBATCH -e slurm_logs/test-%a_%A_%j_err.txt  # file to store standard error

## Initiate the environment
module load conda
conda init bash
conda activate rnaseqmva

mkdir -p slurm_logs

## Define the parameters
RECOUNT_IDS=(SRP035988 SRP042620 SRP056295 SRP057196 SRP061240 SRP062966 SRP066834)
FEATURE_TYPES=(gene transcript)

## Associate RECOUNT_ID and FEATURE_TYPE to the index of the current job
recount_index=$((SLURM_ARRAY_TASK_ID / 2))
RECOUNT_ID=${RECOUNT_IDS[recount_index]}
feature_index=$((SLURM_ARRAY_TASK_ID % 2))
FEATURE_TYPE=${FEATURE_TYPES[feature_index]}

## Define log directory and files
cd /shared/projects/rnaseqmva/RNAseqMVA
WORKSPACE=/shared/projects/rnaseqmva/RNAseqMVA_workspace
LOG_DIR=${WORKSPACE}/logs
mkdir -p ${LOG_DIR}
START_DATE=`date +%Y-%m-%d_%H%M%S`
PREFIX=${RECOUNT_ID}_${FEATURE_TYPE}_${START_DATE}

## The variable $SLURM_ARRAY_TASK_ID takes a different value for each job
## based on the values entered with the argument --array
# srun Rscript --vanilla misc/main_processes.R $SLURM_ARRAY_TASK_ID
echo "${START_DATE} ${SLURM_ARRAY_TASK_ID}  ${SLURM_ARRAY_JOB_ID} ${RECOUNT_ID} ${FEATURE_TYPE} ${LOG_DIR}/${PREFIX}" >> srun_jobs_sent.tsv

## Send a job for the analysis of one RECOUNT_ID and FEATURE_TYPE
srun --mem=32GB --cpus=50 \
  --output ${LOG_DIR}/${PREFIX}_out.txt \
  --error ${LOG_DIR}/${PREFIX}_err.txt \
  Rscript --vanilla misc/main_processes.R ${RECOUNT_ID} ${FEATURE_TYPE}
