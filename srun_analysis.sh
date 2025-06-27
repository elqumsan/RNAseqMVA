#!/bin/bash

## ================================================================
## Script: srun_analysis.sh
## Authors: Julien Seiler and Jacques van Helden
##
## USAGE
##
## ## Initiate the environment
## ## This is necessary before running the sbatch command !
## module load conda
## conda init bash
## conda activate rnaseqmva-2025
##
## ## Send the script to the job scheduler
## sbatch --mem=48GB --cpus=25  --partition=long srun_analysis.sh
## 
## Note: the loading of conda module and environment must apparently be done
## before ruuning sbatch to call this script, and not within the script itself

## Define log directory and files
cd /shared/projects/rnaseqmva/RNAseqMVA
WORKSPACE=/shared/projects/rnaseqmva/RNAseqMVA_workspace
LOG_DIR=${WORKSPACE}/logs
mkdir -p ${LOG_DIR}
START_DATE=`date +%Y-%m-%d_%H%M%S`


## ================================================================
## Run a given analysis
##
## Note: to run all jobs, it is recommended to use job arrays
## (see script srun_jobarray.sh)


#### Parameters ####

## Recount identifier: uncomment one of the following
# RECOUNT_ID=SRP035988
# RECOUNT_ID=SRP042620
RECOUNT_ID=SRP056295
# RECOUNT_ID=SRP057196
# RECOUNT_ID=SRP061240
# RECOUNT_ID=SRP062966
# RECOUNT_ID=SRP066834

## Feature type: uncomment one of the following
# FEATURE=transcript
FEATURE=gene

## Prefix for the log and err files
PREFIX=${RECOUNT_ID}_${FEATURE}_${START_DATE}

## slurm parameters
CPUS=25
MEM=64GB
PARTITION=long

## Report parameters
echo "RECOUNT_ID: ${RECOUNT_ID}"
echo "FEATURE: ${FEATURE}"
echo "PREFIX: ${PREFIX}"
echo "CPUS: ${CPUS}"
echo "MEM: ${MEM}"
echo "PARTITION: ${PARTITION}"

## Submit the job to slurm job scheduler via srun
srun --mem=${MEM} --cpus=${CPUS} --partition=${PARTITION} \
  --output ${LOG_DIR}/${PREFIX}_out.txt \
  --error ${LOG_DIR}/${PREFIX}_err.txt \
  Rscript --vanilla misc/main_processes.R --recountID ${RECOUNT_ID} --feature ${FEATURE}
echo "PREFIX  ${PREFIX}"
echo "  output  ${LOG_DIR}/${PREFIX}_out.txt"
echo "  error   ${LOG_DIR}/${PREFIX}_err.txt"


## Command to send the script to the job scheduler
## sbatch --mem=${MEM} --cpus=${CPUS}  --partition=${PARTITION} srun_analysis.sh
## squeue -u ${UID}

