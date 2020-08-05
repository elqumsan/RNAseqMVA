## SCRIPT NAME:srun_jobarray.mk
##
## AUTHORS: Julien Seiler and Jacques van Helden
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
##
## ## Send the script to the job scheduler
## make -f srun_jobarray.mk [target]

#SBATCH --array=0-13  # Define the IDs for the job array
#SBATCH --mem=32GB # Request 32Gb per job
#SBATCH --cpus=50  # Request 50 CPUs per job
#SBATCH --partition=long  # the long partition (>1 day) is required for some study cases
#SBATCH -o slurm_logs/test_%a_%A_%j_out.txt  # file to store standard output
#SBATCH -e slurm_logs/test-%a_%A_%j_err.txt  # file to store standard error

targets:
	@echo "Targets"
	@echo "	targets		List targets"
	@echo "	list_param	List parameters"
	@echo "	activ_conda	activate conda environment"
	@echo "	one_job		send one job to the queue"
	@echo "	all_jobs	send all jobs to the queue"

## sbatch options
CPUS=50
MEM=32GB
PARTITION=long
#OUT=slurm_logs/test_%a_%A_%j_out.txt
#ERR=slurm_logs/test-%a_%A_%j_err.txt
SBATCH_OPT=--mem=${MEM} --cpus=${CPUS} --partition=${PARTITION}
#-o ${OUT} -e ${ERR}

# Other parameter
RECOUNT_IDS=SRP035988 SRP042620 SRP056295 SRP057196 SRP061240 SRP062966 SRP066834
FEATURE_TYPES=gene transcript
RECOUNT_ID=SRP042620
FEATURE_TYPE=transcript
## Define log directory and files
BASE=/shared/projects/rnaseqmva/RNAseqMVA
WORKSPACE=/shared/projects/rnaseqmva/RNAseqMVA_workspace
LOG_DIR=${WORKSPACE}/logs
START_DATE=`date +%Y-%m-%d_%H%M%S`
PREFIX=${RECOUNT_ID}_${FEATURE_TYPE}_${START_DATE}
list_param:
	@echo "Parameters"
	@echo "	CPUS		${CPUS}"
	@echo "	MEM		${MEM}"
	@echo "	PARTTION	${PARTTION}"
	@echo "	SBATCH_OPT	${SBATCH_OPT}"
	@echo "	RECOUNT_IDS	${RECOUNT_IDS}"
	@echo "	RECOUNT_ID	${RECOUNT_ID}"
	@echo "	FEATURE_TYPES	${FEATURE_TYPES}"
	@echo "	FEATURE_TYPE	${FEATURE_TYPE}"
	@echo "	BASE		${BASE}"
	@echo "	LOG_DIR		${LOG_DIR}"
	@echo "	START_DATE	${START_DATE}"
	@echo "	PREFIX		${PREFIX}"

## The variable $SLURM_ARRAY_TASK_ID takes a different value for each job
## based on the values entered with the argument --array
# srun Rscript --vanilla misc/main_processes.R $SLURM_ARRAY_TASK_ID
##echo "${START_DATE} ${SLURM_ARRAY_TASK_ID}  ${SLURM_ARRAY_JOB_ID} ${RECOUNT_ID} ${FEATURE_TYPE} #${LOG_DIR}/${PREFIX}" >> srun_jobs_sent.tsv

## Initiate the environment
activ_conda:
	module load conda
#	conda init bash
	conda activate rnaseqmva

## Send a job for the analysis of one RECOUNT_ID and FEATURE_TYPE
one_job:
	mkdir -p ${LOG_DIR}
	(cd ${BASE}; srun ${SBATCH_OPT} \
		--output ${LOG_DIR}/${PREFIX}_out.txt \
		--error ${LOG_DIR}/${PREFIX}_err.txt \
		Rscript --vanilla misc/main_processes.R ${RECOUNT_ID} ${FEATURE_TYPE})
