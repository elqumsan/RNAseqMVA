################################################################
## Define some recurrent tasks

targets:
	@echo "Targets"
	@echo "	parameters		display the general parameters"
	@echo " build_and_install	Build and install RNAseqMVA package for R"
	@echo "	convert_pc_plots_all_study_cases	convert PC plots from pdf to ${IMG_FORMAT} for all study cases"
	@echo "	sync_pc_plots_all_study_cases	synchronize PC plots from workspace to the manuscript folder"
	@echo "	ws_dir_to_rsatix	synchronize a directory (default: TO_SYNC=${TO_SYNC}) from your workspace to the shared space on rsatix"
	@echo "	ws_dir_from_rsatix	synchronize a directory from the workspace on rsatix to your workspace"

# SRP035988 SRP042620 SRP056295 SRP057196 SRP061240 SRP062966 SRP066834
IMG_FORMAT=png
WORKSPACE=~/RNAseqMVA_workspace
STUDY_CASES=`ls -1 ~/RNAseqMVA_workspace/results/  | grep SRP | xargs`
parameters:
	@echo "Parameters"
	@echo "	IMG_FORMAT	${IMG_FORMAT}"
	@echo "	STUDY_CASES	${STUDY_CASES}"
	@echo "	STUDY_CASE	${STUDY_CASE}"
	@echo "	PCS		${PCS}"
	@echo "	MANUSCRIPT_DIR	${MANUSCRIPT_DIR}"
	@echo "	FIGURE_ELEMENTS	${FIGURE_ELEMENTS}"

################################################################
## Build and install the RNAseqMVA package for R
build_and_install:
	@echo ""
	@echo "Building RNAseqMVA package"
	(cd ..; R CMD build RNAseqMVA)
	@echo "Installing RNAseqMVA package in R"
	(cd ..; R CMD INSTALL RNAseqMVA_0.1.2.tar.gz)

################################################################
## Convert PC plots from pdf (generated by R) to some other image format
STUDY_CASE=SRP042620
PCS=PC1to4
SOURCE_PLOT=${STUDY_CASE}_log2norm_${PCS}.pdf
SOURCE_PC=${WORKSPACE}/results/${STUDY_CASE}/visualization_of_PCs/${SOURCE_PLOT}
MANUSCRIPT_DIR=~/Dropbox/Mustafa_AbuElQumsan_shared_folder/scientific/research_project_Mustafa/manuscript_supervised_classif_assessment_RNA-seq
FIGURE_ELEMENTS=${MANUSCRIPT_DIR}/figure_elements
TARGET_PC=${FIGURE_ELEMENTS}/${STUDY_CASE}_log2norm_${PCS}.${IMG_FORMAT}
convert_one_pcplot:
	@mkdir -p ${FIGURE_ELEMENTS}
	@convert -flatten -density 300 \
		${SOURCE_PC} ${TARGET_PC}
	@echo "	${TARGET_PC}"

convert_pc_plots_one_study_case:
	@echo "Converting PC plot for study case	${STUDY_CASE}"
	@${MAKE} convert_one_pcplot PCS=PC1-PC2
	@${MAKE} convert_one_pcplot PCS=PC3-PC4
	@${MAKE} convert_one_pcplot PCS=PCplots


convert_pc_plots_all_study_cases:
	@echo "Converting PC plots"
	@for i in ${STUDY_CASES}; do \
		${MAKE} convert_pc_plots_one_study_case STUDY_CASE=$${i}; \
	done


EXCLUDE_OPT=--exclude *~ --exclude .DS_Store --exclude .\#*
sync_one_pcplot:
	@mkdir -p ${FIGURE_ELEMENTS}
	rsync ${EXCLUDE_OPT}${SOURCE_PC} ${FIGURE_ELEMENTS}/


sync_pc_plots_one_study_case:
	@echo "Syncing PC plot for study case	${STUDY_CASE}"
	@${MAKE} sync_one_pcplot PCS=PC1-PC2
	@${MAKE} sync_one_pcplot PCS=PC3-PC4
	@${MAKE} sync_one_pcplot PCS=PCplots


sync_pc_plots_all_study_cases:
	@echo "Syncing PC plots"
	@for i in ${STUDY_CASES}; do \
		${MAKE} sync_pc_plots_one_study_case STUDY_CASE=$${i}; \
	done

RSATIX=rsat-tagc.univ-mrs.fr
LOGIN=rnaseqmva
RSATIX_LOGIN=rnaseqmva@${RSATIX}
RSATIX_WS=/workspace/RNAseqMVA/RNAseqMVA_workspace
LOCAL_WS=~/RNAseqMVA_workspace/
TO_SYNC=memory_images
RSYNC_OPT=
OUT_SOURCE=${LOCAL_WS}/${TO_SYNC}
OUT_TARGET=${LOGIN}@${RSATIX}:${RSATIX_WS}/
ws_dir_to_rsatix:
	@echo "Synchronizing directory from your workspace to shared space on ${RSATIX}"
	@echo "	LOCAL_WS	${LOCAL_WS}"
	@echo "	TO_SYNC		${TO_SYNC}"
	@echo "	RSYNC_OPT	${RSYNC_OPT}"
	@echo "	OUT_SOURCE	${OUT_SOURCE}"
	@echo "	OUT_TARGET	${OUT_TARGET}"
	@echo "	EXCLUDE_OPT	${EXCLUDE_OPT}"
	@rsync -ruptvl ${EXCLUDE_OPT} ${RSYNC_OPT} ${OUT_SOURCE} ${OUT_TARGET}

IN_SOURCE=${LOGIN}@${RSATIX}:${RSATIX_WS}/${TO_SYNC}
IN_TARGET=${LOCAL_WS}
ws_dir_from_rsatix:
	@echo "Synchronizing directory from ${RSATIX} shared space to your workspace"
	@echo "	LOCAL_WS	${LOCAL_WS}"
	@echo "	TO_SYNC		${TO_SYNC}"
	@echo "	RSYNC_OPT	${RSYNC_OPT}"
	@echo "	IN_SOURCE	${IN_SOURCE}"
	@echo "	IN_TARGET	${IN_TARGET}"
	@echo "	EXCLUDE_OPT	${EXCLUDE_OPT}"
	@rsync -ruptvl ${EXCLUDE_OPT} ${RSYNC_OPT} ${IN_SOURCE} ${IN_TARGET}

