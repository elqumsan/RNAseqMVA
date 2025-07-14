# RNAseqMVA: Benchmarking SVM, Random Forest, and KNN for RNA-seq Data: Revisiting Classifier Performance and the Impact of Preprocessing and Hyperparameter Tuning.

## Developers

- Mustafa AbuElQumsan (ORCID [0000-0002-1018-1410](https://orcid.org/0000-0002-1018-1410)) 
- Jacques van Helden (ORCID [0000-0002-8799-8584](https://orcid.org/0000-0002-8799-8584))


## Description

This repository contains the R code (RNAseqMVA package + analysis scripts) required to reproduce the evaluation described in the following manuscript:

- Mustafa AbuElQumsan, Baddih Gathas and Jacques van Helden (2025). Benchmarking SVM, Random Forest, and KNN for RNA-seq Data: Revisiting Classifier Performance and the Impact of Preprocessing and Hyperparameter Tuning. *Submitted*. 

## Dataset Information

The benchmarking was performed based on 6 datasets downloaded from the Recount database ([rna.recount.bio/](https://rna.recount.bio)). 


| Nb | Name | GEO ID | Pubmed ID | RNA-seq Type | Criterion for class labels | Classes | Samples | Description |
|--|-------------------|-----------|---------|----|---------------|----|----|----------------------------------------------|
| 1 | Breast cancer | SRP042620 | 24929677 | B | tissue | 6 | 167 | Recurrent read-through of fusion transcripts in breast cancer |
| 2 | Acute Myeloid Leukemia | SRP056295 | 26237430 | B | tissue | 4 | 263 | Transcriptome analysis of psoriasis in a large case-control sample |
| 3 | Psoriasis | SRP035988 | 24441097 | B | tissue | 2 | 173 | A survey of human brain transcriptome diversity at the single-cell level |
| 4 | Adult & fetal human brain | SRP057196 | 26060301 | SC | tissue + cell type | 15 | 461 | RNA-seq of systemic lupus erythematosus (SLE) whole blood and healthy controls |
| 5 | Cerebral organoids and fetal neocortex | SRP066834 | 26644564 | SC | tissue | 3 | 729 | Human cerebral organoids recapitulate gene expression programs of fetal neocortex development |
| 6 | Lupus | SRP062966 | 26382853 | SC | disease status | 3 | 117 | A comprehensive study of mutations and gene expression in human acute myeloid leukemia (AML) |


## Code Information

The code is written in R and should run on R version >= 3.6.1. 

It is distributed as an R named `RNAseqMVA`, (for Multi-Variate Analysis of RNA-seq data). 
The core of the package is object-oriented, with classes defined in the [`R`](R/) directory, and scripts in the [`misc`](misc/) directory.
 

## Requirements

The file [conda-rnaseqmva_2025.yml](conda-rnaseqmva_2025.yml) specifies the parameters of a conda environment enabling to install all required dependencies (R, RCRAN and Bioconductor libraries).
If conda is not yet installed on your system, follow the [conda installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/). 

## Usage Instructions

### Cloning this github repository

```bash
git clone https://github.com/elqumsan/RNAseqMVA.git
```

### Installing the conda rnaseqmva environment

The file [conda-rnaseqmva_2025.yml](conda-rnaseqmva_2025.yml) specifies all the requirements to run RNAseqMVA. The simplest way to use the RNAseqMVA package is to create a conda environment that will contain all the dependencies. This can be done automatically with the following commands. 

```bash
cd RNAseqMVA
conda env create -f conda-rnaseqmva_2025.yml
```

__Note__: this command needs to be run only once. The next section explains how to update the environment after it has been created. 

### Updating the conda rnaseqmva environment

In case of changes to the RNAseqMVA environment, an update can be
achieved with the following command. 

```bash
cd RNAseqMVA
conda env update -f conda-rnaseqmva_2025.yml
```

### Activating the conda rnaseqmva environemnt

Check that conda version is >= 5

```
conda --version
```

Activate the conda environment. 

```
conda activate rnaseqmva-2025
```


### Compiling the RNAseqMVA package

We assume here that the RNAseqMVA package has been installed in the
home directory. If not, you just need to adapt the path below.

```bash
cd ~/RNAseqMVA
make roxygenise  ## compile the documentation for the RNAseqMVA package
make build_and_install ## build the RNAseqMVA package
```

### Default configuration

All the parameters of an analysis can be specified in a YAML file [misc/00_project_parameters.yml](misc/00_project_parameters.yml). 
Parameters can be changed easily by editing this file with any text editor (nano, gedit, emacs, vi, ...).

```
  #### Usage : uncomment the recount ID you want to use for the analysis ####

#  selected_recount_ids: ['SRP042620'] # Breast cancer
#  selected_recount_ids: ['SRP056295'] # Acute myeloid leukemia
#  selected_recount_ids: ['SRP035988'] # Psoriasis
#  selected_recount_ids: ['SRP057196'] # Adult and fetal brain cells (sc)
#  selected_recount_ids: ['SRP066834'] # Cerebral organoids and fetal neocortex (sc)
  selected_recount_ids: ['SRP062966'] # Lupus (sc)
```

By default, the configuration is setup to analyse a single study case. The default ID is the one uncommented. If several IDs are uncommented, the analysis will run only on the first one. 
An alternative default ID can be set by uncommenting another row of the proposed `selected_recount_ids`. 

The default study case, feature type (gene or transcript) and number of jobs can be overwritten with command arguments (see details below). 

```
Rscript --vanilla misc/main_processes.R --recountID [ID] --feature [gene|transcript] --jobs [job_nb]
```


### Adding study cases

To add study cases, edit the configuration file and add a section following the examples used in the published article. 
Each study case must be specified with the following fields

- classColumn: the metadata field that contains the relevant class labels
- short_label: used for the figure titles or legends

The valid fields for the class labels are the sub-fields of the column "characteristics". Depending on the study, this column contains either one or several metadata fields. 

For example, in the Psoriasis study case (recountID SRP035988), 

```
characteristics
tissue type: normal skin
tissue type: lesional psoriatic skin
tissue type: lesional psoriatic skin
tissue type: normal skin
tissue type: normal skin
...
```

In some studies, the column "characteristics" contains a combination of several fields. This is the case of the single-cell RNA-seq study on adult and fetal brain cells (recountID SRP057196): 

```
characteristics
c("tissue: cortex", "cell type: fetal_replicating", "age: prenatal 16-18 W", "c1 chip id: nochipID10", "experiment_sample_name: FB_S1")
c("tissue: cortex", "cell type: fetal_replicating", "age: prenatal 16-18 W", "c1 chip id: nochipID10", "experiment_sample_name: FB_S1")
c("tissue: cortex", "cell type: astrocytes", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")
c("tissue: cortex", "cell type: astrocytes", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")
c("tissue: cortex", "cell type: astrocytes", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")
c("tissue: cortex", "cell type: neurons", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")
c("tissue: cortex", "cell type: neurons", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")
c("tissue: cortex", "cell type: microglia", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")
c("tissue: cortex", "cell type: neurons", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")
c("tissue: cortex", "cell type: neurons", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")
c("tissue: cortex", "cell type: astrocytes", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")
c("tissue: cortex", "cell type: neurons", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")
c("tissue: cortex", "cell type: astrocytes", "age: postnatal 47 years", "c1 chip id: nochipID9", "experiment_sample_name: AB_S1")![image](https://github.com/user-attachments/assets/05a925c9-3f64-4256-8ef9-0519a9601f7d)

```

Here is an example of the definition of the study case parameters for these two studies: 

```
SRP035988:  ## Psoriasis (bulk)
  short_label: "Psoriasis"
  classColumn: "tissue.type"

SRP057196: ## Adult and fetal brain cells (sc)
  short_label: "Brain cells (sc)"
  classColumn: ["tissue", "cell.type"] ## For this dataset, it is important to combine the two column in order to capture the goal of the study
```

## Running all analyses

The following command will run the analysis for the stydy case selected above. 

```bash
Rscript --vanilla misc/main_processes.R --recountID [ID] --feature [gene|transcript] --jobs [job_nb]
```

For example, to analyse the breast cancer data (recountID SRP042620) with genes as features and running 11 jobs in parallel, type


```bash
Rscript --vanilla misc/main_processes.R --recountID SRP042620 --feature gene --jobs 11
```



The script [`script misc/main_process.R`](misc/main_process.R), calls a series of other scripts to run the successive steps of the analysis in the right order. 


### Running selected analyses

The analyses can also be led step-by-step by opening the project in RStudio (via the project configuration file [`RNAseqMVA.Rproj`](RNAseqMVA.Rproj)). 

Once there, you first need to 

1. Ensure that the dependencies are present in your RStudio environment. 
2. Compile the RNAseqMVA package,
3. Open the file [`script misc/main_process.R`](misc/main_process.R)
4. Run the lines one by one. 

### Specific settings for the core cluster of the Institut Français de Bioinformatique (IFB-core-cluster)

This section is specific to the core cluster of the Institut Français de Bioinformatique (IFB-core-cluster), which was used to run comparative assessment of supervised classification methods for RNA-seq.

If you are working on another infrastructure, you can skip it. 

On the [IFB core cluster](https://www.france-bioinformatique.fr/en/ifb-core-cluster/), conda is loaded via a module, which must be loaded with the following command: 

```bash
## Load required modules for the IFB-core cluster
module load conda ## Load the conda module 
module load r     ## Load the R module
```

After that, the RNAseqMVA environment can be built and activated and the analyses launched in the same way as described in the previous sections, but **each task has to be send to the job scheduler slurm**, via either `sbatch` (sending a script or `srun` (submitting a single-line command).

```bash
srun make roxygenise          ## Generate the user documentation for RNAseqMVA package
srun make build_and_install   ## Build the RNAseqMVA package

## Send the analysis for one study case to slurm job scheduler
##
## Note: the parameters below have to be adapted according to the particular cluster configuration
## - we use the long partition (queue) because for some study cases the analysis can last for more than 24h
## - in the file misc/00_project_parameters.yml, we set jobs: 25
## - we allocate 64Gb of memory in total, which makes >2Gb per CPU, in principle sufficient
srun --partition=long --mem=64GB Rscript --vanilla misc/main_processes.R
```

## Methodology – Steps taken for data processing or modeling.

Data preprocessing is managed by the `RNAseqMVA` class `StudyCase` according to the following steps.

- **Data download**. RNA-seq data (raw counts) and metadata (pheno table) is downloaded from the Recount dabase (the user-selected study case is specified by its Recount ID)
- **Class filtering**. The classes having less than 15 samples are discarded (the min number of samples per class can by modified by editing the configuration file [`misc/00_project_parameters.yml`](misc/00_project_parameters.yml))
- **Feature filtering**. Features (genes or transcript) having a zero variance are suppressed from the data set. 
- **Library size standardization**. Four alternative methods are applied to standardize the library sizes: 

    - upper quartile (q0.75), 
    - weighted trimmed mean of M-values (TMM) (Robinson & Oshlack, 2010), 
    - relative log expression (RLE) from DESeq (Anders & Huber, 2010), 
    - the `estimateSizeFactors() function implemented in the DESeq2 package (Love, Huber & Anders, 2014).  

- **log2 transformation**. Standardized counts are normalized by log2 transformation, after having added an epsilon to avoid problems with zero counts. 
- **Principal component computation**. The principal components are computed from the log2-normalised, library-size standardized, counts.

These steps are described in detail in the manuscript. 

## Limitations

This comparative benchmarking is restricted to 3 supervised classification approaches : Random Forest (RF), Support Vector Machines (SVM) and K-nearest neighbours (KNN). These methods were chosen because they rely on different models and algorithmic approaches to the problem of supervised classification. A systematic evaluation of all existing methods was out of scope for this work, but the R code could be adapted to handle other methods as well. 

## Citations

- Mustafa AbuElQumsan, Baddih Gathas and Jacques van Helden (2025). Benchmarking SVM, Random Forest, and KNN for RNA-seq Data: Revisiting Classifier Performance and the Impact of Preprocessing and Hyperparameter Tuning. *Submitted*.

## License & Contribution Guidelines

[![License](https://img.shields.io/github/license/elqumsan/RNAseqMVA)](LICENSE)

## Materials & Methods

All the datasets used for this evaluation were downloaded from the Recount database ([rna.recount.bio/](https://rna.recount.bio)). 

This code evaluates three alternative methods for the supervised classification of RNA-seq data: 

- Random Forest (RF)
- K nearest neighbours (KNN)
- Support Vector Machines (SVM)

**Evaluation method.** Each classifier is evaluated by measuring the Misclassification Error Rate, with iterative cross-validation (subsampling of the original data set).

**Subsampling.** The evaluation was based on 50 iterations of random subsampling with 2/3 individuals for training and 1/3 for testing. Subsampling was stratified in order to ensure a fair representation of all the classes in both training and testing sets. 

A full description of the methodology is provided in the submitted manuscript mentioned above. 


## Conclusions

- Limitations: Identify limitations in your study


