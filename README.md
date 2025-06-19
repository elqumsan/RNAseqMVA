# RNAseqMVA

## Developers

- Mustafa AbuElQumsan (ORCID [https://orcid.org/0000-0002-1018-1410](0000-0002-1018-1410)) 
- Jacques van Helden (ORCID [https://orcid.org/0000-0002-8799-8584](0000-0002-8799-8584))

## Name of the project

Benchmarking SVM, Random Forest, and KNN for RNA-seq Data: Revisiting Classifier Performance and the Impact of Preprocessing and Hyperparameter Tuning.

## Description

This repository contains the R code (RNAseqMVA package + analysis scripts) required to reproduce the evaluation descrived in the following manuscript:

- Mustafa AbuElQumsan, Baddih Gathas and Jacques van Helden (2025). Benchmarking SVM, Random Forest, and KNN for RNA-seq Data: Revisiting Classifier Performance and the Impact of Preprocessing and Hyperparameter Tuning. _Submitted_. 

## Dataset Information

The benchmarking was performed based on 6 datasets downloaded from the Recount database ([rna.recount.bio/](https://rna.recount.bio)). 

## Code Information

The code is written in R and should run on R version >= 3.6.1. 

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
make build_and_install
```

### Configuring the analysis

All the parameters of an analysis can be specified in a YAML file [misc/00_project_parameters.yml](misc/00_project_parameters.yml). 

Parameters can be changed easily by editing this file with any text editor (nano, gedit, emacs, vi, ...).

## Running all analyses

```bash
Rscript --vanilla misc/main_processes.R
```
This command will run the `script misc/main_process.R`, which will call other scripts in the right order to lead the successive steps of the analysis. 


### Running selected analyses

```bash
R --vanilla
```

Then open the file [misc/main_processes.R](misc/main_processes.R) and identify the scripts you need to run separately.

### Specific settings for the IFB cluster 

This section is specific to the core cluster of the Institut Français de Bioinformatique (IFB-core-cluster), which was used to run comparative assessment of supervised classification methods for RNA-seq.

If you are working on another infrastructure, you can skip it. 

On the [IFB core cluster](https://www.france-bioinformatique.fr/cluster), conda is loaded via a module, which must be loaded with the following command: 

```bash
module load conda ## Load the conda module (for the IFB-core-cluster)
```

after that, the RNAseqMVA environment can be loaded as described above. 

Commands are sent to cluster nodes via srun. 

```bash
srun --mem=32GB Rscript --vanilla misc/main_processes.R
```


## Methodology (if applicable) – Steps taken for data processing or modeling.

## Citations (if applicable) – If this dataset was used in research, provide references.

## License & Contribution Guidelines (if applicable).

## Materials & Methods

Include 3rd party dataset DOI/URL in the main text: Any dataset you have used that has been curated and uploaded by an external source.
Evaluation method: The evaluation method used to evaluate the proposed technique. Evaluation methods (e.g., ablation study, cross-validation, cross-dataset testing) refer to the APPROACH or PROCEDURE used to validate the model’s effectiveness.
Conclusions


