# RNAseqMVA

R package developed to perform an evaluation of multivariate analysis
methods for the supervised classification of RNA-seq data.

**Authors: ** 

- Mustafa AbuElQumsan 
- Jacques van Helden (Jacques.van-Helden@univ-amu.fr](mailto:Jacques.van-Helden@univ-amu.fr))

## Downloading

```
git clone
```

## Conda environment

The file [conda-rnaseqmva.yml](conda-rnaseqmva.yml) specifies the parameters of a conda environment enabling to install all required dependencies (R, RCRAN and Bioconductor libraries). 

### Installing the conda rnaseqmva environment

```
conda env create -f conda-rnaseqmva.yml
```

### Loading the environemnt

```
module load conda ## Load the conda module (for the IFB-core-cluster)
source activate rnaseqmva
```

## Configuring the analysis


All the parameters of an analysis can be specified in a YAML file
[misc/00_project_parameters.yml](misc/00_project_parameters.yml). Parameters
can be changed easily by editing this file with any text editor (nano,
gedit, emacs, vi, ...).

## Running the analysis



