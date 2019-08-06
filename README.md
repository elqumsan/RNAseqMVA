# RNAseqMVA

R package developed to perform an evaluation of multivariate analysis
methods for the supervised classification of RNA-seq data.

## Authors

- Mustafa AbuElQumsan 
- Jacques van Helden (Jacques.van-Helden@univ-amu.fr](mailto:Jacques.van-Helden@univ-amu.fr))

## Downloading

```
git clone https://github.com/elqumsan/RNAseqMVA.git
```

## Conda environment

The file [conda-rnaseqmva.yml](conda-rnaseqmva.yml) specifies the parameters of a conda environment enabling to install all required dependencies (R, RCRAN and Bioconductor libraries).

### Loading conda module 

This step is specific to the IFB core cluster. If you are working on another infrastructure, you can skip it. 

On the [IFB core cluster](https://www.france-bioinformatique.fr/cluster), conda is loaded via a module, so this command must be run before starting conda

```
module load conda ## Load the conda module (for the IFB-core-cluster)
```

### Installing the conda rnaseqmva environment

The file [conda-rnaseqmva.yml](conda-rnaseqmva.yml) specifies all the requirements to run RNAseqMVA. The simplest way to use the RNAseqMVA package is to create a conda environment that will contain all the dependencies. This can be done automatically with the following commands. 

```
cd RNAseqMVA
conda env create -f conda-rnaseqmva.yml
```

### Updating the conda rnaseqmva environment

In case of changes to the RNAseqMVA environment, an update can be
achieved with the following command. 

```
cd RNAseqMVA
conda env update -f conda-rnaseqmva.yml
```

### Loading the environemnt

```
source activate rnaseqmva
```

## Compiling the RNAseqMVA package

We assume here that the RNAseqMVA package has been installed in the
home directory. If not, you just need to adapt the path below.

```
make build_and_install
```

## Configuring the analysis

All the parameters of an analysis can be specified in a YAML file
[misc/00_project_parameters.yml](misc/00_project_parameters.yml). Parameters
can be changed easily by editing this file with any text editor (nano,
gedit, emacs, vi, ...).

## Running all analyses

```
Rscript --vanilla misc/main_processes.R
```

### On the IFB cluster

On the IFB cluster, commands are sent to cluster nodes via srun. 

```
srun --mem=32GB Rscript --vanilla misc/main_processes.R
```


## Running selected analyses

```
R --vanilla
```

Then open the file [misc/main_processes.R](misc/main_processes.R) and
decide whether you need to run each command.
