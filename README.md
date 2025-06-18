# RNAseqMVA

This code contains 
- the R package RNAseqMVA, implemented to perform an evaluation of multivariate analysis
methods for the supervised classification of RNA-seq data.
- the R scripts required to reproduce the evaluation descrived in the manuscript

Mustafa AbuElQumsan, Baddih Gathas and Jacques van Helden (2025). 

## Authors

- Mustafa AbuElQumsan 
- Jacques van Helden ([Jacques.van-Helden@univ-amu.fr](mailto:Jacques.van-Helden@univ-amu.fr))

## Downloading

```bash
git clone https://github.com/elqumsan/RNAseqMVA.git
```

## Conda environment

The file [conda-rnaseqmva.yml](conda-rnaseqmva.yml) specifies the parameters of a conda environment enabling to install all required dependencies (R, RCRAN and Bioconductor libraries).

## Installing miniconda

If conda is not yet installed on your system, follow the [conda installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/). 


### Loading conda module 

This step is specific to the IFB core cluster. If you are working on another infrastructure, you can skip it. 

On the [IFB core cluster](https://www.france-bioinformatique.fr/cluster), conda is loaded via a module, so this command must be run before starting conda

```
module load conda ## Load the conda module (for the IFB-core-cluster)
```

### Installing the conda rnaseqmva environment

The file [conda-rnaseqmva.yml](conda-rnaseqmva.yml) specifies all the requirements to run RNAseqMVA. The simplest way to use the RNAseqMVA package is to create a conda environment that will contain all the dependencies. This can be done automatically with the following commands. 

```bash
cd RNAseqMVA
conda env create -f conda-rnaseqmva.yml
```

__Note__: this command needs to be run only once. The next section explains how to update the environment after it has been created. 

### Updating the conda rnaseqmva environment

In case of changes to the RNAseqMVA environment, an update can be
achieved with the following command. 

```bash
cd RNAseqMVA
conda env update -f conda-rnaseqmva.yml
```

### Loading the environemnt

The next command needs to be adapted depending on your conda version. 

```
conda --version
```

If your version is < 5, use the command `source` below.

```bash
source activate rnaseqmva
```

If you have conda verison >=5, you can run `conda activate`as below.

```
conda activate rnaseqmva
```


## Compiling the RNAseqMVA package

We assume here that the RNAseqMVA package has been installed in the
home directory. If not, you just need to adapt the path below.

```bash
cd ~/RNAseqMVA
make build_and_install
```

## Configuring the analysis

All the parameters of an analysis can be specified in a YAML file
[misc/00_project_parameters.yml](misc/00_project_parameters.yml). Parameters
can be changed easily by editing this file with any text editor (nano,
gedit, emacs, vi, ...).

## Running all analyses

```bash
Rscript --vanilla misc/main_processes.R
```


## Running selected analyses

```bash
R --vanilla
```

Then open the file [misc/main_processes.R](misc/main_processes.R) and
decide whether you need to run each command.

## Specific settings for the IFB cluster 

The cluster of the Institut Français de Bioinformatique (IFB-core-cluster) was used to run comparative assessment of supervised classification methods for RNA-seq.

On this cluster, conda is loaded via a module. 

```bash
module load conda ## Load the conda module (for the IFB-core-cluster)
```

after that, the RNAseqMVA environment can be loaded as described above. 

Commands are sent to cluster nodes via srun. 

```bash
srun --mem=32GB Rscript --vanilla misc/main_processes.R
```

