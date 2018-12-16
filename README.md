# RNAseqMVA

R package developed to perform an evaluation of multivariate analysis
methods for the supervised classification of RNA-seq data.

**Authors: ** 

- Mustafa AbuElQumsan 
- Jacques van Helden (Jacques.van-Helden@univ-amu.fr](mailto:Jacques.van-Helden@univ-amu.fr))

## Downloading

```
git clone https://github.com/elqumsan/RNAseqMVA.git
```

## Conda environment

The file [conda-rnaseqmva.yml](conda-rnaseqmva.yml) specifies the parameters of a conda environment enabling to install all required dependencies (R, RCRAN and Bioconductor libraries). 

### Installing the conda rnaseqmva environment

```
cd RNAseqMVA
conda env create -f conda-rnaseqmva.yml
```

### Loading conda module

On the IFB cluster, conda is loaded via a module

```
module load conda ## Load the conda module (for the IFB-core-cluster)
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
cd ~/RNAseqMVA
Rscript --vanilla misc/main_processes.R
```

## Running selected analyses

```
cd ~/RNAseqMVA
R --vanilla
```

Then open the file [misc/main_processes.R](misc/main_processes.R) and
decide whether you need to run each command.
