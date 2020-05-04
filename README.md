# Snakemake workflow: Zona Incerta Functional Parcellation

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)

This is a Snakemake workflow for performing functional connectivity-based segmentation of the [zona-incerta](https://doi.org/10.1101/2020.03.25.008318). It is developed to work with pre-processed HCP 7T data. Follows the same workflow as the original [Diffusion-based Parcellation](https://github.com/akhanf/zona-diffparc) scheme.

## Authors

* Ali Khan (@akhanf)
* Roy Haast (@royhaast)

## Usage
Workflow depends on Snakemake, [Ciftify](https://github.com/edickie/ciftify) and [Connectome Workbench](https://www.humanconnectome.org/software/get-connectome-workbench) to be installed, and the required HCP data (7T fMRI REST + extended) to be unzipped into the ```hcp1200``` folder specified in the config file.

1. On Graham, create virtual environment

```
module load python/3
virtualenv $HOME/venv/funcparc
source $HOME/venv/funcparc/bin/activate
```

2. Install requirements

```
pip install ciftify
pip install snakemake
```

3. Add transparant singularity modules ([Khanlab](https://github.com/khanlab) specific) location to MODULEPATH (in case it is not already)

```
MODULEPATH=/project/6007967/software/transparentsingularity/modules:$MODULEPATH
export MODULEPATH
```

4. After cloning repository, run Snakemake workflow from within folder

```
module load connectome-workbench
snakemake --use-singuliarty
```