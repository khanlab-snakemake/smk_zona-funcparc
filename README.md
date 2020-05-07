# Snakemake workflow: Zona Incerta Functional Parcellation

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)

This is a Snakemake workflow for performing functional connectivity-based segmentation of the [zona-incerta](https://doi.org/10.1101/2020.03.25.008318). It is developed to work with pre-processed HCP 7T data. Follows the same workflow as the original [Diffusion-based Parcellation](https://github.com/akhanf/zona-diffparc) scheme.

## Authors

* Ali Khan (@akhanf)
* Roy Haast (@royhaast)

## Usage
Workflow depends on Snakemake, and [Ciftify](https://github.com/edickie/ciftify) and [Connectome Workbench](https://www.humanconnectome.org/software/get-connectome-workbench) singularity images to be present. Optimized for HCP data (requires: 7T_REST_1.6mm_preproc + extended, and 3T_Structural_1.6mm_preproc packages).

After cloning repository, run Snakemake workflow from within folder:

```
snakemake --use-singularity --singularity-args '\-e'
```

Or in case of large number of subjects/jobs, you can setup the [cc-slurm profile](https://github.com/khanlab/cc-slurm) for parallel processing of subjects:

```
snakemake --profile cc-slurm
```

## TO DO
- Add rules to evaluate number of clusters
- Increase flexibility to allow other seed ROIs