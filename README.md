# Snakemake workflow: Zona Incerta Functional Parcellation

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)

This is a Snakemake workflow for performing functional connectivity-based segmentation of the zona-incerta. It is developed to work with pre-processed HCP 3T and 7T data. Follows the same workflow as the original [Diffusion-based Parcellation] (https://github.com/akhanf/zona-diffparc) scheme.

## Authors

* Ali Khan (@akhanf)
* Roy Haast (@royhaast)

## Usage

### Running on Graham

#### Step 1: Install neuroglia-helpers

`snakemake`, `snakemake_slurm` and `snakemake_remotebatch` are wrappers in [neuroglia-helpers](http://github.com/khanlab/neuroglia-helpers), make sure you use the `graham_khanlab` cfg file.

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

#### Step 3: 

Test your configuration by performing a dry-run via

    snakemake --use-singularity -n

There are a few different ways to execute the workflow:
  1. Execute the workflow locally using an interactive job
  2. Execute the workflow using `snakemake_slurm`
  3. Execute the workflow using `snakemake_remotebatch` (optimized for slurm/graham)

##### Interactive Job

Execute the workflow locally using an interactive job:
    
    regularInteractive -g  

    snakemake --use-singularity --cores 8 --resources gpus=1

##### snakemake_slurm

To execute the workflow for all subjects, submitting a job for each rule group, use:

    snakemake_slurm

##### snakemake_remotebatch

Alternatively, you can use the `snakemake_remotebatch` wrapper to submit in N batches (e.g. 32), and submit all the batches at once, using the `gather_connmap_group` rule to split batches:

    snakemake_remotebatch gather_connmap_group 32

##### Export to Dropbox

To export files to dropbox, use:
    snakemake -s export_dropbox.smk

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

# Step 4: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.

### Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/snakemake-workflows/zona-diffparc) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.


## Testing

No test cases yet

