from os.path import join
from glob import glob
import pandas as pd

configfile: 'config.yml'

#load participants.tsv file (list of HCP subject IDs),
df = pd.read_table(config['participants_tsv'])
subjects = df.participant_id.to_list() 

wildcard_constraints:
    subject="[0-9]+"

rule all:
    input: expand('funcparc/clustering/seed-ZIR_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True)

rule cifti_separate:
    input: lambda wildcards: glob(config['input_dtseries'].format(**wildcards))
    output: 
        lh = 'funcparc/{subject}/input/rfMRI_REST2_7T_AP.L.59k_fs_LR.func.gii',
        rh = 'funcparc/{subject}/input/rfMRI_REST2_7T_AP.R.59k_fs_LR.func.gii'
    singularity: config['singularity_connectomewb']
    shell:
        'wb_command -cifti-separate {input} COLUMN -metric CORTEX_LEFT {output.lh} -metric CORTEX_RIGHT {output.rh}'

rule prepare_subcort:
    input:
        vol = lambda wildcards: glob(config['input_rsvolume'].format(**wildcards)),
        rois = config['subcort_atlas']
    output:
        out = 'funcparc/{subject}/input/rfMRI_REST2_7T_AP_AtlasSubcortical.nii.gz'
    params:
        sigma = 1.6,
        temp = 'funcparc/{subject}/temp'
    #singularity: config['singularity_connectomewb']
    log: 'logs/prepare_subcort/{subject}.log'
    shell: 'scripts/prep_subcortical.sh {input.vol} {input.rois} {params.temp} {params.sigma} {output.out} &> {log}'

rule create_dtseries:
    input: 
        vol = rules.prepare_subcort.output.out,
        rois = config['subcort_atlas'],
        lh = rules.cifti_separate.output.lh,
        rh = rules.cifti_separate.output.rh
    output: 'funcparc/{subject}/input/rfMRI_REST2_7T_AP.59k_fs_LR.dtseries.nii'
    singularity: config['singularity_connectomewb']
    shell:
        'wb_command -cifti-create-dense-timeseries {output} -volume {input.vol} {input.rois} -left-metric {input.lh} -right-metric {input.rh}'

rule extract_confounds:
    input:
        vol = lambda wildcards: glob(config['input_rsvolume'].format(**wildcards)),
        rois = lambda wildcards: glob(config['input_rois'].format(**wildcards)),
        movreg = lambda wildcards: glob(config['input_movreg'].format(**wildcards))
    output: 'funcparc/{subject}/input/confounds.tsv'
    log: 'logs/extract_confounds/{subject}.log'
    script: 'scripts/extract_confounds.py'

rule clean_tseries:
    input:
        dtseries = rules.create_dtseries.output,
        confounds = rules.extract_confounds.output
    output: 'funcparc/{subject}/input_cleaned/rfMRI_REST2_7T_AP.59k_fs_LR.dtseries.nii'
    #singularity: config['singularity_ciftify']
    log: 'logs/clean_dtseries/{subject}.log'
    shell:
        'ciftify_clean_img --output-file={output} --detrend --standardize --confounds-tsv={input.confounds} --low-pass=0.08 --high-pass=0.009 --tr=1 --verbose {input.dtseries} &> {log}'

rule parcellate_tseries:
    input:
        dtseries = rules.clean_tseries.output,
        rois = config['hcpmmp_atlas']
    output: 'funcparc/{subject}/input_parcellated/rfMRI_REST2_7T_AP.59k_fs_LR.ptseries.nii'
    singularity: config['singularity_connectomewb']
    shell:
        'wb_command -cifti-parcellate {input.dtseries} {input.rois} COLUMN {output}'

rule compute_correlation:
    input:
        ptseries = rules.parcellate_tseries.output,
        vol = rules.clean_tseries.output
    output: 'funcparc/{subject}/output/correlation_matrix.npz'
    script: 'scripts/compute_correlation.py'

rule combine_correlation:
    input: expand('funcparc/{subject}/output/correlation_matrix.npz',subject=subjects,allow_missing=True)
    output: 'funcparc/clustering/correlation_matrix_group.npz'
    run:
        import numpy as np

        data = np.load(input[0])
        nsubjects = len(input)
        combined = np.zeros([nsubjects,data['corr'].shape[0],data['corr'].shape[1]])

        for i,npz in enumerate(input):
            data = np.load(npz)
            combined[i,:,:] = data['corr']

        np.savez(output[0], corr_group=combined,indices=data['indices'])

rule spectral_clustering:
    input:
        correlation = rules.combine_correlation.output,
        rois = config['subcort_atlas'],
    params:
        max_k = config['max_k']
    output:
        niftis = expand('funcparc/clustering/seed-ZIR_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True),
        labels = 'funcparc/clustering/cluster_labels.csv'
    script: 'scripts/spectral_clustering.py'
