from os.path import join
from glob import glob
import pandas as pd

configfile: 'config.yml'


#load participants.tsv file, and strip off sub- from participant_id column
df = pd.read_table(config['participants_tsv'])
subjects = df.participant_id.to_list() 
subjects = [ s.strip('sub-') for s in subjects ]


#get list of ROIs
f = open(config['targets_txt'],'r')
targets = [line.strip() for line in f.readlines()]
f.close()

#get seeds and hemis from config
seeds = config['seeds']
hemis = config['hemis']



wildcard_constraints:
    subject="[a-zA-Z0-9]+"



rule all:
    input: 
        clusters = expand('diffparc/clustering/group_space-{template}_seed-{seed}_hemi-{hemi}_method-spectralcosine_k-{k}_cluslabels.nii.gz',seed=seeds,hemi=hemis,template=config['template'],k=range(2,config['max_k']+1))



rule import_targets:
    input: 
        lh = join(config['targets_seg_dir'],config['targets_seg_lh']),
        rh = join(config['targets_seg_dir'],config['targets_seg_rh'])
    output: 'diffparc/sub-{subject}/masks/lh_rh_targets_native.nii.gz'
    envmodules: 'fsl'
    singularity: config['singularity_neuroglia']
    log: 'logs/import_targets_hcp_mmp_sym/sub-{subject}.log'
    group: 'pre_track'
    shell:
        'fslmaths {input.lh} -max {input.rh} {output} &> {log}'

rule import_template_seed:
    input: join(config['template_seg_dir'],config['template_seg_nii'])
    output: 'diffparc/template_masks/sub-{template}_hemi-{hemi}_desc-{seed}_mask.nii.gz'
    log: 'logs/import_template_seed/{template}_{seed}_{hemi}.log'
    group: 'pre_track'
    shell: 'cp -v {input} {output} &> {log}'


rule transform_to_subject:
    input: 
        seed = rules.import_template_seed.output,
        affine = lambda wildcards: glob(join(config['ants_template_dir'],config['ants_affine_mat'].format(**wildcards))),
        invwarp = lambda wildcards: glob(join(config['ants_template_dir'],config['ants_invwarp_nii'].format(**wildcards))),
        ref = 'diffparc/sub-{subject}/masks/lh_rh_targets_native.nii.gz'
    output: 'diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_{hemi}.nii.gz'
    envmodules: 'ants'
    singularity: config['singularity_neuroglia']
    log: 'logs/transform_to_subject/{template}_sub-{subject}_{seed}_{hemi}.log'
    group: 'pre_track'
    shell:
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.seed} -o {output} -r {input.ref} -t [{input.affine},1] -t {input.invwarp} &> {log}'
    
    
rule resample_targets:
    input: 
        dwi = join(config['prepdwi_dir'],'bedpost','sub-{subject}','mean_S0samples.nii.gz'),
        targets = 'diffparc/sub-{subject}/masks/lh_rh_targets_native.nii.gz'
    params:
        seed_resolution = config['probtrack']['seed_resolution']
    output:
        mask = 'diffparc/sub-{subject}/masks/brain_mask_dwi.nii.gz',
        mask_res = 'diffparc/sub-{subject}/masks/brain_mask_dwi_resampled.nii.gz',
        targets_res = 'diffparc/sub-{subject}/masks/lh_rh_targets_dwi.nii.gz'
    singularity: config['singularity_neuroglia']
    log: 'logs/resample_targets/sub-{subject}.log'
    group: 'pre_track'
    shell:
        'fslmaths {input.dwi} -bin {output.mask} &&'
        'mri_convert {output.mask} -vs {params.seed_resolution} {params.seed_resolution} {params.seed_resolution} {output.mask_res} -rt nearest &&'
        'reg_resample -flo {input.targets} -res {output.targets_res} -ref {output.mask_res} -NN 0  &> {log}'

rule resample_seed:
    input: 
        seed = rules.transform_to_subject.output,
        mask_res = 'diffparc/sub-{subject}/masks/brain_mask_dwi_resampled.nii.gz'
    output:
        seed_res = 'diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_{hemi}_resampled.nii.gz',
    singularity: config['singularity_neuroglia']
    log: 'logs/resample_seed/{template}_sub-{subject}_{seed}_{hemi}.log'
    group: 'pre_track'
    shell:
        'reg_resample -flo {input.seed} -res {output.seed_res} -ref {input.mask_res} -NN 0 &> {log}'

    
    

rule split_targets:
    input: 
        targets = 'diffparc/sub-{subject}/masks/lh_rh_targets_dwi.nii.gz',
    params:
        target_nums = lambda wildcards: [str(i) for i in range(len(targets))]
    output:
        target_seg = expand('diffparc/sub-{subject}/targets/{target}.nii.gz',target=targets,allow_missing=True)
    singularity: config['singularity_neuroglia']
    log: 'logs/split_targets/sub-{subject}.log'
    threads: 32 
    group: 'pre_track'
    shell:
        'parallel  --jobs {threads} fslmaths {input.targets} -thr {{1}} -uthr {{1}} -bin {{2}} &> {log} ::: {params.target_nums} :::+ {output.target_seg}'

rule gen_targets_txt:
    input:
        target_seg = expand('diffparc/sub-{subject}/targets/{target}.nii.gz',target=targets,allow_missing=True)
    output:
        target_txt = 'diffparc/sub-{subject}/target_images.txt'
    log: 'logs/get_targets_txt/sub-{subject}.log'
    group: 'pre_track'
    run:
        f = open(output.target_txt,'w')
        for s in input.target_seg:
            f.write(f'{s}\n')
        f.close()


rule run_probtrack:
    input:
        seed_res = rules.resample_seed.output,
        target_txt = rules.gen_targets_txt.output,
        mask = 'diffparc/sub-{subject}/masks/brain_mask_dwi.nii.gz',
        target_seg = expand('diffparc/sub-{subject}/targets/{target}.nii.gz',target=targets,allow_missing=True)
    params:
        bedpost_merged = join(config['prepdwi_dir'],'bedpost','sub-{subject}','merged'),
        probtrack_opts = config['probtrack']['opts'],
        probtrack_dir = 'diffparc/sub-{subject}/probtrack_{template}_{seed}_{hemi}'
    output:
        target_seg = expand('diffparc/sub-{subject}/probtrack_{template}_{seed}_{hemi}/seeds_to_{target}.nii.gz',target=targets,allow_missing=True)
    threads: 2
    resources: 
        mem_mb = 8000, 
        time = 30, #30 mins
        gpus = 1 #1 gpu
    log: 'logs/run_probtrack/{template}_sub-{subject}_{seed}_{hemi}.log'
    shell:
        'probtrackx2_gpu --samples={params.bedpost_merged}  --mask={input.mask} --seed={input.seed_res} ' 
        '--targetmasks={input.target_txt} --seedref={input.seed_res} --nsamples={config[''probtrack''][''nsamples'']} ' 
        '--dir={params.probtrack_dir} {params.probtrack_opts} -V 2  &> {log}'


rule transform_conn_to_template:
    input:
        connmap_3d = expand('diffparc/sub-{subject}/probtrack_{template}_{seed}_{hemi}/seeds_to_{target}.nii.gz',target=targets,allow_missing=True),
        affine =  lambda wildcards: glob(join(config['ants_template_dir'],config['ants_affine_mat'].format(subject=wildcards.subject))),
        warp =  lambda wildcards: glob(join(config['ants_template_dir'],config['ants_warp_nii'].format(subject=wildcards.subject))),
        ref = join(config['ants_template_dir'],config['ants_ref_nii'])
    output:
        connmap_3d = expand('diffparc/sub-{subject}/probtrack_{template}_{seed}_{hemi}/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True)
    envmodules: 'ants'
    singularity: config['singularity_neuroglia']
    threads: 32
    resources:
        mem_mb = 16000
    log: 'logs/transform_conn_to_template/sub-{subject}_{seed}_{hemi}_{template}.log'
    group: 'post_track'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation Linear -i {{1}} -o {{2}}  -r {input.ref} -t {input.warp} -t {input.affine} &> {log} :::  {input.connmap_3d} :::+ {output.connmap_3d}' 


rule save_connmap_template_npz:
    input:
        connmap_3d = expand('diffparc/sub-{subject}/probtrack_{template}_{seed}_{hemi}/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True),
        mask = 'diffparc/template_masks/sub-{template}_hemi-{hemi}_desc-{seed}_mask.nii.gz'
    output:
        connmap_npz = 'diffparc/sub-{subject}/connmap/sub-{subject}_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz'
    log: 'logs/save_connmap_to_template_npz/sub-{subject}_{seed}_{hemi}_{template}.log'
    group: 'post_track'
    script: 'scripts/save_connmap_template_npz.py'

rule gather_connmap_group:
    input:
        connmap_npz = expand('diffparc/sub-{subject}/connmap/sub-{subject}_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz',subject=subjects,allow_missing=True)
    output:
        connmap_group_npz = 'diffparc/connmap/group_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz'
    log: 'logs/gather_connmap_group/{seed}_{hemi}_{template}.log'
    run:
        import numpy as np
        
        #load first file to get shape
        data = np.load(input['connmap_npz'][0])
        affine = data['affine']
        mask = data['mask']
        conn_shape = data['conn'].shape
        nsubjects = len(input['connmap_npz'])
        conn_group = np.zeros([nsubjects,conn_shape[0],conn_shape[1]])
        
        for i,npz in enumerate(input['connmap_npz']):
            data = np.load(npz)
            conn_group[i,:,:] = data['conn']
            
        #save conn_group, mask and affine
        np.savez(output['connmap_group_npz'], conn_group=conn_group,mask=mask,affine=affine)
     
rule spectral_clustering:
    input:
        connmap_group_npz = 'diffparc/connmap/group_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz'
    params:
        max_k = config['max_k']
    output:
        cluster_k = expand('diffparc/clustering/group_space-{template}_seed-{seed}_hemi-{hemi}_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True)
    script: 'scripts/spectral_clustering.py'
        
     
    
