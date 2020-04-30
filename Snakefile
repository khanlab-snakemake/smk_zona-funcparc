from os.path import join
from glob import glob
import pandas as pd

configfile: 'config.yml'

#load participants.tsv file (list of HCP subject IDs),
df = pd.read_table(config['participants_tsv'])
subjects = df.participant_id.to_list() 

wildcard_constraints:
    subject="[0-9]+"

rule cifti_separate:
    input: join(config['hcp1200_dir'],'/{subject}/MNINonLinear/Results/rfMRI_REST2_7T_AP/rfMRI_REST2_7T_AP_Atlas_1.6mm.dtseries.nii')
    output: 
	lh: 'funcparc/{subject}/input/rfMRI_REST2_7T_AP.L.59k_fs_LR.func.gii'
	rh: 'funcparc/{subject}/input/rfMRI_REST2_7T_AP.R.59k_fs_LR.func.gii'
    singularity: config['singularity_connectomewb']
    log: 'logs/cifti_separate/{subject}.log'
    shell:
        'wb_command -cifti-separate {input} COLUMN -metric CORTEX_LEFT {output.lh} -metric CORTEX_RIGHT {output.rh}'
