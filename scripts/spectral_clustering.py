import numpy as np
import pandas as pd
from sklearn.cluster import SpectralClustering

def save_label_nii (labelimg,affine,header,out_nifti):
    img = nib.Nifti1Image(labelimg,affine=affine,header=header)
    nib.save(img,out_nifti)

data = pd.read_csv(snakemake.input.correlation[0])

cluster_range = range(2,snakemake.params.max_k+1)
labels = np.zeros((data.shape[0],len(cluster_range)))

afile = snakemake.input.rois[0]
atlas = nib.load(afile)
atlas_data = atlas.get_fdata()

# pfile = snakemake.input.ptseries[0]
# p = nib.load(pfile)

# for ip in np.arange(0,n_parcels):
#     parcel = p.header.get_index_map(1)[ip+3]
#     name = parcel.name
#     if name == "BRAIN_STEM":
#         name = "ZIR"
#         zir_indices = parcel.voxel_indices_ijk

# Run spectral clustering
for i,k in enumerate(cluster_range):
    clustering = SpectralClustering(n_clusters=k, assign_labels="discretize",random_state=0,affinity='cosine').fit(abs(data.values))
    labels[:,i] = clustering.labels_
    
    labelimg = np.zeros(atlas_data.shape)
    for j in range(0,len(atlas_data[atlas_data==16])):
        labelimg[zir_indices[j][0],zir_indices[j][1],zir_indices[j][2]] = labels[j,i]
    save_label_nii(labelimg,atlas.affine,atlas.header,snakemake.output.niftis[i])

#df = pd.DataFrame(labels,columns=cluster_range)
#df.to_csv(snakemake.output.labels[0])