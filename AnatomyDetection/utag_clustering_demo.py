from utag import utag
# Use Scanpy to get a h5ad file with provided data
import scanpy as sc

# In[]
adata = sc.read(
    'data/healthy_lung_adata.h5ad',
    backup_url='https://zenodo.org/record/6376767/files/healthy_lung_adata.h5ad?download=1')
# In[]
# Run UTAG on provided data
utag_results = utag(
    adata,
    slide_key="roi",
    max_dist=20,
    normalization_mode='l1_norm',
    apply_clustering=True,
    clustering_method='leiden',
    resolutions=[0.05, 0.1, 0.3]
)
# utag_results = utag(
#     adata,
#     slide_key=None,
#     max_dist=20,
#     normalization_mode='l1_norm',
#     apply_clustering=True,
#     clustering_method='leiden',
#     resolutions=[0.05, 0.1, 0.3]
# )
# In[]
obs = adata.obs
var = adata.var

obs_a = adata.obs[['roi', 'X_centroid', 'Y_centroid']]

adata.obsm["spatial"]
# In[]
for roi in utag_results.obs['roi'].unique():
    result = utag_results[utag_results.obs['roi'] == roi].copy()
    sc.pl.spatial(result, color='UTAG Label_leiden_0.1', spot_size=10)
    print("d")
