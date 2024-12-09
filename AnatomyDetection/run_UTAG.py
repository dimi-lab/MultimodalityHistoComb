import os

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.sparse import csr_matrix
from utag import utag
print(ad.__version__)


def get_specific_feature(df, idx_start, idx_end=None):
    if idx_end is None:
        val = np.array(df.iloc[:, idx_start:])
        names = df.iloc[:, idx_start:].columns
    else:
        val = np.array(df.iloc[:, idx_start:idx_end])
        names = df.iloc[:, idx_start:idx_end].columns
    return names, val


def visualize_rois(adata, n_rois=3, roi_key='roi', color_key='cell type', dpi=200, save_to=None):
    # increase dpi for a more crisp visualization
    fig, axes = plt.subplots(1, n_rois, figsize=(n_rois * 3 + 1, 3), dpi=dpi)

    for i, roi in enumerate(adata.obs[roi_key].unique()[:n_rois]):
        a = adata[adata.obs[roi_key] == roi].copy()
        if n_rois == 1:
            sc.pl.spatial(a, color=color_key, spot_size=10, ax=axes, title='', frameon=False, show=False)
        else:
            if i == n_rois - 1:
                sc.pl.spatial(a, color=color_key, spot_size=10, ax=axes[i], title='', frameon=False, show=False)
            else:
                sc.pl.spatial(a, color=color_key, spot_size=10, ax=axes[i], title='', frameon=False, legend_loc='none',
                              show=False)
    plt.tight_layout()
    if save_to is None:
        plt.show()
    else:
        plt.savefig(save_to)

if __name__ == "__main__":
    output_dir = "/temp/Ovarian_TMA/output/utag_clustering/MxIF"
    data_root_dir = "/temp/Ovarian_TMA/OV_TMA_MxIF_QuPath_proj_AB/export"
    # output_dir = "/temp/Ovarian_TMA/output/utag_clustering/HE"
    # data_root_dir = "/temp/Ovarian_TMA/OV_TMA_HE_QuPath_proj_AB/export"
    roi_list = sorted(os.listdir(data_root_dir))

    resolutions_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for idx, roi_id in enumerate(roi_list):
        print("Processing %s, finished %d/%d" % (roi_id, idx, len(roi_list)))
        csv_fn = os.path.join(data_root_dir, roi_id, roi_id+"_QUANT.tsv")
        org_df = pd.read_csv(csv_fn, sep='\t')
        col = list(org_df.columns)
        df = org_df.dropna(subset=["Image"] + col[6:])
        xy_names, cell_loc = get_specific_feature(df, 6, 8)
        mf_names, morphological_features = get_specific_feature(df, 9, 22)
        pf_names, protein_features = get_specific_feature(df, 22)
        print("debug")

        all_features = np.concatenate([morphological_features, protein_features], axis=1)
        all_feature_names = list(mf_names) + list(pf_names)

        obs = pd.concat([df['Image'], df[xy_names[0]], df[xy_names[1]]], axis=1)

        ov_tma_core_adata = ad.AnnData(
            X=all_features.astype(np.float64),
            obs=obs
        )
        ov_tma_core_adata.obsm['spatial'] = cell_loc



        org_save_to = os.path.join(output_dir, roi_id+"_org.png")
        visualize_rois(ov_tma_core_adata, n_rois=1, roi_key='Image', color_key='Image', save_to=org_save_to)

        utag_results = utag(
            ov_tma_core_adata,
            slide_key=None,
            max_dist=20,
            normalization_mode='l1_norm',
            apply_clustering=True,
            clustering_method='leiden',
            resolutions=resolutions_list
        )
        utag_csv_save_to = os.path.join(output_dir, roi_id + "_utag_csv")
        utag_results.write_csvs(utag_csv_save_to)

        for res in resolutions_list:
            utag_img_save_to = os.path.join(output_dir, roi_id + "_utag_csv", roi_id+"_utag_leiden_{:.1f}.png".format(res))
            visualize_rois(utag_results, n_rois=1, roi_key='Image', color_key='UTAG Label_leiden_{:.1f}'.format(res), save_to=utag_img_save_to)



