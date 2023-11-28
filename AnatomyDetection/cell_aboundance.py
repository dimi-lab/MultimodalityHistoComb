from natsort import natsorted
import os
import pandas as pd
from run_UTAG import get_specific_feature
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import seaborn as sns


def norm(rvalue, newmin, newmax):
    oldmin = min(rvalue)
    oldmax = max(rvalue)
    oldrange = oldmax - oldmin
    newrange = newmax - newmin
    if oldrange == 0:  # Deal with the case where rvalue is constant:
        if oldmin < newmin:  # If rvalue < newmin, set all rvalue values to newmin
            newval = newmin
        elif oldmin > newmax:  # If rvalue > newmax, set all rvalue values to newmax
            newval = newmax
        else:  # If newmin <= rvalue <= newmax, keep rvalue the same
            newval = oldmin
        normal = [newval for _ in rvalue]
    else:
        scale = newrange / oldrange
        normal = [(v - oldmin) * scale + newmin for v in rvalue]
    return normal

def KDE_cell_density(cell_xy):
    xy = np.vstack([cell_xy[:, 0], cell_xy[:, 1]])
    kde_scores = gaussian_kde(xy)(xy)
    norm_kde_scores = norm(kde_scores, 0, 1)
    return norm_kde_scores


if __name__ == "__main__":
    CheckCellsWithNullFeatures = True
    Debug = True
    # output_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/output/utag_clustering/MxIF"
    # data_root_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/OV_TMA_MxIF_QuPath_proj_AB/export"
    output_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/output/utag_clustering/HE"
    data_root_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/OV_TMA_HE_QuPath_proj_AB/export"
    roi_list = natsorted(os.listdir(data_root_dir))
    for idx, roi_id in enumerate(roi_list):
        # roi_id = 'A-1'
        print("Processing: %s, current progress: %d/%d" % (roi_id, idx+1, len(roi_list)))
        csv_fn = os.path.join(data_root_dir, roi_id, roi_id + "_QUANT.tsv")
        org_df = pd.read_csv(csv_fn, sep='\t')

        col = list(org_df.columns)
        df = org_df.dropna(subset=["Image"] + col[6:])
        xy_names, cell_loc = get_specific_feature(df, 6, 8)
        mf_names, morphological_features = get_specific_feature(df, 9, 22)
        pf_names, protein_features = get_specific_feature(df, 22)
        print("Show cells with null features")
        if CheckCellsWithNullFeatures:
            row_idx = org_df.iloc[:, 6:].isnull().any(axis=1)
            df_cell_with_Null_features = org_df.loc[row_idx, :]
            print("Cell counts with missing features: %d" % df_cell_with_Null_features.shape[0])
            if df_cell_with_Null_features.shape[0] > 0 and Debug:
                plt.scatter(cell_loc[:, 0], cell_loc[:, 1], marker=".")
                null_xy_names, null_cell_loc = get_specific_feature(df_cell_with_Null_features, 6, 8)
                plt.scatter(null_cell_loc[:, 0], null_cell_loc[:, 1], marker="*")
                plt.legend(["All cells", 'Cells with Null features'])
                plt.axis('equal')
                plt.gca().invert_yaxis()
                plt.title("%d Cell(s) with null features in %s" % (df_cell_with_Null_features.shape[0], roi_id))
                plt.show()

        all_features = np.concatenate([morphological_features, protein_features], axis=1)
        all_feature_names = list(mf_names) + list(pf_names)

        norm_kde_scores = KDE_cell_density(cell_loc)

        plt.scatter(cell_loc[:, 0], cell_loc[:, 1], c=norm_kde_scores, marker=".", cmap='jet')
        plt.grid()
        plt.title("cell densities")
        # plt.xlim([0, 10])
        # plt.ylim([2, 12])
        plt.colorbar()
        plt.axis('equal')
        plt.gca().invert_yaxis()

        sns.kdeplot(x=cell_loc[:, 0], y=cell_loc[:, 1])

        # plt.show()
        save_to = os.path.join(output_dir, roi_id +"_cell_density.png")
        plt.savefig(save_to)
        plt.close()


        print("###########################################################")





