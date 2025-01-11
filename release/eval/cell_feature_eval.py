import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import pearsonr
from sklearn.preprocessing import normalize
import numpy as np
from scipy.stats import gaussian_kde

Show_plot_level = 1

data_dir = "/path/to/Ovarian_TMA/AlignmentEval/ApplyAlignment/AlignedCellsQuant"
results_dir = "/path/to/Ovarian_TMA/AlignmentEval/ApplyAlignment/AlignedCellsQuant_eval_out"

core_list_fn = "../core_list.txt"
core_list = open(core_list_fn, 'r').readlines()
roi_id_list = []
for core in core_list:
    roi_id = core.split(".")[0]
    roi_id_list.append(roi_id)

morp_features = ["Nucleus: Area µm^2", "Nucleus: Length µm", "Nucleus: Circularity", "Nucleus: Solidity",
                 "Nucleus: Max diameter µm", "Nucleus: Min diameter µm", "Cell: Area µm^2", "Cell: Length µm",
                 "Cell: Circularity", "Cell: Solidity", "Cell: Max diameter µm", "Cell: Min diameter µm",
                 "Nucleus/Cell area ratio"]

HE_H_stain_features = ["Hematoxylin: Nucleus: Mean", "Hematoxylin: Nucleus: Median", "Hematoxylin: Nucleus: Min",
                     "Hematoxylin: Nucleus: Max", "Hematoxylin: Nucleus: Std.Dev."]

MxIF_DAPI_stain_features = ["DAPI_AF_R01: Nucleus: Mean", "DAPI_AF_R01: Nucleus: Median", "DAPI_AF_R01: Nucleus: Min",
                            "DAPI_AF_R01: Nucleus: Max", "DAPI_AF_R01: Nucleus: Std.Dev."]

matplotlib.rcParams.update({'font.size': 12})
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)

morph_r_p_value = np.zeros((len(morp_features), len(roi_id_list), 2))
stain_r_p_value = np.zeros((len(HE_H_stain_features), len(roi_id_list), 2))

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
def KDE_density(x, y):
    xy = np.vstack([x, y])
    kde_scores = gaussian_kde(xy)(xy)
    norm_kde_scores = norm(kde_scores, 0, 1)
    return norm_kde_scores

for idx, roi in enumerate(roi_id_list):
    # if idx > 2:
    #     break
    MxIF_fn = os.path.join(data_dir, roi + "_1on1_MxIF_quant.tsv")
    HE_fn = os.path.join(data_dir, roi + "_1on1_HE_quant.tsv")

    MxIF_df = pd.read_table(MxIF_fn, sep='\t', header=0)
    HE_df = pd.read_table(HE_fn, sep='\t', header=0)

    for m_idx in range(len(morp_features)):
        he_mf = HE_df[morp_features[m_idx]]
        mxif_mf = MxIF_df[morp_features[m_idx]]

        # r_m = pearsonr(normalize([list(he_mf)]), normalize([list(mxif_mf)]))
        r_m = pearsonr(normalize([list(he_mf)]).flatten().T, normalize([list(mxif_mf)]).flatten().T)
        morph_r_p_value[m_idx, idx, :] = [r_m.correlation, r_m.pvalue]

        if Show_plot_level == 0:
            fig = plt.figure(figsize=(5, 4), dpi=300)
            plt.title(roi)
            plt.plot(he_mf, mxif_mf, ".")
            y = max(list(mxif_mf)) * 0.9
            x = min(list(he_mf)) * 1.1
            if r_m.pvalue < 0.001:
                p_str = "p<0.001"
            elif r_m.pvalue > 0.001 and r_m.pvalue < 0.05:
                p_str = "p<0.05"
            else:
                p_str = "p=%.4f" % r_m.pvalue
            plt.text(x, y, "r=%.4f\n%s" % (r_m.correlation, p_str))
            plt.xlabel("H&E %s" % morp_features[m_idx])
            plt.ylabel("MxIF %s" % morp_features[m_idx])
            plt.tight_layout()
            plt.show()

    for s_idx in range(len(HE_H_stain_features)):
        he_sf = HE_df[HE_H_stain_features[s_idx]]
        mxif_sf = MxIF_df[MxIF_DAPI_stain_features[s_idx]]

        # r_s = pearsonr(normalize(he_sf), normalize(mxif_sf))
        r_s = pearsonr(normalize([list(he_sf)]).flatten().T, normalize([list(mxif_sf)]).flatten().T)
        stain_r_p_value[s_idx, idx, :] = [r_s.correlation, r_s.pvalue]

        if Show_plot_level == 0:
            fig2 = plt.figure(figsize=(5, 4), dpi=300)
            plt.plot(he_sf, mxif_sf, ".")
            plt.title(roi)
            y = max(list(mxif_sf)) * 0.9
            x = min(list(he_sf)) * 1.1
            if r_s.pvalue < 0.001:
                p_str = "p<0.001"
            elif r_s.pvalue > 0.001 and r_s.pvalue < 0.05:
                p_str = "p<0.05"
            else:
                p_str = "p=%.4f" % r_s.pvalue
            plt.text(x, y, "r=%.4f\n%s" % (r_s.correlation, p_str))
            plt.xlabel("H&E %s" % HE_H_stain_features[s_idx])
            plt.ylabel("MxIF %s" % MxIF_DAPI_stain_features[s_idx])
            plt.tight_layout()
            plt.show()

r = morph_r_p_value[:, :, 0]
p = morph_r_p_value[:, :, 1]


total_HE_df_list = []
total_MxIF_df_list = []
for idx, roi in enumerate(roi_id_list):
    # if idx > 2:
    #     break
    MxIF_fn = os.path.join(data_dir, roi + "_1on1_MxIF_quant.tsv")
    HE_fn = os.path.join(data_dir, roi + "_1on1_HE_quant.tsv")

    MxIF_df = pd.read_table(MxIF_fn, sep='\t', header=0)
    HE_df = pd.read_table(HE_fn, sep='\t', header=0)

    total_HE_df_list.append(HE_df)
    total_MxIF_df_list.append(MxIF_df)

total_HE_df = pd.concat(total_HE_df_list)
total_MxIF_df = pd.concat(total_MxIF_df_list)

for m_idx in range(len(morp_features)):
    he_mf = total_HE_df[morp_features[m_idx]]
    mxif_mf = total_MxIF_df[morp_features[m_idx]]

    # r_m = pearsonr(he_mf, mxif_mf)
    r_m = pearsonr(normalize([list(he_mf)]).flatten().T, normalize([list(mxif_mf)]).flatten().T)
    if Show_plot_level == 1:
    # if Show_plot_level == 1 and m_idx == 0:
        norm_kde_scores = KDE_density(list(he_mf), list(mxif_mf))
        fig = plt.figure(figsize=(5, 4), dpi=300)
        # plt.plot(he_mf, mxif_mf, ".")
        plt.scatter(he_mf, mxif_mf, c=norm_kde_scores, marker=".", cmap='jet')
        y = max(list(mxif_mf)) * 0.9
        x = min(list(he_mf)) * 1.1
        if r_m.pvalue < 0.001:
            p_str = "p<0.001"
        elif r_m.pvalue > 0.001 and r_m.pvalue < 0.05:
            p_str = "p<0.05"
        else:
            p_str = "p=%.4f" % r_m.pvalue
        plt.text(x, y, "r=%.4f\n%s" % (r_m.correlation, p_str))
        ele = morp_features[m_idx].split(":")
        if len(ele) == 2:
            x_label = morp_features[m_idx].split(":")[0]
            y_label = morp_features[m_idx].split(":")[0]
            title = morp_features[m_idx].split(":")[1]
        else:
            x_label = morp_features[m_idx].replace("/", " ")
            y_label = morp_features[m_idx].replace("/", " ")
            title = morp_features[m_idx].replace("/", " ")
        plt.title(title)
        plt.xlabel("H&E %s" % x_label)
        plt.ylabel("MxIF %s" % y_label)
        # plt.colorbar()
        plt.tight_layout()
        fn = os.path.join(results_dir, "morph_feature_corr_%s.png" % title)
        plt.savefig(fn)
        # plt.show()

for s_idx in range(len(HE_H_stain_features)):
    he_sf = total_HE_df[HE_H_stain_features[s_idx]]
    mxif_sf = total_MxIF_df[MxIF_DAPI_stain_features[s_idx]]
    # r_s = pearsonr(he_sf, mxif_sf)
    r_s = pearsonr(normalize([list(he_sf)]).flatten().T, normalize([list(mxif_sf)]).flatten().T)
    if Show_plot_level == 1:
    # if Show_plot_level == 1 and s_idx == 0:
        norm_kde_scores = KDE_density(list(he_sf), list(mxif_sf))
        fig2 = plt.figure(figsize=(5, 4), dpi=300)
        # plt.plot(he_sf, mxif_sf, ".")
        plt.scatter(he_sf, mxif_sf, c=norm_kde_scores, marker=".", cmap='jet')
        y = max(list(mxif_sf)) * 0.9
        x = min(list(he_sf)) * 1.1
        if r_s.pvalue < 0.001:
            p_str = "p<0.001"
        elif r_s.pvalue > 0.001 and r_s.pvalue < 0.05:
            p_str = "p<0.05"
        else:
            p_str = "p=%.4f" % r_s.pvalue
        plt.text(x, y, "r=%.4f\n%s" % (r_s.correlation, p_str))
        x_label = HE_H_stain_features[s_idx].split(":")[0]
        y_label = MxIF_DAPI_stain_features[s_idx].split(":")[0]
        title = HE_H_stain_features[s_idx].split(":")[1] + ":" + HE_H_stain_features[s_idx].split(":")[2]
        plt.title(title)
        plt.xlabel("H&E %s" % x_label)
        plt.ylabel("MxIF %s" % "DAPI")
        # plt.colorbar()
        plt.tight_layout()
        fn = os.path.join(results_dir, "stain_feature_corr_%s.png" % title)
        plt.savefig(fn)
        # plt.show()



print("debug")

print("Done")
