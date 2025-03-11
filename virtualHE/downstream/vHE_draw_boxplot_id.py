import os
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


csv_fn = "/data_folder/Ovarian_TMA/virtual_HE/pytorch_model_testing/vHE_pix2pix/metrics/psnr_ssim_vif_all.csv"
sv_dir = "/data_folder/Ovarian_TMA/virtual_HE/pytorch_model_testing/vHE_pix2pix/metrics"

df = pd.read_csv(csv_fn)

psnr_cols = ["PSNR", "PSNR_r", "PSNR_rr"]
ssim_cols = ["SSIM", "SSIM_r", "SSIM_rr"]
vif_cols = ["VIF", "VIF_r",	"VIF_rr"]

metrics = ['PSNR', 'SSIM', 'VIF']
metrics_group = [psnr_cols, ssim_cols, vif_cols]

sns.set(rc={'figure.dpi': 250})
sns.set(font_scale=1.8)
for idx, metric in enumerate(metrics):
    df_temp = pd.DataFrame(data=df, columns=metrics_group[idx])
    ax = sns.boxplot(data=df_temp, legend="full")
    # ax.set_xticklabels(['PSNR1', 'PSNR2', 'PSNR3'])
    # ax.set_xticklabels(['R_vHE/MxIF', 'R_vHE/O_vHE', 'R_vHE/V_vHE'])
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles, ['R_vHE/MxIF', 'R_vHE/O_vHE', 'R_vHE/V_vHE'], loc="best")
    # ax.legend(handles, ['r_vHE vs. MxIF', 'r_vHE vs. our_vHE', 'r_vHE vs. vendor_vHE'], loc="best")
    plt.title("%s by group" % metric)
    plt.show()


print("DONE!")

