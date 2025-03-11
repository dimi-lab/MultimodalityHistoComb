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


# def convert_data_frame(df, metric, cols):
#     new_df = pd.DataFrame(columns=[metric, 'Value'])
#     group = ['Real_vs_MxIF', 'Real_vs_Ours', 'Real_vs_Vendor']
#     for index, row in df[cols].iterrows():
#         # row = next(df[cols].iterrows())[1]
#         print(row)
#         for idx, g in enumerate(group):
#             new_row = {metric: g, 'Value': row[cols[idx]]}
#             new_df.loc[len(new_df)] = new_row
#     return new_df

sns.set(rc={'figure.dpi': 250})
sns.set(font_scale=1.8)
# fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for idx, metric in enumerate(metrics):
    df_temp = pd.DataFrame(data=df, columns=metrics_group[idx])
    ax = sns.boxplot(data=df_temp, legend="full")
    # ax.set_xticklabels(['PSNR1', 'PSNR2', 'PSNR3'])
    # ax.set_xticklabels(['Real_vs_MxIF', 'Real_vs_Ours', 'Real_vs_Vendor'])
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles, ['Real_vs_MxIF', 'Real_vs_Ours', 'Real_vs_Vendor'], loc="best")
    # ax.legend(handles, ['r_vHE vs. MxIF', 'r_vHE vs. our_vHE', 'r_vHE vs. vendor_vHE'], loc="best")
    plt.title("%s by group" % metric)
    plt.show()





#
# sns.set(font_scale=1.5)
# # fig, axes = plt.subplots(1, 3, figsize=(15, 4))
#
# new_df = convert_data_frame(df, "PSNR", psnr_cols)
# sns.boxplot(data=new_df, x="PSNR", y="Value")

# ax.set_xticklabels(['PSNR1', 'PSNR2', 'PSNR3'])
# axes[idx].set_xticklabels(['Real_vs_MxIF', 'Real_vs_Our', 'Real_vs_Vendor'])
# axes[idx].set_xticklabels(['R_vHE/MxIF', 'R_vHE/O_vHE', 'R_vHE/V_vHE'])
# axes[0].set_xticklabels(['', '', ''])
# handles, _ = axes[idx].get_legend_handles_labels()
# axes[idx].legend(handles, ['Real_vs_MxIF', 'Real_vs_Ours', 'Real_vs_Vendor'], loc="best")
# ax.legend(handles, ['r_vHE vs. MxIF', 'r_vHE vs. our_vHE', 'r_vHE vs. vendor_vHE'], loc="best")
# axes[0].set_title("%s" % metrics[0])
#
# df_temp = pd.DataFrame(data=df, columns=metrics_group[1])
# sns.boxplot(data=df_temp, ax=axes[1], legend="full")
# # axes[1].set_xticklabels(['', '', ''])
# axes[1].set_title("%s" % metrics[1])
#
# df_temp = pd.DataFrame(data=df, columns=metrics_group[2])
# sns.boxplot(data=df_temp, ax=axes[2], legend="full")
# # axes[2].set_xticklabels(['', '', ''])
# axes[2].set_title("%s" % metrics[2])

plt.tight_layout()
plt.show()



#
# sns.set(rc={'figure.dpi': 250})
# sns.set(font_scale=1.5)
# fig, axes = plt.subplots(1, 3, figsize=(15, 4))
#
# df_temp = pd.DataFrame(data=df, columns=metrics_group[0])
#
#
#
# sns.boxplot(data=df_temp, ax=axes[0], legend="full")
# # ax.set_xticklabels(['PSNR1', 'PSNR2', 'PSNR3'])
# # axes[idx].set_xticklabels(['Real_vs_MxIF', 'Real_vs_Our', 'Real_vs_Vendor'])
# # axes[idx].set_xticklabels(['R_vHE/MxIF', 'R_vHE/O_vHE', 'R_vHE/V_vHE'])
# # axes[0].set_xticklabels(['', '', ''])
# # handles, _ = axes[idx].get_legend_handles_labels()
# # axes[idx].legend(handles, ['Real_vs_MxIF', 'Real_vs_Ours', 'Real_vs_Vendor'], loc="best")
# # ax.legend(handles, ['r_vHE vs. MxIF', 'r_vHE vs. our_vHE', 'r_vHE vs. vendor_vHE'], loc="best")
# axes[0].set_title("%s" % metrics[0])
#
# df_temp = pd.DataFrame(data=df, columns=metrics_group[1])
# sns.boxplot(data=df_temp, ax=axes[1], legend="full")
# # axes[1].set_xticklabels(['', '', ''])
# axes[1].set_title("%s" % metrics[1])
#
# df_temp = pd.DataFrame(data=df, columns=metrics_group[2])
# sns.boxplot(data=df_temp, ax=axes[2], legend="full")
# # axes[2].set_xticklabels(['', '', ''])
# axes[2].set_title("%s" % metrics[2])
#
# plt.tight_layout()
# plt.show()


print("DONE!")

