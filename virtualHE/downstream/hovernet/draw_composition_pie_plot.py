import seaborn as sns
import matplotlib.pyplot as plt
import os
import pandas as pd

# tsv_dir = "/data_folder/Ovarian_TMA/virtual_HE/dataset/Img_split/testB_pred/qupath"
tsv_root_dir = "/data_folder/Ovarian_TMA/virtual_HE/pytorch_model_testing/vHE_pix2pix/test_latest"


models = [ "MoNuSAC", "CoNSeP", "PanNuke"]

def merge_csvs(folder_path, postfix=".tsv"):
    dataframes = []
    for file_name in os.listdir(folder_path):
        if file_name.endswith(postfix):
            file_path = os.path.join(folder_path, file_name)
            df = pd.read_csv(file_path, sep='\t')
            dataframes.append(df)
    merged_df = pd.concat(dataframes, ignore_index=True)
    return merged_df


sns.set(rc={'figure.dpi': 250})
sns.set(font_scale=2)

for m in models:
    tsv_dir = os.path.join(tsv_root_dir, "images_for_HoverNet_%s_pred" % m.lower(), "qupath")

    merged_df_real = merge_csvs(tsv_dir, "real_B.tsv")
    merged_df_fake = merge_csvs(tsv_dir, "fake_B.tsv")

    # fig, axes = plt.subplots(2, 1, figsize=(5, 10))
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    sns.set_style("whitegrid")
    merged_df_real.replace('necros', 'other', inplace=True)
    merged_df_real.replace('connec', 'other', inplace=True)
    category_counts = merged_df_real['name'].value_counts()
    axes[0].pie(category_counts, labels=category_counts.index, autopct='%1.1f%%')
    axes[0].set_title('Real H&E ')

    merged_df_fake.replace('necros', 'other', inplace=True)
    merged_df_fake.replace('connec', 'other', inplace=True)
    category_counts = merged_df_fake['name'].value_counts()
    axes[1].pie(category_counts, labels=category_counts.index, autopct='%1.1f%%')
    axes[1].set_title('Virtual H&E')

    # fig.suptitle('Cell detection model: %s' %m)
    fig.suptitle("Cell Composition")
    plt.tight_layout()
    plt.show()

print("DONE!")
