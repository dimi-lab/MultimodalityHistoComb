import os
from sklearn.model_selection import train_test_split


input_dir = "/data_folder/Ovarian_TMA/virtual_HE/dataset/ImgPair"  # Define an output directory
output_dir = "/data_folder/Ovarian_TMA/virtual_HE/dataset/ImgPairs_split"  # Define an output directory

sub_dir = ["train", "val", "test"]
core_list = os.listdir(input_dir)
all_fn_list = []
for c in core_list:
    core_img_fn_list = os.listdir(os.path.join(input_dir, c))
    all_fn_list += core_img_fn_list

X_temp, X_val = train_test_split(all_fn_list, test_size=0.2, random_state=42)
X_train, X_test = train_test_split(X_temp, test_size=0.2, random_state=42)




print("debug")