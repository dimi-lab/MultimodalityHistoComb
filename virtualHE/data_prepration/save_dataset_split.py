import os
from sklearn.model_selection import train_test_split
import shutil



he_input_dir = "/data_folder/Ovarian_TMA/virtual_HE/dataset/HE"  # Define an output directory
mxif_input_dir = "/data_folder/Ovarian_TMA/virtual_HE/dataset/MxIF"
output_dir = "/data_folder/Ovarian_TMA/virtual_HE/dataset/Img_split"  # Define an output directory

sub_dir = ["train", "test"]
sub_dir2 = ["A", "B"]
core_list = os.listdir(he_input_dir)
he_fn_list = []
for c in core_list:
    he_core_img_fn_list = os.listdir(os.path.join(he_input_dir, c))
    he_fn_list += he_core_img_fn_list

X_train, X_test = train_test_split(he_fn_list, test_size=0.2, random_state=42)

train_A_dir = os.path.join(output_dir, "trainA")
if not os.path.exists(train_A_dir):
    os.makedirs(train_A_dir)

train_B_dir = os.path.join(output_dir, "trainB")
if not os.path.exists(train_B_dir):
    os.makedirs(train_B_dir)

test_A_dir = os.path.join(output_dir, "testA")
if not os.path.exists(test_A_dir):
    os.makedirs(test_A_dir)

test_B_dir = os.path.join(output_dir, "testB")
if not os.path.exists(test_B_dir):
    os.makedirs(test_B_dir)

for x in X_train:
    core_id = x.split("_")[0]
    src_train_A = os.path.join(mxif_input_dir, core_id, x)
    dst_train_A = os.path.join(train_A_dir, x)
    shutil.copyfile(src_train_A, dst_train_A)

    src_train_B = os.path.join(he_input_dir, core_id, x)
    dst_train_B = os.path.join(train_B_dir, x)
    shutil.copyfile(src_train_B, dst_train_B)

for x in X_test:
    core_id = x.split("_")[0]
    src_test_A = os.path.join(mxif_input_dir, core_id, x)
    dst_test_A = os.path.join(test_A_dir, x)
    shutil.copyfile(src_test_A, dst_test_A)

    src_test_B = os.path.join(he_input_dir, core_id, x)
    dst_test_B = os.path.join(test_B_dir, x)
    shutil.copyfile(src_test_B, dst_test_B)


print("debug")