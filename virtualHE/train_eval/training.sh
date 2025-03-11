#!/usr/bin/bash

# export CUDA_VISIBLE_DEVICES=0
patch_root_dir=/data_folder/Ovarian_TMA/virtual_HE/dataset/ImgPair   # directory where saves image patch pairs for training

# train with 50 cases
training_case_txt_50=/data_folder/Ovarian_TMA/virtual_HE/dataset/training_cases.txt  # directory where saves WSI training cases
Train_output=/data_folder/Ovarian_TMA/virtual_HE/model_training  # Output directory to save trained models.

python pix2pix.py \
  --mode train \
  --output_dir ${Train_output} \
  --max_epochs 25 \
  --input_patch_root_dir ${patch_root_dir} \
  --input_case_list_txt ${training_case_txt_50} \
  --which_direction AtoB


