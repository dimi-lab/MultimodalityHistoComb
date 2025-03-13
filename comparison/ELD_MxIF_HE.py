# Refer to example code within ELD
# https://github.com/ekvall93/ELD/blob/main/notebooks/multimodal.ipynb
# Pretrained models were downloaded from https://figshare.com/articles/dataset/models/23059625?file=40799408
# Multimodel->SMA->model_10.fan.pth
# Tested on Python=3.11.11
# Other dependencies are described in ELD original repo.
# 
# Please Note: This method mainly for key points detection. 
# Stricky parts in preprocessing steps:
# 1. mask_background does not work if the H&E image arrary is saved in float data type. for MxIF images

import os
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"]="1"

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from skimage.io import imread, imshow
import tifffile as tf
import xml
import cv2
from glob import glob
import pandas as pd
import json
import os, io
import numpy as np

import torch
import seaborn as sns
from ELD.utils import (toImg, preprocess, predict_landmarks, create_target_landmarks, 
                       create_target_images, download_images_urls, downscale_images, plot_images, 
                       mask_background, padImg, crop_non_tissue, downsize_and_save, 
                       rescale_landmarks, pad_image_and_adjust_landmarks, corr,plot_warped_images)
from ELD.model import loadFan, crop, toGrey
from ELD.warp import Homo, Rigid, TPS







def get_tiff_channel_names(fn):
    ome_tiff = tf.TiffFile(fn)
    # option 1
    # omexml_string = ome_tiff.pages[0].description
    # root = xml.etree.ElementTree.parse(io.StringIO(omexml_string))
    # namespaces = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
    # channels = root.findall('ome:Image[1]/ome:Pixels/ome:Channel', namespaces)
    # channel_names = [c.attrib['Name'] for c in channels]
    # option 2
    channel_names = ome_tiff.imagej_metadata.get("Labels")

    return channel_names

def get_channel_idx(channel_name_list, channel_query_list):
    channel_index = []
    for c in channel_query_list:
        for idx, i in enumerate(channel_name_list):
            if i.lower() == c.lower():
                channel_index.append(idx)
    return channel_index

def read_img_arr(img_fn):
    ext = os.path.splitext(img_fn)[1]
    if ext == ".tif" or ext == ".tiff":
        img_arr = tf.imread(img_fn)
        return img_arr
    else:
        print("Currently, doesn't support this format")
        return None

def get_mxif_img(mxif_wsi_fn, selected_channels):
    channel_names = get_tiff_channel_names(mxif_wsi_fn)
    selected_channel_idx = get_channel_idx(channel_names, selected_channels)
    assert len(selected_channel_idx) == 3
    mxif_img_arr = read_img_arr(mxif_wsi_fn)
    mxif_img_arr = np.swapaxes(mxif_img_arr, 0, 2)
    mxif_img_arr = np.swapaxes(mxif_img_arr, 0, 1)[:, :, selected_channel_idx]
    # mxif_img_arr[np.all(mxif_img_arr == 0, axis=-1)] = [255, 255, 255]
    mxif_img_arr[:, :, 1] = mxif_img_arr[:, :, 0] + mxif_img_arr[:, :, 1]
    mxif_img_arr[:, :, 2] = mxif_img_arr[:, :, 0] + mxif_img_arr[:, :, 2]
    return mxif_img_arr[:, :, 0:3]

def normalize_img_arr(img_arr):
    normalized_img_arr = np.zeros_like(img_arr)
    for i in range(img_arr.shape[-1]):
        channel = img_arr[:, :, i]
        min_val = np.min(channel)
        max_val = np.max(channel)
        normalized_img_arr[:, :, i] = 255 * (channel - min_val) / (max_val - min_val)
    return normalized_img_arr.astype(np.uint8)


if __name__ == "__main__":

    # HE_img_fn = "/Users/jjiang10/Data/OV_TMA/HE_H-13.tif"
    # MxIF_img_fn = "/Users/jjiang10/Data/OV_TMA/MxIF_H-13.ome.tif"

    HE_img_fn = "/data/jjiang10/Data/OV_TMA/HE_H-13.tif"
    MxIF_img_fn = "/data/jjiang10/Data/OV_TMA/MxIF_H-13.ome.tif"

    selected_channels = ["DAPI_AF_R01", "NAK", "PANCK"]

    he_img = tf.TiffFile(HE_img_fn).pages[0].asarray()  # read the images
    # mxif_img = get_mxif_img(MxIF_img_fn, selected_channels)
    mxif_img_arr = read_img_arr(MxIF_img_fn)[0:3,:, :]
    mxif_img_arr = np.swapaxes(mxif_img_arr, 0, 2)
    mxif_img = np.swapaxes(mxif_img_arr, 0, 1)
    mxif_img_norm = normalize_img_arr(mxif_img)


    imageList = []
    imageList.append(he_img)
    imageList.append(mxif_img_norm)
    # imageList.append(he_img)

    print(len(imageList))
    imageList = downscale_images(imageList)
    imageList = mask_background(imageList)
    imageList = crop_non_tissue(imageList)
    small_imgs = downsize_and_save(imageList, "/data/jjiang10/Data/ELD_test/")

    '''
    Train the model
    ``` bash
    python -m visdom.server -port 9006

    eld-train --elastic_sigma 5 --cuda 1 --port 9006 --data_path /data/ekvall/tutorial/ --npts 14 --o scratch --step_size 5 --ws 0 --gamma 0.9 --angle 8 --model unimodal
    ``` 
    '''


    image = torch.stack([preprocess(img) for img in small_imgs])

    fan = loadFan(npoints=14,n_channels=3,path_to_model="../Exp_2/model_158.fan.pth")
    #predict landmarks
    pts = predict_landmarks(fan, image)

    #combine landmarks and image
    np_img = toImg(image.cuda()[:,:3], pts, 128)

    fig, axs = plt.subplots(3, 4, figsize=(15, 10))  # adjust the size as needed
    axs = axs.ravel()

    for i in range(len(np_img)):
        img = np_img[i]
        axs[i].imshow(img)
        axs[i].set_title(f"Image {i+1}")
        axs[i].axis('off')  # to hide the axis

    plt.tight_layout()
    plt.show()

    scaled_pts = rescale_landmarks(pts, imageList)
    padded_images_torch, adjusted_landmarks = pad_image_and_adjust_landmarks(imageList, scaled_pts)

    np_img = toImg(padded_images_torch.cuda()[:,:3], adjusted_landmarks, 5 * 128)

    fig, axs = plt.subplots(3, 4, figsize=(15, 10))  # adjust the size as needed
    axs = axs.ravel()

    for i in range(len(np_img)):
        img = np_img[i]
        axs[i].imshow(img)
        axs[i].set_title(f"Image {i+1}")
        axs[i].axis('off')  # to hide the axis

    plt.tight_layout()
    plt.show()

    image = padded_images_torch
    dst_image = create_target_images(image, 0)

    pts = adjusted_landmarks
    dst_pts = create_target_landmarks(pts, 0)


    # Rigid transformation
    rigid_transform = Rigid()
    #warp images
    mapped_imgs = rigid_transform.warp_img(image.cuda(), pts, dst_pts, (863, 863))
    #warp landmarks
    mapped_pts = rigid_transform.warp_pts(pts, dst_pts, pts)

    #rigid_loss = corr(*crop(mapped_imgs, dst_image.cuda())).cpu().numpy()[1:]
    rigid_loss = corr(mapped_imgs, dst_image.cuda()).cpu().numpy()[1:]

    plot_warped_images(mapped_imgs, mapped_pts, rigid_loss, 5 * 128, 'Rigid')


