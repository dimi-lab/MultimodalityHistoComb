import os
import tifffile as tf
import imreg_dft as ird
from skimage import color
from skimage import io
import numpy as np
import cv2
import matplotlib.pyplot as plt
from path_config import my_path

HE_FOV_dir = my_path.HE_FOV_dir     # data directory to save HE FOV images
MxIF_FOV_dir = my_path.MxIF_FOV_dir  # data directory to save MxIF FOV images

FOV_id_list = ["F-10", "D-17", "J-26", "I-3", "H-1"]

#################################################################
HE_fn = os.path.join(HE_FOV_dir, FOV_id_list[0] + ".tif")
MxIF_fn = os.path.join(MxIF_FOV_dir, FOV_id_list[0] + ".tif")

HE_gray = color.rgb2gray(io.imread(HE_fn))
# plt.imshow(HE_gray, cmap="gray")
# plt.show()

MxIF_wsi_obj = tf.TiffFile(MxIF_fn)
MxIF_DAPI = MxIF_wsi_obj.pages[0].asarray().squeeze().astype(float)
# min(MxIF_DAPI.flatten())
# max(MxIF_DAPI.flatten())
MxIF_DAPI = (255 * (MxIF_DAPI - np.amin(MxIF_DAPI)) / np.ptp(MxIF_DAPI)).astype(float)
# MxIF_DAPI[MxIF_DAPI < 50] = 0
# MxIF_DAPI[MxIF_DAPI > 100] = 255
#################################################################
HE_gray = (255 * (HE_gray - np.amin(HE_gray)) / np.ptp(HE_gray)).astype(float)
# HE_gray[HE_gray < 50] = 0
# HE_gray[HE_gray > 100] = 255
Img_fix = HE_gray # H&E image
Img_float = np.zeros(Img_fix.shape) + 255 # create an array to save MxIF image
start = int((Img_fix.shape[0] - MxIF_DAPI.shape[0])/2)
Img_float[start: start+MxIF_DAPI.shape[0], start:start+MxIF_DAPI.shape[1]] = np.abs(MxIF_DAPI-255)

plt.imshow(Img_fix, cmap="gray")
plt.show()
plt.imshow(Img_float, cmap="gray")
plt.show()

con_s = dict(angle=[0, 0], scale=[1.3, 1.3])
sim = ird.similarity(Img_fix, Img_float, numiter=10, order=1, filter_pcorr=20, constraints=con_s)

new_img = ird.transform_img_dict(Img_float, sim)
plt.imshow(new_img, cmap='gray')
plt.show()


# plt.imshow(sim["timg"], cmap='gray')
# plt.show()


tvec = sim["tvec"].round(4)
angle = sim["angle"]
score = sim["success"]
offset = [tvec[1],tvec[0]]


# https://github.com/matejak/imreg_dft/blob/master/src/imreg_dft/imreg.py
# ird.similarity_matrix(sim['scale'], angle, vector)

print(offset)







