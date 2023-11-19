import os
import glob
import pandas as pd
import numpy as np
import cv2
from aicspylibczi import CziFile
from pathlib import Path
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
import tifffile as tf
from openslide import open_slide

DEBUG = True
########################
# load human annotated landmarks
HE_landmarks = [[68835.5, 17866.8], [33810.8, 35933.6], [104541.3, 44922.6]]
MxIF_landmarks = [[46675.1, 11409.4], [22632.0, 23613.8], [70447.5, 30261.4]]

# calculate affine transformation
M = cv2.getAffineTransform(np.array(HE_landmarks).astype(np.float32), np.array(MxIF_landmarks).astype(np.float32))

def get_transformed_coord(M, src_x, src_y):
    '''
    Apply the affine transformation to a set of source coordinates
     to get the corresponding coordinate in destination coordinate.
    :param M: An affine transformation matrix from cv2.getAffineTransform, with shape of 2x2
    :param xy:
    :return:
    '''
    # src_x = xy[:, 0]
    # src_y = xy[:, 1]
    xy = np.column_stack((np.transpose(src_x), np.transpose(src_y)))
    MM = np.vstack((M, np.array([0, 0, 1])))  # extend into 3x3 matrix
    new_XY = np.vstack((np.transpose(xy), np.ones((1, len(xy)))))

    transformed_loc = np.dot(MM, new_XY)
    transformed_xy = np.transpose(transformed_loc[0:2, :])
    return transformed_xy

##################################################
# FOV locations from QuPath dearrayer
csv_dir = r"\\mfad\researchmn\HCPR\HCPR-GYNECOLOGICALTUMORMICROENVIRONMENT\Archive\WSIs\Ovarian_TMA"
csv_fn = os.path.join(csv_dir, "Rab_Spine-22R919-A-SERBSVG-5X-08.czi_FOV_loc.csv")
df = pd.read_csv(csv_fn)
src_x = np.array(df["X"])
src_y = np.array(df["Y"])

FOV_W = int(list(df["W"])[0])
FOV_H = int(list(df["H"])[0])
FOV_id_list = list(df["FOV"])


trans_xy = get_transformed_coord(M, src_x, src_y)

wrt_str = "FOV,X,Y,W,H,Missing\n"
MxIF_df = df.copy(deep=True)
MxIF_df["X"] = trans_xy[:, 0]
MxIF_df["Y"] = trans_xy[:, 1]
csv_fn = os.path.join(csv_dir, "Slide2050_24Plex.ome.tif - SLIDE-2050_c1.tif_FOV_loc_affine.csv")
MxIF_df.to_csv(csv_fn, sep=",", index=False)
















