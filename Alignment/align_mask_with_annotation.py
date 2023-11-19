import os
import tifffile as tf
import imreg_dft as ird
from skimage import color
from skimage import io
import numpy as np
import cv2
import matplotlib.pyplot as plt
from path_config import my_path


def load_annotation():
    pass

def calculate_affine():
    pass

def apply_affine():
    pass


if __name__ == "__main__":
    HE_export_dir = my_path.HE_FOV_export_dir     # data directory to save HE FOV images
    MxIF_export_dir = my_path.MxIF_FOV_export_dir  # data directory to save MxIF FOV images

    FOV_id_list = ["A-8", "A-9", "A-12", "A-14", "B-24"]

    for fov in FOV_id_list:
        he_dir = os.path.join(HE_export_dir, fov, "*-mask.png")
        mxif_dir = os.path.join(MxIF_export_dir, fov, "*-mask.png")


















