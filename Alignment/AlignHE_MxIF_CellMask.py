import os
import imreg_dft as ird
from skimage import color
from skimage import io
import numpy as np
import cv2
from PIL import Image
import matplotlib.pyplot as plt
from path_config import my_path


if __name__ == "__main__":
    HE_export_dir = my_path.HE_FOV_export_dir     # data directory to save HE FOV images
    MxIF_export_dir = my_path.MxIF_FOV_export_dir  # data directory to save MxIF FOV images

    fov = "A-22"
    he_cell_mask_fn = os.path.join(HE_export_dir, fov, "A-22_(1.00,0,0,4108,4108)-mask.png")
    mxif_cell_mask_fn = os.path.join(MxIF_export_dir, fov, "A-22_(1.00,0,0,3000,3000)-mask.png")

    Img_fix = np.array(Image.open(he_cell_mask_fn))
    Img_float_org = np.array(Image.open(mxif_cell_mask_fn))

    Img_float = np.zeros(Img_fix.shape) + 255  # create an array to save MxIF image
    start = int((Img_fix.shape[0] - Img_float_org.shape[0]) / 2)
    Img_float[start: start + Img_float_org.shape[0], start:start + Img_float_org.shape[1]] = Img_float_org

    plt.imshow(Img_fix[:,:,0], cmap="gray")
    plt.show()
    plt.imshow(Img_float[:,:,0], cmap="gray")
    plt.show()

    con_s = dict(angle=[0, 2], scale=[1, 4])
    sim = ird.similarity(Img_fix[:,:,0], Img_float[:,:,0], numiter=3, order=3, filter_pcorr=0, constraints=con_s)

    new_img = ird.transform_img_dict(Img_float, sim)
    plt.imshow(new_img, cmap='gray')
    plt.show()

    # plt.imshow(sim["timg"], cmap='gray')
    # plt.show()

    tvec = sim["tvec"].round(4)
    angle = sim["angle"]
    score = sim["success"]
    offset = [tvec[1], tvec[0]]

    print(offset)







