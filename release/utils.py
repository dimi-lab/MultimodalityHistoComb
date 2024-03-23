import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import tifffile as tf
import cv2
from skimage.measure import label, regionprops

def get_cell_loc(HE_quant_fn, col_name_x='Centroid X µm', col_name_y='Centroid Y µm'):
    HE_quant_df = pd.read_csv(HE_quant_fn, sep='\t')
    he_x = HE_quant_df[col_name_x]
    he_y = HE_quant_df[col_name_y]
    source = np.array([he_x, he_y]).T
    return source

def get_cell_loc_from_StarDist_pred(pred_lbl):
    cell_centroids = []
    for r in regionprops(pred_lbl):
        cell_centroids.append(list(r.centroid))
    return cell_centroids

def apply_aff_trans2points(source, M):
    '''
    source: cell centroids from H&E segmentation (N*2), N is the number of cell
    '''
    new_XY = np.vstack((np.transpose(source), np.ones((1, len(source)))))
    transformed_loc = np.dot(M, new_XY)
    affined_points = np.transpose(transformed_loc[0:2, :])
    return affined_points

def plot_centroids_trans(s_points, t_points, t_s_points, legend, title, out_dir, fn):
    plt.figure(dpi=300)
    plt.scatter(s_points[:, 1], s_points[:, 0], c='r', s=1)
    plt.scatter(t_points[:, 1], t_points[:, 0],  c='g', s=1)
    plt.scatter(t_s_points[:, 1], t_s_points[:, 0], c='b', s=1)
    r_patch = mpatches.Patch(color='red', label=legend[0])
    g_patch = mpatches.Patch(color='green', label=legend[1])
    b_patch = mpatches.Patch(color='blue', label=legend[2])
    plt.axis('image')
    plt.legend(handles=[r_patch, g_patch, b_patch])
    plt.title(title)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, fn))

def plot_cell_detection(img, cell_centroids, out_dir, fn):
    plt.figure(dpi=300)
    plt.imshow(img)
    plt.scatter(cell_centroids[:, 1], cell_centroids[:, 0], c='r', s=2)
    plt.axis('image')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, fn))

def get_M_from_cpd(tf_param, source_pix_size, target_pix_size):
    R = tf_param.rot
    T = tf_param.t
    real_S = source_pix_size/target_pix_size
    M = np.array([[real_S*R[0, 0], R[0, 1], T[1]/target_pix_size],
                  [R[1, 0], real_S*R[1, 1], T[0]/target_pix_size]]).astype(float)

    theta = math.atan(R[1, 0] / R[0, 0])
    degrees = theta * (180.0 / math.pi)
    s = tf_param.scale
    delta = math.sqrt(M[0, 2] ** 2 + M[1, 2] ** 2)

    return [theta, degrees, s, delta, M]

def save_transformed_HE(he_img, M, target_shape, target_pixel_size, affine_HE_fn):
    affined_image = cv2.warpAffine(src=he_img, M=M, dsize=target_shape)
    res_img = affined_image.astype(np.uint8)
    res_img[np.all(res_img == 0, axis=-1)] = [240, 242, 240]  # replace black background with gray
    resolution = (10000 / target_pixel_size, 10000 / target_pixel_size)  # convert um to cm for saving
    tf.imwrite(affine_HE_fn, res_img, photometric='rgb', resolution=resolution,
               resolutionunit="CENTIMETER")
    return res_img

def save_transformed_HE_and_MxIF(trans_he_img, mxif_img_3_channel, target_shape, target_pixel_size, s_fn):
    res_img = np.zeros((target_shape[0], target_shape[1]*2, 3), dtype=np.uint8)
    res_img[:,0:target_shape[1],:] = mxif_img_3_channel
    res_img[:, target_shape[1]:, :] = trans_he_img
    resolution = (10000 / target_pixel_size, 10000 / target_pixel_size)  # convert um to cm for saving
    tf.imwrite(s_fn, res_img, photometric='rgb', resolution=resolution,
               resolutionunit="CENTIMETER")
    return res_img
