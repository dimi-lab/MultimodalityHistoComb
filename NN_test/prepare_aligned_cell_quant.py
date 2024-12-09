import math
from sklearn.metrics import pairwise_distances
from math import cos, sin, asin, acos
import numpy as np
import os
import copy
from natsort import natsorted
from glob import glob
import matplotlib.pyplot as plt
import random
import pandas as pd
from probreg import cpd, gmmtree, filterreg, bcpd, math_utils, l2dist_regs
from probreg import callbacks
from visualize_cell import *
import networkx as nx
import pygmtools as pygm
from matplotlib.patches import ConnectionPatch
from shapely.geometry.polygon import Polygon
from shapely.geometry import box
from shapely.geometry import Point
import cv2
import tifffile as tf
from scipy.spatial import KDTree

FOV_img_dir = "/temp/Ovarian_TMA"
data_root_dir = "/temp/Ovarian_TMA/AlignmentEval"
output_quant_dir = "/temp/Ovarian_TMA/AlignmentEval/ApplyAlignment/AlignedCellsQuant"
output_data_dir = "/temp/Ovarian_TMA/AlignmentEval/ApplyAlignment/AlignedCells"
DEBUG_Plot = 0  #
SAVE_TRANSFORMED_HE = 0

def get_TMA_core_list(img_path: str) -> list:
    img_fn_list = os.listdir(img_path)
    roi_list = [i.split("_")[0] for i in img_fn_list]
    return list(set(roi_list))

def get_cell_loc(HE_quant_fn):
    HE_quant_df = pd.read_csv(HE_quant_fn, sep='\t')
    he_x = HE_quant_df["Centroid X µm"]
    he_y = HE_quant_df["Centroid Y µm"]
    source = np.array([he_x, he_y]).T
    return source

def get_selected_cell_quant(quant_fn, selected_idx):
    df = pd.read_csv(quant_fn, sep='\t')
    select_df = df.loc[selected_idx]
    return select_df


def apply_aff_trans2points(source, M):
    '''
    source: cell centroids from H&E segmentation (N*2), N is the number of cell
    '''
    new_XY = np.vstack((np.transpose(source), np.ones((1, len(source)))))
    transformed_loc = np.dot(M, new_XY)
    affined_points = np.transpose(transformed_loc[0:2, :])
    return affined_points

def get_cells_in_box(cell_locations, xmin, xmax, ymin,ymax):
    if not isinstance(cell_locations, np.ndarray):
        cell_locations = np.array(cell_locations)
    x_idx = (cell_locations[:, 0] < xmax) & (cell_locations[:, 0] >= xmin)
    y_idx = (cell_locations[:, 1] < ymax) & (cell_locations[:, 1] >= ymin)
    selected_cells = cell_locations[x_idx & y_idx]
    idx_list = np.array(np.where(x_idx & y_idx)).flatten() # TODO:
    return selected_cells, idx_list

def save_selected_points_to_csv(points, idx_list, output_fn):
    fp = open(output_fn, "w")
    wrt_str = "x_pix,y_pix\n"
    for p in idx_list:
        wrt_str += str(points[p][0]) + "," + str(points[p][1]) + "\n"
    fp.write(wrt_str)
    fp.close()

def save_transformed_HE(he_img, M, target_shape, target_pixel_size, affine_HE_fn):
    affined_image = cv2.warpAffine(src=he_img, M=M, dsize=target_shape)
    res_img = affined_image.astype(np.uint8)
    resolution = (10000 / target_pixel_size, 10000 / target_pixel_size)  # convert um to cm for saving
    tf.imwrite(affine_HE_fn, res_img, photometric='rgb', resolution=resolution,
               resolutionunit="CENTIMETER")

def calculate_KDTree_loss_K(qr_points, target, k=5):
    target_tree = KDTree(target, leafsize=10)
    rmse = sum(target_tree.query(qr_points, k)[0].flatten()) / qr_points.shape[0]
    return rmse

# euclidean_distance
def calculate_AED(a, b):
    dist = np.linalg.norm(a-b)/a.shape[0]
    return dist

def calculate_sigma2(X, Y):
    """
       Initialize the variance (sigma2).

       Attributes
       ----------
       X: numpy array
           NxD array of points for target.

       Y: numpy array
           MxD array of points for source.

       Returns
       -------
       sigma2: float
           Initial variance.
       """
    (N, D) = X.shape
    (M, _) = Y.shape
    diff = X[None, :, :] - Y[:, None, :]
    err = abs(diff)
    return np.sum(err) / (D * M * N)

def corr2_coeff(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)

    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))

def align_cell_segmentation(source, target):
    # compute cpd registration
    # set update_scale to False, as the cell locations were denoted with microns (calculated based on pixel size)
    # TODO: But how to apply the affine transformation obtained from this registration if the pixel sizes are different?
    # TODO: sx=sy=source_pixel_size/target_pixel_size
    tf_param, sigma2, q = cpd.registration_cpd(source, target, maxiter=150, update_scale=False)

    result = copy.deepcopy(source)
    r_points = tf_param.transform(result)

    return tf_param, r_points, sigma2, q

#############################################
sec_id_list = ["Sec1", "Sec2"]
seg_method_list = ["StarDist", "Watershed"]
all_counts = []
legend_list = []

HE_pixel_size = 0.2201  # unit micron
MxIF_pixel_size = 0.325

save_str = "sec,roi,gt_degree,cpd_degree,gt_shift,cpd_shift,before_trans_dist,gt_N1_dist,cpd_N1_dist,gt_kdtree1_loss,cpd_kdtree1_loss" \
           ",gt_kdtree2_loss,cpd_kdtree2_loss,gt_kdtree4_loss,cpd_kdtree4_loss," \
           "gt_N1_sigma,cpd_N1_sigma,gt_sigma,cpd_sigma,q\n"

for seg_method in seg_method_list:
    if seg_method == "StarDist":
        mxif_cell_quant_dir = os.path.join(data_root_dir, "QuPathAnnoProj_MxIF", "export")
    elif seg_method == "Watershed":
        mxif_cell_quant_dir = os.path.join(data_root_dir, "QuPathAnnoProj_MxIF_watershed", "export")
    else:
        raise Exception("Segmentation method is not defined")

    for sec_id in sec_id_list:
        print("processing %s" % sec_id)
        if seg_method == "StarDist":
            HE_cell_quant_dir = os.path.join(data_root_dir, "QuPathAnnoProj_HE_" + sec_id, "export")
        elif seg_method == "Watershed":
            HE_cell_quant_dir = os.path.join(data_root_dir, "QuPathAnnoProj_HE_" + sec_id +"_watershed", "export")
        else:
            raise Exception("Segmentation method is not defined")

        # compare the cell number differences between MxIF and H&E
        ROI_list = natsorted(get_TMA_core_list(HE_cell_quant_dir))
        counts = []
        for roi in ROI_list:
            print("processing %s" % roi)

            MxIF_quant_fn = glob(os.path.join(mxif_cell_quant_dir, roi + "_*_QUANT.tsv"))[0]
            HE_quant_fn = glob(os.path.join(HE_cell_quant_dir, roi + "*_QUANT.tsv"))[0]

            # unit: um
            HE_centroids = get_cell_loc(HE_quant_fn)
            MxIF_centroids = get_cell_loc(MxIF_quant_fn)

            ground_truth_affine_dir = os.path.join(data_root_dir, sec_id+"GroundTruth")
            [gt_theta, gt_degrees, gt_s, gt_delta, gt_M] = load_trans(os.path.join(ground_truth_affine_dir, roi + "_trans.dat"))

            print("Ground Truth")
            print("\t Rotate: %f degrees %f" % (gt_theta, gt_degrees))
            print("\t Scale: %f." % gt_s)
            print("\t Translation: (x=%f, y=%f). Shift distance: %f" % (gt_M[0, 2], gt_M[1, 2], gt_delta))

            MxIF_contour_fn = os.path.join(mxif_cell_quant_dir, roi + "_cell_poly.csv")
            HE_contour_fn = os.path.join(HE_cell_quant_dir, roi + "_cell_poly.csv")

            # unit: pixel
            HE_contour_pix = np.array(load_QuPath_cell_poly(HE_contour_fn))
            MxIF_contour_pix = np.array(load_QuPath_cell_poly(MxIF_contour_fn))

            transformed_HE_contour = np.array(apply_aff_trans2poly(HE_contour_pix * HE_pixel_size, gt_M))  # unit: um
            transformed_HE_centroids = np.array(apply_aff_trans2points(HE_centroids, gt_M))                # unit: um
            transformed_HE_centroids_pix = transformed_HE_centroids/MxIF_pixel_size               # unit: pixel
            MxIF_centroids_pix = MxIF_centroids / MxIF_pixel_size                                 # unit: pixel
            HE_centroids_pix = HE_centroids / HE_pixel_size                                       # unit: pixel
            if DEBUG_Plot == 1:
                plot_two_contours(MxIF_contour_pix * MxIF_pixel_size, transformed_HE_contour, roi + "_contour")
                plot_contours_and_centroids(transformed_HE_contour, MxIF_centroids, ["trans_HE_contour", "MxIF_centroid"], roi + "_contour_and_points")
                plot_contours_and_centroids(MxIF_contour_pix * MxIF_pixel_size, transformed_HE_centroids, ["MxIF_contour", "trans_HE_centroid"], roi + "_contour_and_points")
                plot_two_centroids(MxIF_centroids, transformed_HE_centroids, ["MxIF_centroid", "trans_HE_centroid"],roi + "_aligned_points")
                plot_two_centroids(MxIF_centroids_pix, transformed_HE_centroids_pix, ["MxIF_centroid", "trans_HE_centroid"],roi + "_aligned_points")

            N0_cell_idx_MxIF = {}  # MxIF_idx: [],
            N1_cell_idx_MxIF = {}  # MxIF_idx: [HE_idx]
            N2_cell_idx_MxIF = {}  # MxIF_idx: [HE_idx_list]
            for idx, MxIF_cell_poly_points in enumerate(MxIF_contour_pix):
                poly = Polygon(MxIF_cell_poly_points)
                x_min, y_min, x_max, y_max = poly.bounds
                selected_cells, selected_cell_idx_HE = get_cells_in_box(transformed_HE_centroids_pix, x_min, x_max, y_min, y_max)

                if len(selected_cells) == 0:
                    N0_cell_idx_MxIF[idx] = []
                    # print("N=0")
                else:
                    if DEBUG_Plot == 2:
                        plt.figure(dpi=300)
                        contour = MxIF_cell_poly_points
                        plt.fill(contour[:, 1], contour[:, 0], edgecolor='r', linewidth=2, fill=False)
                        plt.scatter(selected_cells[:, 1], selected_cells[:, 0], c='g', s=40)
                        r_patch = mpatches.Patch(color='red', label="MxIF")
                        g_patch = mpatches.Patch(color='green', label="H&E")
                        plt.axis('image')
                        plt.legend(handles=[r_patch, g_patch])
                        plt.title("single cell alignment evaluation")
                        plt.tight_layout()
                        plt.show()

                    within_list = []
                    for sc in selected_cells:
                        p = Point(sc)
                        if p.within(poly):
                            within_list.append(True)
                        else:
                            within_list.append(False)
                    N = sum(within_list)
                    HE_idx_within = np.array(selected_cell_idx_HE)[within_list]
                    if N == 0:
                        N0_cell_idx_MxIF[idx] = []
                        # print("N=0")
                    elif N == 1:
                        N1_cell_idx_MxIF[idx] = np.array(selected_cell_idx_HE)[within_list]
                        # print("N=1")
                    else:
                        N2_cell_idx_MxIF[idx] = np.array(selected_cell_idx_HE)[within_list]
                        # print("N>1")

            debug_str = '''
            ROI: %s
            MxIF cell number: %d; H&E cell number: %d; Delta=%d
            N0 cells: %d
            N1 cells: %d
            N1+ cells: %d
            '''
            print(debug_str % (roi, len(MxIF_centroids), len(HE_centroids), len(MxIF_centroids) - len(HE_centroids),
                               len(N0_cell_idx_MxIF.keys()), len(N1_cell_idx_MxIF.keys()), len(N2_cell_idx_MxIF.keys())))

            selected_MxIF_cell_idx_N0 = N0_cell_idx_MxIF.keys()
            selected_MxIF_contour_N0 = [MxIF_contour_pix[c] for c in selected_MxIF_cell_idx_N0]

            HE_cell_idx_N1 = [i for k in N1_cell_idx_MxIF.keys() for i in N1_cell_idx_MxIF.get(k)]
            selected_MxIF_cell_idx_N1 = N1_cell_idx_MxIF.keys()
            selected_MxIF_contour_N1 = [MxIF_contour_pix[c] for c in selected_MxIF_cell_idx_N1]

            HE_cell_idx_N2 = [i for k in N2_cell_idx_MxIF.keys() for i in N2_cell_idx_MxIF.get(k)]
            selected_MxIF_cell_idx_N2 = N2_cell_idx_MxIF.keys()
            selected_MxIF_contour_N2 = [MxIF_contour_pix[c] for c in selected_MxIF_cell_idx_N2]

            selected_HE_quant_df = get_selected_cell_quant(HE_quant_fn, HE_cell_idx_N1)
            selected_MxIF_quant_df = get_selected_cell_quant(MxIF_quant_fn, selected_MxIF_cell_idx_N1)
            save_fn = os.path.join(output_quant_dir, roi + "_1on1_HE_quant.tsv")
            selected_HE_quant_df.reset_index(inplace=False).to_csv(save_fn, sep='\t')
            save_fn = os.path.join(output_quant_dir, roi + "_1on1_MxIF_quant.tsv")
            selected_MxIF_quant_df.reset_index(inplace=False).to_csv(save_fn, sep='\t')

            if DEBUG_Plot == 3:
                plot_contors(selected_MxIF_contour_N0, roi + "_N0")
                plot_contours_and_centroids(selected_MxIF_contour_N1, transformed_HE_centroids_pix[HE_cell_idx_N1],
                                        ["MxIF_contour", "trans_HE_centroid"], roi + "_N1")
                plot_contours_and_centroids(selected_MxIF_contour_N2, transformed_HE_centroids_pix[HE_cell_idx_N2],
                                            ["MxIF_contour", "trans_HE_centroid"], roi + "_N2+")

            cell_loc_fn = os.path.join(output_data_dir, roi + "_N0_MxIF.csv")
            save_selected_points_to_csv(MxIF_centroids, selected_MxIF_cell_idx_N0, cell_loc_fn)

            cell_loc_fn = os.path.join(output_data_dir, roi + "_N1_MxIF.csv")
            save_selected_points_to_csv(MxIF_centroids, selected_MxIF_cell_idx_N1, cell_loc_fn)
            cell_loc_fn = os.path.join(output_data_dir, roi + "_N1_HE.csv")
            save_selected_points_to_csv(transformed_HE_centroids_pix, HE_cell_idx_N1, cell_loc_fn)

            cell_loc_fn = os.path.join(output_data_dir, roi + "_N2_MxIF.csv")
            save_selected_points_to_csv(MxIF_centroids, selected_MxIF_cell_idx_N2, cell_loc_fn)
            cell_loc_fn = os.path.join(output_data_dir, roi + "_N2_HE.csv")
            save_selected_points_to_csv(transformed_HE_centroids_pix, HE_cell_idx_N2, cell_loc_fn)

            '''
            ###############################################
            # Do it outside of Python. Show the saved points in QuPath. refer to code: load_locations.groovy
            ###############################################
            '''
            '''
            ###############################################
            # calculate loss
            ###############################################
            '''

            # TODO: # calculate loss, save later
            # All_init_loss = calculate_KDTree_loss_K(HE_centroids_pix, MxIF_centroids_pix)
            # # Ground Truth loss
            # gt_loss = calculate_KDTree_loss_K(transformed_HE_centroids_pix, MxIF_centroids_pix)
            # gt_N1_loss = calculate_KDTree_loss_K(transformed_HE_centroids_pix[HE_cell_idx_N1], MxIF_centroids_pix[list(selected_MxIF_cell_idx_N1)])
            # gt_N1_loss_total = calculate_KDTree_loss_K(transformed_HE_centroids_pix[HE_cell_idx_N1], MxIF_centroids_pix)
            cpd_affine_dir = os.path.join(data_root_dir, sec_id + "_stardist_CPD")
            [cpd_theta, cpd_degrees, cpd_s, cpd_delta, cpd_M] = load_trans(
                os.path.join(cpd_affine_dir, roi + "_trans.dat"))
            print("CPD")
            print("\t Rotate: %f degrees %f" % (cpd_theta, cpd_degrees))
            print("\t Scale: %f." % cpd_s)
            print("\t Translation: (x=%f, y=%f). Shift distance: %f" % (cpd_M[0, 2], cpd_M[1, 2], cpd_delta))

            tf_param, r_points, sigma2, q = align_cell_segmentation(HE_centroids, MxIF_centroids)

            cpd_transformed_HE_centroids = apply_aff_trans2points(HE_centroids, cpd_M)

            if DEBUG_Plot == 3:
                plot_two_centroids(MxIF_centroids_pix[list(selected_MxIF_cell_idx_N1)], transformed_HE_centroids_pix[HE_cell_idx_N1], ["MxIF_centroid", "trans_HE_centroid"],
                                   roi + " Ground Truth")
                plot_two_centroids(MxIF_centroids[list(selected_MxIF_cell_idx_N1)],
                                   cpd_transformed_HE_centroids[HE_cell_idx_N1], ["MxIF_centroid", "trans_HE_centroid"],
                                   roi + " CPD")




            dist_before_trans = calculate_sigma2(HE_centroids_pix, MxIF_centroids_pix)

            # absolute differences between one-on-one matched cells
            abs_gt_dist_loss = calculate_AED(transformed_HE_centroids_pix[HE_cell_idx_N1],
                                             MxIF_centroids_pix[list(selected_MxIF_cell_idx_N1)])
            abs_cpd_dist_loss = calculate_AED(cpd_transformed_HE_centroids[HE_cell_idx_N1],
                                              MxIF_centroids_pix[list(selected_MxIF_cell_idx_N1)])

            # absolute differences between a cell to the closest cell
            abs_gt_dist_loss_all = calculate_KDTree_loss_K(transformed_HE_centroids_pix, MxIF_centroids_pix, k=1)
            abs_cpd_dist_loss_all = calculate_KDTree_loss_K(cpd_transformed_HE_centroids, MxIF_centroids_pix, k=1)

            # absolute differences between a cell to the closest 4 cell
            abs_gt_dist_loss_all_2 = calculate_KDTree_loss_K(transformed_HE_centroids_pix, MxIF_centroids_pix, k=2) / 2
            abs_cpd_dist_loss_all_2 = calculate_KDTree_loss_K(cpd_transformed_HE_centroids, MxIF_centroids_pix,
                                                              k=2) / 2

            # absolute differences between a cell to the closest 4 cell
            abs_gt_dist_loss_all_4 = calculate_KDTree_loss_K(transformed_HE_centroids_pix, MxIF_centroids_pix, k=4) / 4
            abs_cpd_dist_loss_all_4 = calculate_KDTree_loss_K(cpd_transformed_HE_centroids, MxIF_centroids_pix,
                                                              k=4) / 4

            # absolute differences between a cell to the closest 4 cell
            abs_gt_dist_loss_all_8 = calculate_KDTree_loss_K(transformed_HE_centroids_pix, MxIF_centroids_pix, k=8) / 8
            abs_cpd_dist_loss_all_8 = calculate_KDTree_loss_K(cpd_transformed_HE_centroids, MxIF_centroids_pix,
                                                              k=8) / 8

            gt_N1_sigma = calculate_sigma2(transformed_HE_centroids_pix[HE_cell_idx_N1], MxIF_centroids_pix[list(selected_MxIF_cell_idx_N1)])
            cpd_N1_sigma = calculate_sigma2(cpd_transformed_HE_centroids[HE_cell_idx_N1], MxIF_centroids_pix[list(selected_MxIF_cell_idx_N1)])
            gt_sigma = calculate_sigma2(transformed_HE_centroids_pix, MxIF_centroids_pix)
            cpd_sigma = calculate_sigma2(cpd_transformed_HE_centroids, MxIF_centroids_pix)

            save_str += ",".join([sec_id, roi, str(gt_degrees), str(cpd_degrees), str(gt_delta), str(cpd_delta),
                                  str(dist_before_trans), str(abs_gt_dist_loss), str(abs_cpd_dist_loss),
                                           str(abs_gt_dist_loss_all), str(abs_cpd_dist_loss_all),
                                           str(abs_gt_dist_loss_all_2), str(abs_cpd_dist_loss_all_2),
                                           str(abs_gt_dist_loss_all_4), str(abs_cpd_dist_loss_all_4),
                                           str(gt_N1_sigma), str(cpd_N1_sigma),
                                           str(gt_sigma), str(cpd_sigma), str(q)]) + "\n"

            # out1 = corr2_coeff(transformed_HE_centroids_pix[HE_cell_idx_N1],
            #                                   MxIF_centroids_pix[list(selected_MxIF_cell_idx_N1)])
            # out2 = corr2_coeff(cpd_transformed_HE_centroids_pix[HE_cell_idx_N1],
            #                   MxIF_centroids_pix[list(selected_MxIF_cell_idx_N1)])

            # print("debug")
            # All_loss_init = calculate_KDTree_loss_K(HE_centroids_pix, MxIF_centroids_pix)
            # All_loss = calculate_KDTree_loss_K(transformed_HE_centroids_pix, MxIF_centroids_pix)

            '''
            sigma2: refer: https://github.com/siavashk/pycpd/blob/e5ca02d2501fb4b633c1664de939c336ffa2349e/pycpd/emregistration.py#L205
            (N, D) = X.shape
            (M, _) = Y.shape
            diff = X[None, :, :] - Y[:, None, :]
            err = diff ** 2
            return np.sum(err) / (D * M * N)

            q: float
            The objective function value that represents the misalignment between source
            and target point clouds.
            '''

            # N=0

            # N=1

            # N>1



            '''
            ###############################################
            # Save transformed H&E images
            ###############################################
            '''
            if SAVE_TRANSFORMED_HE:
                if sec_id == "Sec1":
                    HE_img_data_dir = os.path.join(FOV_img_dir, "HE_FOVs", "same_section")
                    Aff_HE_img_data_dir = os.path.join(FOV_img_dir, "AlignmentEval", "ApplyAlignment", "Sec1Output")
                else:
                    HE_img_data_dir = os.path.join(FOV_img_dir, "HE_FOVs", "Rab_Spine-22R919-A-SERBSVG-5X-08")
                    Aff_HE_img_data_dir = os.path.join(FOV_img_dir, "AlignmentEval", "ApplyAlignment", "Sec2Output")
                MxIF_img_data_dir = os.path.join(FOV_img_dir, "MxIF_FOVs", "Slide2050_24Plex")
                he_fn = os.path.join(HE_img_data_dir, roi + ".tif")
                he_img = tf.TiffFile(he_fn).pages[0].asarray().astype(np.float)
                mxif_fn = os.path.join(MxIF_img_data_dir, roi + ".tif")
                mxif_img = tf.TiffFile(mxif_fn)
                dapi_img = mxif_img.pages[0].asarray().astype(np.float)
                affine_HE_fn = os.path.join(Aff_HE_img_data_dir, roi + "_HE_gt_aff.tif")
                save_transformed_HE(he_img, gt_M, dapi_img.shape, MxIF_pixel_size, affine_HE_fn)

        break
    break
fp = open(os.path.join(data_root_dir, "eval_metrics.csv"), 'w')
fp.write(save_str)
fp.close()

print("Debug")














