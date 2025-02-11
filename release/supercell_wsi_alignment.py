import os
import argparse
import random

from probreg import cpd
# from utils import *
import copy
# from visual_utils import *
from random import sample
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.cluster import AffinityPropagation
from itertools import cycle, islice
import pandas as pd
import numpy as np

def get_cell_loc(HE_quant_fn, sample_rate, col_name_x='Centroid X µm', col_name_y='Centroid Y µm'):
    '''
    get
    :param HE_quant_fn: cell quantification file
    :param col_name_x:  column name save the x coordinate
    :param col_name_y:  column name save the y coordinate
    :return: cell locations
    '''
    HE_quant_df = pd.read_csv(HE_quant_fn, sep='\t')
    he_x = HE_quant_df[col_name_x]
    he_y = HE_quant_df[col_name_y]
    tmp = np.array([he_x, he_y]).T
    cell_centroids = np.array(sample(list(tmp), int(len(he_x) * sample_rate)), dtype=np.float32)
    return cell_centroids


def get_cell_samples(quant_fn, spacing=100, points_per_cell=10, col_name_x='Centroid X µm', col_name_y='Centroid Y µm'):
    '''
    get
    :param quant_fn: cell quantification file
    :param col_name_x:  column name save the x coordinate
    :param col_name_y:  column name save the y coordinate
    :return: cell locations
    '''
    quant_df = pd.read_csv(quant_fn, sep='\t')
    loc_x = quant_df[col_name_x]
    loc_y = quant_df[col_name_y]

    x_min, x_max = min(loc_x), max(loc_x)  # x-coordinate range
    y_min, y_max = min(loc_y), max(loc_y)  # y-coordinate range

    # Create the edges of the grid cells
    x_edges = np.arange(x_min, x_max + spacing, spacing)
    y_edges = np.arange(y_min, y_max + spacing, spacing)

    sampled_points = []
    # Loop over each grid cell to sample points
    for x in x_edges[:-1]:
        for y in y_edges[:-1]:
            # Sample random points within the current grid cell
            idx_x = np.logical_and(loc_x < x + spacing, loc_x >= x)
            idx_y = np.logical_and(loc_y >= y, loc_y < y + spacing)
            idx_xy = np.logical_and(idx_x, idx_y)
            if sum(idx_xy) >= points_per_cell:

                selected_samples = random.sample(list(np.column_stack([loc_x[idx_xy], loc_y[idx_xy]])), points_per_cell)
                # # y_samples = random.sample(list(loc_y[idx_y]), points_per_cell)
                #
                # # Combine x and y samples into (x, y) pairs
                # cell_samples = np.column_stack((x_samples, y_samples))
                sampled_points += selected_samples

    sampled_points = np.array(sampled_points)     # Concatenate all sampled points into a single array
    return sampled_points

def update_cell_loc(HE_quant_fn, updated_loc_xy, output_fn, col_name_x='Centroid X µm', col_name_y='Centroid Y µm'):
    '''
    update the transformed cell locations and save to tsv file.
    :param HE_quant_fn: original cell quantification file
    :param updated_loc_xy: transformed/updated cell locations
    :param output_fn: file name to save the updated cell quantification
    :param col_name_x:  column name save the x coordinate
    :param col_name_y:  column name save the y coordinate
    :return:
    '''
    HE_quant_df = pd.read_csv(HE_quant_fn, sep='\t')
    HE_quant_df[col_name_x] = updated_loc_xy[:, 0]
    HE_quant_df[col_name_y] = updated_loc_xy[:, 1]
    HE_quant_df.to_csv(output_fn, sep='\t')


def sample_and_plot_centroids_trans(s_points, t_points, t_s_points, legend, title, out_dir, fn, sample_rate=0.3, size=1):
    '''
    create scatter plot to show the points before and after transformation
    :param s_points: source points
    :param t_points: target points
    :param t_s_points: Transformed points
    :param legend: legend of the plot
    :param title: title of the plot
    :param out_dir: output directory to save the plot
    :param fn: output file name
    :param sample_rate: sample the points to show if there are too many. Useful for whole slide image, because cell number is large
    :param size: dot size in the scatter plot
    :return:
    '''
    plt.figure(dpi=300)

    if sample_rate!=1:
        s_points = np.array(sample(list(s_points), int(len(s_points)*sample_rate)))
        t_points = np.array(sample(list(t_points), int(len(t_points)*sample_rate)))
        t_s_points = np.array(sample(list(t_s_points), int(len(t_s_points) * sample_rate)))

    plt.scatter(s_points[:, 1], s_points[:, 0], c='r', s=size)
    plt.scatter(t_points[:, 1], t_points[:, 0], c='g', s=size)
    plt.scatter(t_s_points[:, 1], t_s_points[:, 0], c='b', s=size)
    r_patch = mpatches.Patch(color='red', label=legend[0])
    g_patch = mpatches.Patch(color='green', label=legend[1])
    b_patch = mpatches.Patch(color='blue', label=legend[2])
    plt.axis('image')
    plt.legend(handles=[r_patch, g_patch, b_patch])
    plt.title(title)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, fn))

def save_clustering_plot(cell_locations, y_pred, title, results_dir, fn):
    '''
    create scatter plot to show the clustering results
    :param cell_locations:  cell locations
    :param y_pred: class label of the clustering results
    :param results_dir: directory to save the plot
    :param fn:  file name to save the plot
    :return:
    '''
    colors = np.array(
        list(
            islice(
                cycle(
                    [
                        "#377eb8",
                        "#ff7f00",
                        "#4daf4a",
                        "#f781bf",
                        "#a65628",
                        "#984ea3",
                        "#999999",
                        "#e41a1c",
                        "#dede00",
                    ]
                ),
                int(max(y_pred) + 1),
            )
        )
    )
    colors = np.append(colors, ["#000000"])
    plt.figure(1, dpi=300)
    plt.scatter(cell_locations[:, 1], cell_locations[:, 0], s=2, color=colors[y_pred])
    plt.axis('image')
    plt.title(title + " Affinity Propagation Clustering")
    plt.tight_layout()
    # plt.show()
    sv_fn = os.path.join(results_dir, fn + ".png")
    plt.savefig(sv_fn)
    plt.close()

def save_cell_sampling_plot(cell_locations, title, results_dir, fn):
    plt.figure(1, dpi=300)
    plt.scatter(cell_locations[:, 1], cell_locations[:, 0], s=1)
    plt.axis('image')
    plt.title(title + " Affinity Propagation Clustering")
    plt.tight_layout()
    # plt.show()
    sv_fn = os.path.join(results_dir, fn + ".png")
    plt.savefig(sv_fn)
    plt.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--he_quant_fn", # cell quantification from source image. We usually treat H&E image as
                        # source image, because there is less image channels, easier to write the transformed image.
                        required=True,
                        dest="he_quant_fn",
                        help="source, H&E Cell quantification from QuPath. Cell coordinate need to be saved in µm")
    parser.add_argument("-t", "--mxif_quant_fn",  # cell quantification from target image.
                        required=True,
                        dest='mxif_quant_fn',
                        help="target, Cell quantification from QuPath. Cell coordinate need to be saved in µm")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Save debug information"
                        )
    parser.add_argument("-o", "--output_dir",
                        default=os.getcwd(),
                        dest='output_dir',
                        help="Directory for transformed quant file, and other output ")

    # args = parser.parse_args()
    # MxIF_quant_fn = args.mxif_quant_fn
    # HE_quant_fn = args.he_quant_fn
    # output_dir = args.output_dir
    # verbose = args.verbose
    #
    # MxIF_quant_fn = "/projects/wangc/raymond/markovic.bms.gastric/SLIDE-3026_sec01.ome.tif_QUANT.tsv"
    # HE_quant_fn = "/projects/wangc/raymond/markovic.bms.gastric/SLIDE-3026_sec1_HnE_QUANT.tsv"

    MxIF_quant_fn = "/projects/wangc/raymond/markovic.bms.gastric/SLIDE-3026_sec02.ome.tif_QUANT.tsv"
    HE_quant_fn = "/projects/wangc/raymond/markovic.bms.gastric/SLIDE-3026_sec2_HnE_QUANT.tsv"
    #
    # MxIF_quant_fn = "/projects/wangc/raymond/markovic.bms.gastric/SLIDE-3026_sec03.ome.tif_QUANT.tsv"
    # HE_quant_fn = "/projects/wangc/raymond/markovic.bms.gastric/SLIDE-3026_sec3_HnE_QUANT.tsv"
    output_dir = "/projects/wangc/raymond/markovic.bms.gastric"
    # MxIF_quant_fn = "/projects/wangc/jun/data/Raymond_WSI_alignment/data/MxIF_I-17.tif_QUANT.tsv"
    # HE_quant_fn = "/projects/wangc/jun/data/Raymond_WSI_alignment/data/HE_I-17.tif_QUANT.tsv"
    # output_dir = "/projects/wangc/jun/data/Raymond_WSI_alignment/output"
    verbose = True

    HE_tif_fn = os.path.split(HE_quant_fn)[1].replace(".tif_QUANT.tsv", "")
    MxIF_tif_fn = os.path.split(MxIF_quant_fn)[1].replace(".tif_QUANT.tsv", "")

    All_HE_cell_locations = get_cell_loc(HE_quant_fn, sample_rate=1)
    All_MxIF_cell_locations = get_cell_loc(MxIF_quant_fn, sample_rate=1)
    HE_cell_locations = get_cell_samples(HE_quant_fn, spacing=200, points_per_cell=1)   # you may need to specify column names to get the cell coordinates
    MxIF_cell_locations = get_cell_samples(MxIF_quant_fn, spacing=200, points_per_cell=1)
    # print("H&E cell clustering")
    # HE_clustering = AffinityPropagation(damping=0.5, random_state=5).fit(HE_cell_locations)
    # HE_y_pred = HE_clustering.labels_.astype(int)
    # print("MxIF cell clustering")
    # MxIF_clustering = AffinityPropagation(damping=0.5, random_state=5).fit(MxIF_cell_locations)
    # MxIF_y_pred = MxIF_clustering.labels_.astype(int)
    #
    # print("Saving clustering results into plots")
    # save_clustering_plot(HE_cell_locations, HE_y_pred, "HE", output_dir, HE_tif_fn)
    # save_clustering_plot(MxIF_cell_locations, MxIF_y_pred, "MxIF", output_dir, MxIF_tif_fn)
    #

    save_cell_sampling_plot(HE_cell_locations, "HE", output_dir, HE_tif_fn)
    save_cell_sampling_plot(MxIF_cell_locations, "MxIF", output_dir, MxIF_tif_fn)

    # MxIF_cluster_centroids = np.array([MxIF_cell_locations[MxIF_y_pred == i].mean(axis=0) for i in set(MxIF_y_pred)])
    # HE_cluster_centroids = np.array([HE_cell_locations[HE_y_pred == i].mean(axis=0) for i in set(HE_y_pred)])

    # tf_param, sigma2, q = cpd.registration_cpd(HE_cell_locations.T, MxIF_cell_locations.T, maxiter=40, tol=1.0e-11,
    #                                            update_scale=False)
    tf_param, sigma2, q = cpd.registration_cpd(HE_cell_locations.T, MxIF_cell_locations.T, maxiter=40, tol=1.0e-11,
                                               update_scale=True)
    # tf_param, sigma2, q = cpd.registration_cpd(HE_cluster_centroids.T, MxIF_cluster_centroids.T, maxiter=20, tol=1.0e-11,
    #                                            update_scale=False)
    print("Apply the transformation to H&E cell coordinates.")
    # r_points = tf_param.transform(copy.deepcopy(All_HE_cell_locations.T))
    r_points = tf_param.transform(copy.deepcopy(HE_cell_locations.T))

    print("Update the transformed H&E cell coordinates to the table.")
    # output_fn = os.path.join(output_dir, "transformed_" + os.path.split(HE_quant_fn)[1])
    # update_cell_loc(HE_quant_fn, r_points.T, output_fn)
    if verbose:
        title = "Cell Centroids Before and after alignment"
        legend = ["HE cells", "MxIF cells", "trans_HE"]
        fn = "log_" + HE_tif_fn + "_centroids_alignment.png"
        # sample_and_plot_centroids_trans(All_HE_cell_locations, All_MxIF_cell_locations, r_points.T, legend, title, output_dir, fn)
        sample_and_plot_centroids_trans(HE_cell_locations, MxIF_cell_locations, r_points.T, legend, title,
                                        output_dir, fn)

    ## TODO: export affine transformation, maybe need to convert um to pixel coordinates


    # TODO: add refinement after CPD  (assigned to Jun)

    # TODO: Currently, unable to write whole slide image, so skip the step of saving the transformed image.
    # HE_pixel_size = 0.2201  # unit micron
    # MxIF_pixel_size = 0.325
    # theta, degrees, s, delta, cpd_M = get_M_from_cpd(tf_param, HE_pixel_size, MxIF_pixel_size)



    # affine_HE_fn = os.path.join(output_dir, tif_fn)
    # he_img = tf.TiffFile(HE_img_fn).pages[0].asarray().astype(np.float)  # read the images
    # mxif_img = tf.TiffFile(MxIF_img_fn)
    # dapi_img = mxif_img.pages[0].asarray().astype(np.float)
    # trans_he_img = save_transformed_HE(he_img, cpd_M, dapi_img.shape, MxIF_pixel_size, affine_HE_fn)
    #
    # if verbose:
    #     s_fn = os.path.join(output_dir, "side_by_side_"+tif_fn)
    #     mxif_img_3_channel1 = mxif_img.pages[0].asarray().astype(np.uint8)
    #     mxif_img_3_channel2 = mxif_img.pages[1].asarray().astype(np.uint8)
    #     mxif_img_3_channel3 = mxif_img.pages[2].asarray().astype(np.uint8)
    #     mxif_img_3_channel = np.stack([mxif_img_3_channel1, mxif_img_3_channel2, mxif_img_3_channel3], axis=2)
    #     mxif_img_3_channel = trans_img_range(mxif_img_3_channel)
    #     save_transformed_HE_and_MxIF(trans_he_img, mxif_img_3_channel, dapi_img.shape, MxIF_pixel_size, s_fn)


























