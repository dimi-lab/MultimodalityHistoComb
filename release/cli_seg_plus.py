import numpy as np
import os
import argparse
from probreg import cpd
from utils import *
import tifffile as tf
import copy
from stardist.models import StarDist2D, Config2D
from csbdeep.utils import Path, normalize
from gm_utils import *
import matplotlib.pyplot as plt
from LPM import *
import functools
# import tensorflow
# gpus = tensorflow.config.list_physical_devices('GPU')
# tensorflow.config.experimental.set_memory_growth(gpus[0], True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--he_img_fn",
                        required=True,
                        dest="he_img_fn",
                        help="source image, H&E image file, supposed to be ome tiff")
    parser.add_argument("-t", "--mxif_img_fn",
                        required=True,
                        dest="mxif_img_fn",
                        help="target image (MxIF image file), supposed to be ome tiff")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Save debug information"
                        )
    parser.add_argument("-o", "--output_dir",
                        default=os.getcwd(),
                        dest='output_dir',
                        help="Transformed H&E output directory")

    args = parser.parse_args()

    HE_img_fn = args.he_img_fn
    MxIF_img_fn = args.mxif_img_fn
    output_dir = args.output_dir
    verbose = args.verbose

    # HE_img_fn = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/HE_FOVs/same_section/A-8.tif"
    # MxIF_img_fn = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/MxIF_FOVs/Slide2050_24Plex/A-8.tif"
    # # HE_img_fn = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/HE_FOVs/same_section/B-9.tif"
    # # MxIF_img_fn = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/MxIF_FOVs/Slide2050_24Plex/B-9.tif"
    # output_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/test_align_wsi"
    # verbose = True

    tif_fn = os.path.split(HE_img_fn)[1]
    print("Processing %s" % tif_fn)
    Aff_M_fn = os.path.join(output_dir, tif_fn + "_final_M.npy")
    if not os.path.exists(Aff_M_fn):
        he_cell_cent_fn = os.path.join(output_dir, tif_fn + "_he_cell_centroids.npy")
        mxif_cell_cent_fn = os.path.join(output_dir, tif_fn + "_mxif_cell_centroids.npy")

        he_cell_feature_fn = os.path.join(output_dir, tif_fn + "_he_cell_features.npy")
        mxif_cell_feature_fn = os.path.join(output_dir, tif_fn + "_mxif_cell_features.npy")

        # TODO: read pixel size from the image tags
        HE_pixel_size = 0.2201  # unit micron
        MxIF_pixel_size = 0.325

        #################################### Do cell segmentation #################################################
        if not os.path.exists(he_cell_cent_fn) or (not os.path.exists(mxif_cell_cent_fn)):
            print("Getting models from the list")
        StarDist2D.from_pretrained()

        # read the images
        he_img = tf.TiffFile(HE_img_fn).pages[0].asarray().astype(float)
        norm_he_img = normalize(he_img)
        if os.path.exists(he_cell_cent_fn) and os.path.exists(he_cell_feature_fn):
            HE_centroids_pix = np.load(he_cell_cent_fn)
            HE_morph_features = np.load(he_cell_feature_fn)
        else:  # if cell location and feature files don't exist, do cell segmentation
            he_model = StarDist2D.from_pretrained('2D_versatile_he')
            he_lbl, _ = he_model.predict_instances_big(norm_he_img, prob_thresh=0.5, axes='YXC', block_size=2048,
                                                       min_overlap=128, context=128)
            HE_centroids_pix, HE_morph_features = get_cell_info_from_StarDist_pred(he_lbl)
            np.save(he_cell_cent_fn, HE_centroids_pix)
            np.save(he_cell_feature_fn, HE_morph_features)

        mxif_img = tf.TiffFile(MxIF_img_fn)

        mxif_img_3_channel1 = mxif_img.pages[0].asarray().astype(np.uint8)  # DAPI channel
        mxif_img_3_channel2 = mxif_img.pages[1].asarray().astype(np.uint8)
        mxif_img_3_channel3 = mxif_img.pages[2].asarray().astype(np.uint8)

        mxif_img_3_channel = np.stack([mxif_img_3_channel1, mxif_img_3_channel2, mxif_img_3_channel3], axis=2)
        if os.path.exists(mxif_cell_cent_fn) and os.path.exists(mxif_cell_feature_fn):
            MxIF_centroids_pix = np.load(mxif_cell_cent_fn)
            MxIF_morph_features = np.load(mxif_cell_feature_fn)
        else: # if cell location and feature files don't exist, do cell segmentation
            mxif_model = StarDist2D.from_pretrained('2D_versatile_fluo')
            norm_mxif_img = normalize(mxif_img_3_channel1.astype(float), 1, 99.8, axis=(0, 1))
            mxif_lbl, _ = mxif_model.predict_instances_big(norm_mxif_img, axes='YX', block_size=2048, min_overlap=128,
                                                           context=128)
            MxIF_centroids_pix, MxIF_morph_features = get_cell_info_from_StarDist_pred(mxif_lbl)
            np.save(mxif_cell_cent_fn, MxIF_centroids_pix)
            np.save(mxif_cell_feature_fn, MxIF_morph_features)

        if verbose:
            fig, (a, b) = plt.subplots(1, 2, figsize=(16, 8))
            a.scatter(HE_centroids_pix[:, 1], HE_centroids_pix[:, 0], c='r', s=2)
            b.scatter(MxIF_centroids_pix[:, 1], MxIF_centroids_pix[:, 0], c='r', s=2)
            # a.axis('off')
            # b.axis('off')
            plt.show()
            print("H&E Cell Count: %d" % len(HE_centroids_pix))
            print("MxIF Cell Count: %d" % len(MxIF_centroids_pix))

        ##################################### Align based on cell segmentation ########################################
        print("Aligning based on cell segmentation.")
        # if the pixel size is Null, use 1. That means directly use pixel coordinate to align.
        update_scale = False
        if HE_pixel_size is None or MxIF_pixel_size is None:
            HE_pixel_size = 1
            MxIF_pixel_size = 1
            update_scale = True

        HE_centroids_um = np.array(HE_centroids_pix) * HE_pixel_size
        MxIF_centroids_um = np.array(MxIF_centroids_pix) * MxIF_pixel_size
        tf_param, sigma2, q = cpd.registration_cpd(HE_centroids_um, MxIF_centroids_um, maxiter=20, tol=1.0e-11,
                                                   update_scale=update_scale)
        trans_HE_centroids_um = tf_param.transform(copy.deepcopy(HE_centroids_um))
        if verbose:
            title = "Cell Centroids Before and after alignment"
            legend = ["HE cells", "MxIF cells", "trans_HE"]
            fn = "log_" + tif_fn + "_centroids_alignment.png"
            plot_centroids_trans(HE_centroids_um, MxIF_centroids_um, trans_HE_centroids_um, legend, title, output_dir, fn)

        if HE_pixel_size is None or MxIF_pixel_size is None:
            theta, degrees, s, delta, CPD_M = get_M_from_cpd_pix(tf_param)
        else:
            theta, degrees, s, delta, CPD_M = get_M_from_cpd(tf_param, HE_pixel_size, MxIF_pixel_size)

        #############################
        if verbose:
            fig, (a, b) = plt.subplots(1, 2, figsize=(16, 8))
            a.scatter(HE_centroids_um[:, 1], HE_centroids_um[:, 0], c='r', s=2)
            b.scatter(MxIF_centroids_um[:, 1], MxIF_centroids_um[:, 0], c='r', s=2)
            # a.axis('off')
            # b.axis('off')
            plt.show()
            print("H&E Cell Count: %d" % len(HE_centroids_pix))
            print("MxIF Cell Count: %d" % len(MxIF_centroids_pix))

        HE_selection, HE_nodes_loc = get_points_in_window(trans_HE_centroids_um, location=(520, 620), window_size=50)
        HE_nodes_features = HE_morph_features[HE_selection]
        HE_graph, HE_Adj_M = create_graph_for_point_set(HE_nodes_loc)
        search_range = 1
        q_MxIF_nodes_dis, q_MxIF_idx = create_graph_use_query(HE_nodes_loc, MxIF_centroids_um, K=search_range)

        MxIF_selection = list(set(np.array(q_MxIF_idx).flatten()))
        MxIF_nodes_loc = MxIF_centroids_um[MxIF_selection]
        MxIF_nodes_features = MxIF_morph_features[MxIF_selection]

        plt.figure(1)
        plt.scatter(HE_nodes_loc[:, 1], HE_nodes_loc[:, 0], c='r', s=10)
        plt.scatter(MxIF_nodes_loc[:, 1], MxIF_nodes_loc[:, 0], c='b', s=10)
        plt.legend(["HE", "MxIF"])
        plt.show()

        ###############################################3
        ## use CPD again, dosen't work very well
        tf_param_A, sigma2_A, q_A = cpd.registration_cpd(HE_nodes_loc, MxIF_nodes_loc, maxiter=20, tol=1.0e-11,
                                                   update_scale=False)
        trans_HE_nodes_loc = tf_param_A.transform(copy.deepcopy(HE_nodes_loc))
        plt.figure(1)
        plt.scatter(trans_HE_nodes_loc[:, 1], trans_HE_nodes_loc[:, 0], c='r', s=10)
        plt.scatter(MxIF_nodes_loc[:, 1], MxIF_nodes_loc[:, 0], c='b', s=10)
        plt.legend(["HE", "MxIF"])
        plt.show()
        ###############################################3
        print("Fine tune with graph matching")
        MxIF_graph, MxIF_Adj_M = create_graph_for_point_set(MxIF_nodes_loc)

        pos1 = {_: HE_nodes_loc[_, :] for _ in range(len(HE_nodes_loc))}
        pos2 = {_: MxIF_nodes_loc[_, :] for _ in range(len(MxIF_nodes_loc))}
        color1 = ['#FF5733' for _ in range(len(HE_nodes_loc))]
        color2 = ['#1f78b4' for _ in range(len(MxIF_nodes_loc))]

        draw_two_graphs(HE_graph, MxIF_graph, pos1, pos2, color1, color2)
        ################## show graph in the same plot
        nx.draw_networkx(HE_graph, pos=pos1, edge_color='r', node_color=color1, with_labels=True, node_size=15)
        nx.draw_networkx(MxIF_graph, pos=pos2, edge_color='b', node_color=color2, with_labels=True, node_size=5)
        r_patch = mpatches.Patch(color='red', label="HE")
        b_patch = mpatches.Patch(color='blue', label="MxIF")
        plt.legend(handles=[r_patch, b_patch])
        plt.show()
        ##################
        X = GM_match(HE_Adj_M, MxIF_Adj_M, HE_nodes_features, MxIF_nodes_features)
        ### get matching pairs
        source_match_idx = range(len(HE_nodes_loc))
        sorted_MxIF_nodes_loc, sorted_MxIF_nodes_idx, target_match_idx = get_matching_pairs(source_match_idx, X, MxIF_nodes_loc, MxIF_selection)
        ###  filtering the matching pairs
        mask = LPM_filter(HE_nodes_loc, sorted_MxIF_nodes_loc)

        ### Draw the matching links before filtering
        draw_graph_with_matching_links(HE_graph, MxIF_graph,
                                       pos1, pos2, color1, color2,
                                       source_match_idx, target_match_idx, filter=None)

        ## draw the links between matching points after filtering
        draw_graph_with_matching_links(HE_graph, MxIF_graph,
                                       pos1, pos2, color1, color2,
                                       source_match_idx, target_match_idx, filter=mask)

        src_points = HE_centroids_um[HE_selection][mask]
        dst_points = sorted_MxIF_nodes_loc[mask]
        M = cv2.estimateAffinePartial2D(src_points, dst_points)[0]
        GM_M = get_M_from_cv2_affine(M, HE_pixel_size, MxIF_pixel_size)
        use_GM_results = check_CPD_GM(CPD_M, GM_M)  # check the affine transformation calculated from graph matching, if it's too different from the result from CPD
        if sum(mask) > 2:
            if use_GM_results:
                final_M = GM_M
            else:
                final_M = CPD_M  # unable to fine tune, use CPD transformation
        else:
            final_M = CPD_M
        Aff_CPD_M_fn = os.path.join(output_dir, tif_fn + "_CPD_M.npy")
        np.save(Aff_CPD_M_fn, CPD_M)
        np.save(Aff_M_fn, final_M)  # save the final affine transformation matrix
        ################################################################
        affine_HE_fn = os.path.join(output_dir, tif_fn)

        dapi_img = mxif_img.pages[0].asarray().astype(float)
        trans_he_img = save_transformed_HE(he_img, final_M, dapi_img.shape, MxIF_pixel_size, affine_HE_fn)

        if verbose:
            s_fn = os.path.join(output_dir, "side_by_side_" + tif_fn)
            save_transformed_HE_and_MxIF(trans_he_img, mxif_img_3_channel, dapi_img.shape, MxIF_pixel_size, s_fn)

    print("Done")
