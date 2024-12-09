import os
import argparse
import random
from probreg import cpd
from utils import *
import tifffile as tf
import copy
from gm_utils import *
import matplotlib.pyplot as plt
from LPM import *

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
    parser.add_argument("-is", "--he_img_fn",
                        required=True,
                        dest="he_img_fn",
                        help="source image, H&E image file, supposed to be ome tiff")
    parser.add_argument("-it", "--mxif_img_fn",
                        required=True,
                        dest="mxif_img_fn",
                        help="target image (MxIF image file), supposed to be ome tiff")
    parser.add_argument("-v", "--verbose", type=str2bool, nargs='?',
                        const=True, default=False,
                        help="Show debug information"
                        )
    parser.add_argument("-o", "--output_dir",
                        default=os.getcwd(),
                        dest='output_dir',
                        help="Transformed H&E output directory")

    args = parser.parse_args()
    MxIF_quant_fn = args.mxif_quant_fn
    HE_quant_fn = args.he_quant_fn
    HE_img_fn = args.he_img_fn
    MxIF_img_fn = args.mxif_img_fn
    output_dir = args.output_dir
    verbose = args.verbose

    # HE_img_fn = "/temp/Ovarian_TMA/HE_FOVs/same_section/A-8.tif"
    # MxIF_img_fn = "/temp/Ovarian_TMA/MxIF_FOVs/Slide2050_24Plex/A-8.tif"
    # MxIF_quant_fn = "/temp/Ovarian_TMA/AlignmentEval/QuPathAnnoProj_MxIF/export/A-8_StarDist_QUANT.tsv"
    # HE_quant_fn = "/temp/Ovarian_TMA/AlignmentEval/QuPathAnnoProj_HE_Sec1/export/A-8_StarDist_QUANT.tsv"
    # output_dir = "/temp/Ovarian_TMA/test_align_seg_qupath"
    # verbose = True

    # TODO: read pixel size from the image tags
    HE_pixel_size = 0.2201  # unit micron
    MxIF_pixel_size = 0.325

    tif_fn = os.path.split(HE_img_fn)[1]
    print("Processing %s" % tif_fn)
    Aff_M_fn = os.path.join(output_dir, tif_fn + "_final_M.npy")
    Aff_CPD_M_fn = os.path.join(output_dir, tif_fn + "_CPD_M.npy")
    if not os.path.exists(Aff_M_fn) or not os.path.exists(Aff_CPD_M_fn):
        print("\t Loading cell segmentation.")
        HE_centroids_um, HE_morph_features = get_cell_info_from_QuPath(HE_quant_fn)   # you may need to specify column names to get the cell coordinates
        MxIF_centroids_um, MxIF_morph_features = get_cell_info_from_QuPath(MxIF_quant_fn)
        print("\t Aligning with CPD.")
        tf_param, sigma2, q = cpd.registration_cpd(HE_centroids_um, MxIF_centroids_um, maxiter=20, tol=1.0e-11,
                                                   update_scale=False)
        trans_HE_centroids_um = tf_param.transform(copy.deepcopy(HE_centroids_um))
        if verbose:
            title = "Cell Centroids Before and after alignment"
            legend = ["HE cells", "MxIF cells", "trans_HE"]
            fn = "log_" + tif_fn + "_centroids_alignment.png"
            plot_centroids_trans(HE_centroids_um, MxIF_centroids_um, trans_HE_centroids_um, legend, title, output_dir, fn)

        theta, degrees, s, delta, CPD_M = get_M_from_cpd(tf_param, HE_pixel_size, MxIF_pixel_size)

        #  refinement after CPD
        sample_wind = 50
        dense_scores = KDE_cell_density(trans_HE_centroids_um)
        sample_pos_options_idx = np.array(dense_scores) > 0.75 * max(dense_scores)
        sample_pos_options = trans_HE_centroids_um[sample_pos_options_idx]
        sample_i = random.choice(range(len(sample_pos_options)))
        sample_loc = tuple(sample_pos_options[sample_i] - sample_wind / 2)

        HE_selection, HE_nodes_loc = get_points_in_window(trans_HE_centroids_um, location=sample_loc,
                                                          window_size=sample_wind)
        HE_nodes_features = HE_morph_features[HE_selection]
        HE_graph, HE_Adj_M = create_graph_for_point_set(HE_nodes_loc)
        search_range = 1
        q_MxIF_nodes_dis, q_MxIF_idx = create_graph_use_query(HE_nodes_loc, MxIF_centroids_um, K=search_range)

        MxIF_selection = list(set(np.array(q_MxIF_idx).flatten()))
        MxIF_nodes_loc = MxIF_centroids_um[MxIF_selection]
        MxIF_nodes_features = MxIF_morph_features[MxIF_selection]

        print("\t Fine tune with graph matching")
        MxIF_graph, MxIF_Adj_M = create_graph_for_point_set(MxIF_nodes_loc)

        pos1 = {_: HE_nodes_loc[_, :] for _ in range(len(HE_nodes_loc))}
        pos2 = {_: MxIF_nodes_loc[_, :] for _ in range(len(MxIF_nodes_loc))}
        color1 = ['#FF5733' for _ in range(len(HE_nodes_loc))]
        color2 = ['#1f78b4' for _ in range(len(MxIF_nodes_loc))]

        X = GM_match(HE_Adj_M, MxIF_Adj_M, HE_nodes_features, MxIF_nodes_features)
        ### get matching pairs
        source_match_idx = range(len(HE_nodes_loc))
        sorted_MxIF_nodes_loc, sorted_MxIF_nodes_idx, target_match_idx = get_matching_pairs(source_match_idx, X,
                                                                                            MxIF_nodes_loc, MxIF_selection)
        ###  filtering the matching pairs
        mask = LPM_filter(HE_nodes_loc, sorted_MxIF_nodes_loc)

        if sum(mask) > 2:
            src_points = HE_centroids_um[HE_selection][mask]
            dst_points = sorted_MxIF_nodes_loc[mask]
            M = cv2.estimateAffinePartial2D(src_points, dst_points)[0]
            GM_M = get_M_from_cv2_affine(M, HE_pixel_size, MxIF_pixel_size)
            # check the affine transformation calculated from graph matching, if it's too different from the result from CPD
            use_GM_results = check_CPD_GM(CPD_M, GM_M)
            if use_GM_results:
                final_M = GM_M
            else:
                final_M = CPD_M  # unable to fine tune, use CPD transformation
        else:
            final_M = CPD_M

        np.save(Aff_CPD_M_fn, CPD_M)
        np.save(Aff_M_fn, final_M)  # save the final affine transformation matrix

        if verbose:
            plt.figure(1)
            plt.scatter(HE_nodes_loc[:, 1], HE_nodes_loc[:, 0], c='r', s=10)
            plt.scatter(MxIF_nodes_loc[:, 1], MxIF_nodes_loc[:, 0], c='b', s=10)
            plt.legend(["HE", "MxIF"])
            plt.show()

        affine_HE_fn = os.path.join(output_dir, tif_fn)
        he_img = tf.TiffFile(HE_img_fn).pages[0].asarray().astype(float)  # read the images
        mxif_img = tf.TiffFile(MxIF_img_fn)
        dapi_img = mxif_img.pages[0].asarray().astype(float)

        trans_he_img = save_transformed_HE(he_img, final_M, dapi_img.shape, MxIF_pixel_size, affine_HE_fn)

        if verbose:
            s_fn = os.path.join(output_dir, "side_by_side_"+tif_fn)
            mxif_img_3_channel1 = mxif_img.pages[0].asarray().astype(np.uint8)
            mxif_img_3_channel2 = mxif_img.pages[1].asarray().astype(np.uint8)
            mxif_img_3_channel3 = mxif_img.pages[2].asarray().astype(np.uint8)
            mxif_img_3_channel = np.stack([mxif_img_3_channel1, mxif_img_3_channel2, mxif_img_3_channel3], axis=2)
            save_transformed_HE_and_MxIF(trans_he_img, mxif_img_3_channel, dapi_img.shape, MxIF_pixel_size, s_fn)


    print("Done")




















