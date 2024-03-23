import os
import argparse
from probreg import cpd
from utils import *
import tifffile as tf
import copy
from stardist.models import StarDist2D, Config2D
from csbdeep.utils import Path, normalize

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    #
    # parser.add_argument("-is", "--he_img_fn",
    #                     required=True,
    #                     dest="he_img_fn",
    #                     help="source image, H&E image file, supposed to be ome tiff")
    # parser.add_argument("-it", "--mxif_img_fn",
    #                     required=True,
    #                     dest="mxif_img_fn",
    #                     help="target image (MxIF image file), supposed to be ome tiff")
    # parser.add_argument("-v", "--verbose",
    #                     dest="verbose",
    #                     default=False,
    #                     help="Save debug information"
    #                     )
    # parser.add_argument("-o", "--output_dir",
    #                     default=os.getcwd(),
    #                     dest='output_dir',
    #                     help="Transformed H&E output directory")
    #
    # args = parser.parse_args()
    #
    # HE_img_fn = args.he_img_fn
    # MxIF_img_fn = args.mxif_img_fn
    # output_dir = args.output_dir
    # verbose = args.verbose

    # MxIF_model_fn = "/infodev1/non-phi-data/junjiang/model/StarDist/dsb2018_heavy_augment.pb"
    # HE_model_fn = "/infodev1/non-phi-data/junjiang/model/StarDist/he_heavy_augment.pb"

    HE_img_fn = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/HE_FOVs/same_section/A-8.tif"
    MxIF_img_fn = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/MxIF_FOVs/Slide2050_24Plex/A-8.tif"
    output_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/test_align_wsi"
    verbose = True

    tif_fn = os.path.split(HE_img_fn)[1]
    # TODO: read pixel size from the image tags (assigned to Raymond)
    HE_pixel_size = 0.2201  # unit micron
    MxIF_pixel_size = 0.325

    #################################### Do cell segmentation #################################################
    he_cell_cent_fn = os.path.join(output_dir, tif_fn + "_he_cell_centroids.npy")
    mxif_cell_cent_fn = os.path.join(output_dir, tif_fn + "_mxif_cell_centroids.npy")

    if not os.path.exists(he_cell_cent_fn) or (not os.path.exists(mxif_cell_cent_fn)):
        print("Getting models from the list")
        StarDist2D.from_pretrained()

    # read the images
    he_img = tf.TiffFile(HE_img_fn).pages[0].asarray().astype(float)
    norm_he_img = normalize(he_img)
    if not os.path.exists(he_cell_cent_fn):
        he_model = StarDist2D.from_pretrained('2D_versatile_he')
        he_lbl, _ = he_model.predict_instances_big(norm_he_img, prob_thresh=0.5, axes='YXC', block_size=2048,
                                                   min_overlap=128, context=128)
        HE_centroids_pix = np.array(get_cell_loc_from_StarDist_pred(he_lbl))
        np.save(he_cell_cent_fn, HE_centroids_pix)
    else:
        HE_centroids_pix = np.load(he_cell_cent_fn)

    mxif_img = tf.TiffFile(MxIF_img_fn)

    mxif_img_3_channel1 = mxif_img.pages[0].asarray().astype(np.uint8)  # DAPI channel
    mxif_img_3_channel2 = mxif_img.pages[1].asarray().astype(np.uint8)
    mxif_img_3_channel3 = mxif_img.pages[2].asarray().astype(np.uint8)

    mxif_img_3_channel = np.stack([mxif_img_3_channel1, mxif_img_3_channel2, mxif_img_3_channel3], axis=2)
    if not os.path.exists(mxif_cell_cent_fn):
        mxif_model = StarDist2D.from_pretrained('2D_versatile_fluo')
        norm_mxif_img = normalize(mxif_img_3_channel1.astype(float), 1, 99.8, axis=(0, 1))
        mxif_lbl, _ = mxif_model.predict_instances_big(norm_mxif_img, axes='YX', block_size=2048, min_overlap=128,
                                                       context=128)
        MxIF_centroids_pix = np.array(get_cell_loc_from_StarDist_pred(mxif_lbl))
        np.save(mxif_cell_cent_fn, MxIF_centroids_pix)
    else:
        MxIF_centroids_pix = np.load(mxif_cell_cent_fn)

    if verbose:
        fig, (a, b) = plt.subplots(1, 2, figsize=(16, 8))
        a.scatter(HE_centroids_pix[:, 1], HE_centroids_pix[:, 0], c='r', s=2)
        b.scatter(MxIF_centroids_pix[:, 1], MxIF_centroids_pix[:, 0], c='r', s=2)
        a.axis('off')
        b.axis('off')
        plt.show()
        print("H&E Cell Count: %d" % len(HE_centroids_pix))
        print("MxIF Cell Count: %d" % len(MxIF_centroids_pix))

    ##################################### Align based on cell segmentation ########################################
    print("Aligning based on cell segmentation.")
    HE_centroids_um = np.array(HE_centroids_pix) * HE_pixel_size
    MxIF_centroids_um = np.array(MxIF_centroids_pix) * MxIF_pixel_size
    # HE_centroids_um = np.array(HE_centroids_pix) * MxIF_pixel_size
    # MxIF_centroids_um = np.array(MxIF_centroids_pix) * HE_pixel_size
    # HE_centroids_um = np.array(HE_centroids_pix)
    # MxIF_centroids_um = np.array(MxIF_centroids_pix)
    tf_param, sigma2, q = cpd.registration_cpd(HE_centroids_um, MxIF_centroids_um, maxiter=20, tol=1.0e-11,
                                               update_scale=False)
    r_points = tf_param.transform(copy.deepcopy(HE_centroids_um))
    if verbose:
        title = "Cell Centroids Before and after alignment"
        legend = ["HE cells", "MxIF cells", "trans_HE"]
        fn = "log_" + tif_fn + "_centroids_alignment.png"
        plot_centroids_trans(HE_centroids_um, MxIF_centroids_um, r_points, legend, title, output_dir, fn)

    # real_S = HE_pixel_size/MxIF_pixel_size
    # real_S = 1
    theta, degrees, s, delta, gt_M = get_M_from_cpd(tf_param, HE_pixel_size, MxIF_pixel_size)

    # TODO: add refinement after CPD  (assigned to Jun)
    affine_HE_fn = os.path.join(output_dir, tif_fn)

    dapi_img = mxif_img.pages[0].asarray().astype(np.float)
    trans_he_img = save_transformed_HE(he_img, gt_M, dapi_img.shape, MxIF_pixel_size, affine_HE_fn)

    if verbose:
        s_fn = os.path.join(output_dir, "side_by_side_" + tif_fn)

        save_transformed_HE_and_MxIF(trans_he_img, mxif_img_3_channel, dapi_img.shape, MxIF_pixel_size, s_fn)
