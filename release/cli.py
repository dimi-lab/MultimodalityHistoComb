import os
import argparse
from probreg import cpd
from utils import *
import tifffile as tf
import copy

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
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        default=False,
                        help="Save debug information"
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

    # MxIF_quant_fn = "/opt/export/A-8_StarDist_QUANT.tsv"
    # HE_quant_fn = "/opt/export/A-8_StarDist_QUANT.tsv"
    # HE_img_fn = "/opt/HE_FOVs/same_section/G-8.tif"
    # MxIF_img_fn = "/opt/MxIF_FOVs/Slide2050_24Plex/G-8.tif"
    # output_dir = "/opt/test_align_wsi"
    # verbose = True

    tif_fn = os.path.split(HE_img_fn)[1]

    HE_centroids = get_cell_loc(HE_quant_fn)   # you may need to specify column names to get the cell coordinates
    MxIF_centroids = get_cell_loc(MxIF_quant_fn)

    tf_param, sigma2, q = cpd.registration_cpd(HE_centroids, MxIF_centroids, maxiter=20, tol=1.0e-11,
                                               update_scale=False)
    r_points = tf_param.transform(copy.deepcopy(HE_centroids))
    if verbose:
        title = "Cell Centroids Before and after alignment"
        legend = ["HE cells", "MxIF cells", "trans_HE"]
        fn = "log_" + tif_fn + "_centroids_alignment.png"
        plot_centroids_trans(HE_centroids, MxIF_centroids, r_points, legend, title, output_dir, fn)

    # TODO: read pixel size from the image tags (assigned to Raymond)
    HE_pixel_size = 0.2201  # unit micron
    MxIF_pixel_size = 0.325
    real_S = HE_pixel_size/MxIF_pixel_size
    theta, degrees, s, delta, gt_M = get_M_from_cpd(tf_param, real_S)

    # TODO: add refinement after CPD  (assigned to Jun)
    affine_HE_fn = os.path.join(output_dir, tif_fn)
    he_img = tf.TiffFile(HE_img_fn).pages[0].asarray().astype(np.float)  # read the images
    mxif_img = tf.TiffFile(MxIF_img_fn)
    dapi_img = mxif_img.pages[0].asarray().astype(np.float)
    trans_he_img = save_transformed_HE(he_img, gt_M, dapi_img.shape, MxIF_pixel_size, affine_HE_fn)

    if verbose:
        s_fn = os.path.join(output_dir, "side_by_side_"+tif_fn)
        mxif_img_3_channel1 = mxif_img.pages[0].asarray().astype(np.uint8)
        mxif_img_3_channel2 = mxif_img.pages[1].asarray().astype(np.uint8)
        mxif_img_3_channel3 = mxif_img.pages[2].asarray().astype(np.uint8)
        mxif_img_3_channel = np.stack([mxif_img_3_channel1, mxif_img_3_channel2, mxif_img_3_channel3], axis=2)
        save_transformed_HE_and_MxIF(trans_he_img, mxif_img_3_channel, dapi_img.shape, MxIF_pixel_size, s_fn)























