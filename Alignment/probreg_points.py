import os
from path_config import my_path
from probreg import cpd, gmmtree, filterreg, bcpd, math_utils
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pickle
from scipy.spatial import KDTree
import tifffile as tf
import cv2

def get_cell_loc_from_df(HE_quant_df, MxIF_quant_df):
    he_x = HE_quant_df["Centroid X µm"]
    he_y = HE_quant_df["Centroid Y µm"]
    mxif_x = MxIF_quant_df["Centroid X µm"]
    mxif_y = MxIF_quant_df["Centroid Y µm"]

    source = np.array([he_x, he_y]).T
    target = np.array([mxif_x, mxif_y]).T
    return source, target

def norm(rvalue, newmin, newmax):
    oldmin = min(rvalue)
    oldmax = max(rvalue)
    oldrange = oldmax - oldmin
    newrange = newmax - newmin
    if oldrange == 0:  # Deal with the case where rvalue is constant:
        if oldmin < newmin:  # If rvalue < newmin, set all rvalue values to newmin
            newval = newmin
        elif oldmin > newmax:  # If rvalue > newmax, set all rvalue values to newmax
            newval = newmax
        else:  # If newmin <= rvalue <= newmax, keep rvalue the same
            newval = oldmin
        normal = [newval for _ in rvalue]
    else:
        scale = newrange / oldrange
        normal = [(v - oldmin) * scale + newmin for v in rvalue]
    return normal

def KDE_cell_density(cell_xy):
    xy = np.vstack([cell_xy[:, 0], cell_xy[:, 1]])
    kde_scores = gaussian_kde(xy)(xy)
    norm_kde_scores = norm(kde_scores, 0, 1)
    return norm_kde_scores
def get_specific_feature(df, idx_start, idx_end=None):
    if idx_end is None:
        val = np.array(df.iloc[:, idx_start:])
        names = df.iloc[:, idx_start:].columns
    else:
        val = np.array(df.iloc[:, idx_start:idx_end])
        names = df.iloc[:, idx_start:idx_end].columns
    return names, val

def get_cell_loc_with_density(HE_quant_df, MxIF_quant_df):
    col = list(HE_quant_df.columns)
    df = HE_quant_df.dropna(subset=["Image"] + col[7:])
    xy_names, he_cell_loc = get_specific_feature(df, 7, 9)
    he_norm_kde_scores = KDE_cell_density(he_cell_loc)

    col = list(MxIF_quant_df.columns)
    df = MxIF_quant_df.dropna(subset=["Image"] + col[6:])
    xy_names, mxif_cell_loc = get_specific_feature(df, 6, 8)
    mxif_norm_kde_scores = KDE_cell_density(mxif_cell_loc)

    source = np.array([he_cell_loc[:, 0], he_cell_loc[:, 1], np.array(he_norm_kde_scores)*300]).T
    target = np.array([mxif_cell_loc[:, 0], mxif_cell_loc[:, 1], np.array(mxif_norm_kde_scores)*300]).T
    return source, target

def align_cell_segmentation(source, target):
    # compute cpd registration
    # TODO: tf_param, _, _ = cpd.registration_cpd(source, target, update_scale=False)
    # TODO: set update_scale to False, as the cell locations were denoted with microns (calculated based on pixel size)
    # TODO: But how to apply the affine transformation obtained from this registration if the pixel sizes are different?
    # TODO: sx=sy=source_pixel_size/target_pixel_size
    # TODO:
    # tf_param, sigma2, q = cpd.registration_cpd(source, target, maxiter=150)
    tf_param, sigma2, q = cpd.registration_cpd(source, target, maxiter=150, update_scale=False)
    # tf_param = bcpd.registration_bcpd(source, target, maxiter=150)

    result = copy.deepcopy(source)
    r_points = tf_param.transform(result)

    return tf_param, r_points, sigma2, q

def generate_core_list(row_id_range:list, column_id_range: list):
    core_list = []
    for r in row_id_range:
        for c in column_id_range:
            core_list.append(r + "-" + str(c))
    return core_list

def get_TMA_core_list(img_path:str) -> list:
    img_fn_list = os.listdir(img_path)
    roi_list = [i.split(".")[0] for i in img_fn_list]
    return roi_list

if __name__ == "__main__":
    output_alignment_dir = my_path.output_alignment_same_sec
    # HE_data_dir = my_path.HE_FOV_export_dir
    same_HE_data_dir = my_path.same_HE_FOV_export_dir
    MxIF_data_dir = my_path.MxIF_FOV_export_dir

    #TODO: how to get the transformed image based on transformation.
    HE_img_data_dir = my_path.he_image_data_sec1
    MxIF_img_data_dir = my_path.mxif_image_data_sec1
    # TODO: read MxIF and H&E image
    # TODO:

    # ROI_list = generate_core_list(["A", "B"], list(range(1, 23)))
    ROI_list = get_TMA_core_list(MxIF_img_data_dir)

    HE_pxiel_size = 0.2201  # unit micron
    MxIF_pixel_size = 0.325

    pix_scale = HE_pxiel_size/MxIF_pixel_size
    target_pixel_size = MxIF_pixel_size

    eval_fn = os.path.join(output_alignment_dir, "eval_metrics.csv")
    wrt_str2eval = "roi_id,source,target,source_point_num,target_point_num,sigma2,q,KDTree_rmse_before,KDTree_rmse_after\n"
    eval_fp = open(eval_fn, 'w')
    for roi_id in ROI_list:
        # roi_id ='A-22'
        # roi_id = 'A-8'
        trans_fn = os.path.join(output_alignment_dir, roi_id+"_trans.dat")
        HE_quant_fn = os.path.join(same_HE_data_dir, roi_id, roi_id + "_QUANT.tsv")
        MxIF_quant_fn = os.path.join(MxIF_data_dir, roi_id, roi_id + "_QUANT.tsv")
        print(HE_quant_fn)
        print(MxIF_quant_fn)

        if not os.path.exists(HE_quant_fn) or not os.path.exists(MxIF_quant_fn):
            if os.path.exists(HE_quant_fn):
                print(HE_quant_fn)
            if os.path.exists(MxIF_quant_fn):
                print(MxIF_quant_fn)
            print("file does not exist: %s" % roi_id)
            continue
        else:
            print("processing: %s" % roi_id)
            HE_quant_df = pd.read_csv(HE_quant_fn, sep='\t')
            MxIF_quant_df = pd.read_csv(MxIF_quant_fn, sep='\t')
            source, target = get_cell_loc_from_df(HE_quant_df, MxIF_quant_df) # without cell density
            print("Number of points: (%d, %d)" % (len(source), len(target)))

            # Get affine transformation from points (cell centroids)
            # source, target = get_cell_loc_with_density(HE_quant_df, MxIF_quant_df) # with cell density append to cell locations
            tf_param, r_points, sigma2, q = align_cell_segmentation(source, target)

            # TODO: apply affine transformation to the tissue core
            #
            he_fn = os.path.join(HE_img_data_dir, roi_id + ".tif")
            he_img = tf.TiffFile(he_fn).pages[0].asarray().astype(np.float)
            mxif_fn = os.path.join(MxIF_img_data_dir, roi_id + ".tif")
            mxif_img = tf.TiffFile(mxif_fn)

            dapi_img = mxif_img.pages[0].asarray().astype(np.float)

            R = tf_param.rot
            T = tf_param.t
            M = np.array([[pix_scale*R[0, 0], R[0,1], T[0]/target_pixel_size],
                          [R[1, 0], pix_scale*R[1, 1], T[1]/target_pixel_size]]).astype(float)


            affined_image = cv2.warpAffine(src=he_img, M=M, dsize=dapi_img.shape)

            # plt.imshow(he_img.astype(np.uint8))
            # plt.show()
            #
            # plt.imshow(affined_image.astype(np.uint8))
            # plt.show()
            #
            # plt.imshow(dapi_img, cmap="gray")
            # plt.show()

            res_img = affined_image.astype(np.uint8)
            affine_HE_sec1_fn = os.path.join(my_path.HE_affine_output_sec1, roi_id+"_trans.tif")

            resolution = (10000 / target_pixel_size, 10000 / target_pixel_size)  # convert um to cm for saving
            tf.imwrite(affine_HE_sec1_fn, res_img, photometric='rgb', resolution=resolution, resolutionunit="CENTIMETER")

            # calculate loss, save later
            target_tree = KDTree(target, leafsize=10)
            rmse_after = math_utils.compute_rmse(r_points, target_tree)
            rmse_before = math_utils.compute_rmse(source, target_tree)
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
            wrt_str2eval += roi_id + ",HE,MxIF," + str(len(source)) + "," + str(len(target)) + "," + str(sigma2) + "," + str(q) + "," + str(rmse_before) + "," + str(rmse_after) + "\n"

            data = [tf_param.rot, tf_param.scale, tf_param.t]
            with open(trans_fn, "wb") as f:
                pickle.dump(data, f)
            # with open(trans_fn, "rb") as f:
            #     a = pickle.load(f)
            #     print(a)


            # Y = np.dot(R, X) + t
            # new_xy = (np.dot(tf_param.rot, xy) + tf_param.t) * tf_param.scale
            ## https://opencv-tutorial.readthedocs.io/en/latest/trans/transform.html#rotation

            # create plot to show the points before and after alignment
            plt.figure(0, figsize=(8, 8))
            plt.scatter(source[:,0], source[:,1], s=1)
            plt.scatter(target[:,0], target[:,1], s=1)
            plt.axis('equal')
            plt.title("Before Transformation")
            plt.xlim((0, 1000))
            plt.ylim((0, 1000))
            plt.xlabel("x location: µm")
            plt.ylabel("y location: µm")
            plt.legend(["H&E cells", "MxIF cells"])
            plt.tight_layout()
            sv_fn = os.path.join(output_alignment_dir, roi_id+"_before_trans.png")
            plt.savefig(sv_fn, dpi=200)
            # plt.show()
            plt.close()

            plt.figure(1, figsize=(8, 8))
            plt.scatter(r_points[:,0], r_points[:,1], s=1)
            plt.scatter(target[:,0], target[:,1], s=1)
            plt.axis('equal')
            plt.title("After Transformation")
            plt.xlim((0, 1000))
            plt.ylim((0, 1000))
            plt.legend(["H&E cells", "MxIF cells"])
            plt.xlabel("x location: µm")
            plt.ylabel("y location: µm")
            plt.tight_layout()
            sv_fn = os.path.join(output_alignment_dir, roi_id+"_after_trans.png")
            plt.savefig(sv_fn, dpi=200)
            # plt.show()
            plt.close()


    eval_fp.write(wrt_str2eval)
    eval_fp.close()

print("x,y")

# plt.scatter(source[:,0], source[:,1], s=1)
# plt.axis('equal')
# plt.title("H&E cells")
# plt.show()
#
# plt.scatter(target[:,0], target[:,1], s=1)
# plt.axis('equal')
# plt.title("MxIF cells")
# plt.show()