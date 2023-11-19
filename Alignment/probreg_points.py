import os
from path_config import my_path
from probreg import cpd, gmmtree, filterreg, bcpd
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_cell_loc_from_df(HE_quant_df, MxIF_quant_df):
    he_x = HE_quant_df["Centroid X µm"]
    he_y = HE_quant_df["Centroid Y µm"]
    mxif_x = MxIF_quant_df["Centroid X µm"]
    mxif_y = MxIF_quant_df["Centroid Y µm"]

    source = np.array([he_x, he_y]).T
    target = np.array([mxif_x, mxif_y]).T
    return source, target

def align_cell_segmentation(source, target):
    # compute cpd registration
    tf_param, a, b = cpd.registration_cpd(source, target, maxiter=150)
    # tf_param = bcpd.registration_bcpd(source, target, maxiter=150)

    result = copy.deepcopy(source)
    r_points = tf_param.transform(result)

    return tf_param, r_points

def generate_core_list(row_id_range:list, column_id_range: list):
    core_list = []
    for r in row_id_range:
        for c in column_id_range:
            core_list.append(r + "-" + str(c))
    return core_list

if __name__ == "__main__":
    output_alignment_dir = my_path.output_alignment
    HE_data_dir = my_path.HE_FOV_export_dir
    MxIF_data_dir = my_path.MxIF_FOV_export_dir

    ROI_list = generate_core_list(["A", "B"], list(range(1, 22)))
    for roi_id in ROI_list:
        HE_quant_fn = os.path.join(HE_data_dir, roi_id, roi_id+"_QUANT.tsv")
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
            source, target = get_cell_loc_from_df(HE_quant_df, MxIF_quant_df)
            tf_param, r_points = align_cell_segmentation(source, target)

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