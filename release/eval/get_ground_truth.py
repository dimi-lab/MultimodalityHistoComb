import os
import numpy as np
import cv2
import math
import matplotlib.pyplot as plt
import pickle
import glob
import itertools
from eval_utils import *
from sys import platform


Sec = 1

if platform == "linux" or platform == "linux2":
    data_root_dir = "/temp/Ovarian_TMA/AlignmentEval"
    if Sec == 1:
        # Sec1
        ground_truth_output_dir = os.path.join(data_root_dir, "GroundTruthEvaluation", "GT_HE_Sec1_MxIF")
        HE_export_dir = os.path.join(data_root_dir, "GroundTruthAnnotation", "HE_Sec1")
    elif Sec == 2:
        # Sec2
        ground_truth_output_dir = os.path.join(data_root_dir, "GroundTruthEvaluation", "GT_HE_Sec2_MxIF")
        HE_export_dir = os.path.join(data_root_dir, "GroundTruthAnnotation", "HE_Sec2")
    else:
        raise Exception("Undefined Section")
elif platform == "win32":
    data_root_dir = r"\\mfad\researchmn\HCPR\HCPR-GYNECOLOGICALTUMORMICROENVIRONMENT\Archive\WSIs\Ovarian_TMA\AlignmentEval"
    if Sec == 1:
        # Sec1
        ground_truth_output_dir = os.path.join(data_root_dir, "GroundTruthEvaluation", "GT_HE_Sec1_MxIF")
        HE_export_dir = os.path.join(data_root_dir, "GroundTruthAnnotation", "HE_Sec1")
    elif Sec == 2:
        # Sec2
        ground_truth_output_dir = os.path.join(data_root_dir, "GroundTruthEvaluation", "GT_HE_Sec2_MxIF")
        HE_export_dir = os.path.join(data_root_dir, "GroundTruthAnnotation", "HE_Sec2")
    else:
        raise Exception("Undefined Section")
else:
    raise Exception("not defined")

MxIF_export_dir = os.path.join(data_root_dir, "GroundTruthAnnotation", "MxIF_Sec1")

# file_name_list_with_path = glob.glob(os.path.join(HE_export_dir, "*align_anno_8.csv"))
# file_name_list = [os.path.split(i)[1] for i in file_name_list_with_path]

core_list_fn = "../core_list.txt"
core_list = open(core_list_fn, 'r').readlines()
file_name_list = []
for i in core_list:
    file_name = os.path.join(HE_export_dir, i.strip().replace(".tif", "_align_anno_8.csv"))
    file_name_list.append(file_name)

# the export landmarks were pixel-based rather than um-based
HE_pixel_size = 0.2201  # unit micron
MxIF_pixel_size = 0.325
#########################################
# for debug
# DEBUG = True
DEBUG = False
DEBUG_1 = True
#########################################

for fn in file_name_list:
    print("Processing %s" % fn)
    csv_fn = os.path.split(fn)[1]
    roi_id = csv_fn.split("_")[0]
    # Please note: The point location is exported as pixels, not microns
    HE_landmarks_fn = os.path.join(HE_export_dir, csv_fn)
    MxIF_landmarks_fn = os.path.join(MxIF_export_dir, csv_fn)
    point_sort_fn = os.path.join(data_root_dir, "GroundTruthAnnotation", "Landmark_correspondence", "Sec%d" %Sec, roi_id+".csv")

    sorted_HE_landmarks, sorted_MxIF_landmarks = get_sorted_annotation_landmark_pairs(HE_landmarks_fn, MxIF_landmarks_fn, point_sort_fn)

    if DEBUG_1:
        fig, ax = plt.subplots(1, 2)
        ax[0].scatter(sorted_HE_landmarks[:, 0], sorted_HE_landmarks[:, 1], c='r', s=8)
        for i in range(len(sorted_HE_landmarks)):
            ax[0].annotate(str(i), (sorted_HE_landmarks[i, 0], sorted_HE_landmarks[i, 1]))
        plt.title("HE landmarks:%s" % roi_id)
        # sorted_MxIF_landmarks = MxIF_landmarks

        ax[1].scatter(sorted_MxIF_landmarks[:, 0], sorted_MxIF_landmarks[:, 1], c='b', s=8)
        for i in range(len(sorted_MxIF_landmarks)):
            ax[1].annotate(str(i), (sorted_MxIF_landmarks[i, 0], sorted_MxIF_landmarks[i, 1]))
        plt.title("MxIF landmarks:%s" % roi_id)
        plt.show()

    # src_landmarks = sorted_HE_landmarks.astype(np.float32)
    # dst_landmarks = sorted_MxIF_landmarks.astype(np.float32)
    wrt_str = "theta,degrees,s,delta_x,delta_y,dis\n"
    MM = np.zeros((2, 3))
    comb = list(itertools.combinations(range(len(sorted_HE_landmarks)), 2))
    for c in comb:
        src_landmarks = sorted_HE_landmarks[c, :].astype(np.float32) * HE_pixel_size
        dst_landmarks = sorted_MxIF_landmarks[c, :].astype(np.float32) * MxIF_pixel_size

        M = cv2.estimateAffinePartial2D(src_landmarks, dst_landmarks)[0]

        theta = math.atan(M[1, 0] / M[0, 0])
        degrees = theta * (180.0 / math.pi)
        s = M[0, 0] / math.cos(theta)
        delta = math.sqrt(M[0, 2] ** 2 + M[1, 2] ** 2)

        # For debug
        if DEBUG:
            print("From H&E to MxIF (target)")
            print("Rotate: %f degrees %f" % (theta, degrees))
            print("Scale: %f." % s)
            print("Translation: (x=%f, y=%f). Shift distance: %f" % (M[0, 2], M[1, 2], delta))

        wrt_str += ",".join([str(theta), str(degrees), str(s), str(M[0, 2]), str(M[1, 2]), str(delta)])
        wrt_str += "\n"
        # MM = np.vstack((M, np.array([0, 0, 1])))
        MM += M
        new_XY = np.vstack((np.transpose(src_landmarks), np.ones((1, len(src_landmarks)))))
        transformed_loc = np.dot(M, new_XY)
        transformed_landmarks = np.transpose(transformed_loc[0:2, :])

        if DEBUG:
            plt.plot(src_landmarks[:, 0], src_landmarks[:, 1], 'r*')
            plt.plot(dst_landmarks[:, 0], dst_landmarks[:, 1], 'bo')
            plt.plot(transformed_landmarks[:, 0], transformed_landmarks[:, 1], 'g*')
            plt.axis('equal')
            plt.legend(["HE_landmarks", "MxIF_landmarks", "trans_HE_landmarks"])
            plt.show()

    trans_fn = os.path.join(ground_truth_output_dir, roi_id + "_gt_trans.npy")
    data = MM/len(comb)
    np.save(trans_fn, data)

    trans_csv_fn = os.path.join(ground_truth_output_dir, roi_id + "_gt_trans_comb.csv")
    fp = open(trans_csv_fn, 'w')
    fp.write(wrt_str)
    fp.close()
