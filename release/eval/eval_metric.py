import os
from eval_utils import *
import matplotlib.pyplot as plt


Sec = 2

HE_pixel_size = 0.2201  # unit micron
MxIF_pixel_size = 0.325

core_list_fn = "../core_list.txt"
core_list = open(core_list_fn, 'r').readlines()
roi_id_list = []
for core in core_list:
    roi_id = core.split(".")[0]
    roi_id_list.append(roi_id)

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

MxIF_export_dir = os.path.join(data_root_dir, "GroundTruthAnnotation", "MxIF_Sec1")

gt_data_root = "/temp/Ovarian_TMA/AlignmentEval/GroundTruthEvaluation"

qupath_results_root = "/temp/Ovarian_TMA/test_align_seg_qupath"
results_root = "/temp/Ovarian_TMA/test_align_wsi"

def draw_points(source_points, target_points, title):
    plt.figure(1)
    plt.scatter(source_points[:, 1], source_points[:, 0], c='r', s=10)
    plt.scatter(target_points[:, 1], target_points[:, 0], c='b', s=10)
    plt.legend(["HE", "MxIF"])
    plt.title(title)
    plt.axis('equal')
    plt.show()

for roi in roi_id_list:
    print("Processing %s" % roi)
    HE_landmarks_fn = os.path.join(HE_export_dir, roi+"_align_anno_8.csv")
    MxIF_landmarks_fn = os.path.join(MxIF_export_dir, roi+"_align_anno_8.csv")
    point_sort_fn = os.path.join(data_root_dir, "GroundTruthAnnotation", "Landmark_correspondence", "Sec%d" %Sec, roi+".csv")

    gt_trans_fn = os.path.join(gt_data_root, "GT_HE_Sec%d_MxIF" % Sec, roi+"_gt_trans.npy")
    final_trans_fn = os.path.join(results_root, "Sec%d" % Sec, roi+".tif_final_M.npy")
    cpd_trans_fn = os.path.join(results_root, "Sec%d" % Sec, roi + ".tif_CPD_M.npy")

    sorted_HE_landmarks, sorted_MxIF_landmarks = get_sorted_annotation_landmark_pairs(HE_landmarks_fn, MxIF_landmarks_fn, point_sort_fn)
    sorted_HE_landmarks_um = sorted_HE_landmarks * HE_pixel_size
    sorted_MxIF_landmarks_um = sorted_MxIF_landmarks * MxIF_pixel_size

    gt_M = np.load(gt_trans_fn)
    final_M = np.load(final_trans_fn)
    CPD_M = np.load(cpd_trans_fn)

    # check if CPD_M is the same as final_M
    a = final_M - CPD_M
    if sum(a.flatten()) != 0:
        print("\t %s CPD results were fine tuned" % roi)

    # check distance between transformed source landmark points and target landmark points
    gt_trans_HE_landmarks_um = apply_aff_trans2points(sorted_HE_landmarks_um, gt_M)
    gt_dist = calculate_transformed_landmark_dist(gt_trans_HE_landmarks_um, sorted_MxIF_landmarks_um)

    trans_HE_landmarks_um = apply_aff_trans2points(sorted_HE_landmarks, final_M)*MxIF_pixel_size
    final_dist = calculate_transformed_landmark_dist(trans_HE_landmarks_um, sorted_MxIF_landmarks_um)

    cpd_trans_HE_landmarks_um = apply_aff_trans2points(sorted_HE_landmarks, CPD_M)*MxIF_pixel_size
    cpd_dist = calculate_transformed_landmark_dist(cpd_trans_HE_landmarks_um, sorted_MxIF_landmarks_um)

    # debug
    if sum(a.flatten()) != 0:
    # if roi == "C-10":
        draw_points(gt_trans_HE_landmarks_um, sorted_MxIF_landmarks_um, "Landmarks with ground truth transformation")
        draw_points(cpd_trans_HE_landmarks_um, sorted_MxIF_landmarks_um, "Landmarks with CPD transformation")
        draw_points(trans_HE_landmarks_um, sorted_MxIF_landmarks_um, "Landmarks with final transformation")

    print("\t Distance differences: %f\t %f\t %f" % (gt_dist, cpd_dist, final_dist))

    # check theta and delta differences











