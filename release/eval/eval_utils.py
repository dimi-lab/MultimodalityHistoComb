import os
import numpy as np

def get_sorted_annotation_landmark_pairs(HE_points_fn, MxIF_points_fn, sort_fn):
    with open(HE_points_fn) as f:
        lines = (line for line in f if not line.startswith('x'))
        HE_points = np.loadtxt(lines, delimiter=',', skiprows=0)
    with open(MxIF_points_fn) as f:
        lines = (line for line in f if not line.startswith('x'))
        MxIF_points = np.loadtxt(lines, delimiter=',', skiprows=0)
    assert len(HE_points) == len(MxIF_points)

    HE_points_tmp = HE_points[HE_points[:, 0].argsort()]
    MxIF_points_tmp = MxIF_points[MxIF_points[:, 0].argsort()]

    with open(sort_fn) as f:
        lines = f.readlines()
        orders = np.loadtxt(lines, delimiter=',').astype(int)
    sorted_HE_landmarks = np.zeros(HE_points_tmp.shape)
    sorted_MxIF_landmarks = np.zeros(MxIF_points_tmp.shape)
    for i, o in enumerate(orders):
        sorted_HE_landmarks[i, :] = HE_points_tmp[o[0], :]
        sorted_MxIF_landmarks[i, :] = MxIF_points_tmp[o[1], :]
    return sorted_HE_landmarks, sorted_MxIF_landmarks