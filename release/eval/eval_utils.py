import os
import numpy as np
from numpy import dot
from numpy.linalg import norm

def read_landmarks(fn):
    with open(fn) as f:
        lines = (line for line in f if not line.startswith('x'))
        HE_points = np.loadtxt(lines, delimiter=',', skiprows=0)
    return HE_points


def get_sorted_annotation_landmark_pairs(HE_points_fn, MxIF_points_fn, sort_fn):
    HE_points = read_landmarks(HE_points_fn)
    MxIF_points = read_landmarks(MxIF_points_fn)
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


def apply_aff_trans2points(source, M):
    '''
    source: cell centroids from H&E segmentation (N*2), N is the number of cell
    '''
    new_XY = np.vstack((np.transpose(source), np.ones((1, len(source)))))
    transformed_loc = np.dot(M, new_XY)
    affined_points = np.transpose(transformed_loc[0:2, :])
    return affined_points

def calculate_transformed_landmark_dist(trans_landmarks, target_landmarks):
    dist_arr = trans_landmarks - target_landmarks
    a = np.sum(np.multiply(dist_arr, dist_arr), axis=1)
    dist = np.average(np.sqrt(a))
    return dist

