import math

import numpy as np
from sklearn.metrics import pairwise_distances
import networkx as nx
import matplotlib.pyplot as plt
import pygmtools as pygm
from matplotlib.patches import ConnectionPatch
from scipy.spatial import KDTree
import functools
from csbdeep.utils import Path, normalize


def get_points_in_window(all_points, location=(500, 600), window_size=100):
    x_idx2 = np.logical_and(all_points[:, 0] > location[0], all_points[:, 0] < location[0] + window_size)
    y_idx2 = np.logical_and(all_points[:, 1] > location[1], all_points[:, 1] < location[1] + window_size)
    selected_idx = np.logical_and(x_idx2, y_idx2)
    selected_points = all_points[selected_idx, :]
    return selected_idx, selected_points


def create_graph_for_point_set(points, dis=15):
    dm1 = pairwise_distances(np.array(points))
    Adj_Matrix = dm1 * (dm1 < dis)
    np.fill_diagonal(Adj_Matrix, 0)
    G1 = nx.from_numpy_array(Adj_Matrix)
    return G1, Adj_Matrix


def draw_two_graphs(G1, G2, pos1, pos2, color1, color2):
    plt.figure(figsize=(8, 4))
    ax1 = plt.subplot(1, 2, 1)
    plt.title('Source')
    plt.gca().margins(0.4)
    nx.draw_networkx(G1, pos=pos1, node_color=color1, with_labels=True, node_size=5)
    plt.subplot(1, 2, 2, sharex=ax1, sharey=ax1)
    plt.title('Target')
    nx.draw_networkx(G2, pos=pos2, node_color=color2, with_labels=True, node_size=5)
    plt.tight_layout()
    plt.show()


def create_graph_use_query(source_points, target_points, K=5):
    target_tree = KDTree(target_points, leafsize=5)
    # q_targets, q_loc = target_tree.query(source_points, k=K, distance_upper_bound=30)
    q_targets, q_loc = target_tree.query(source_points, k=K)

    return q_targets, q_loc


def GM_match(source_Adj_M, target_Adj_M, source_node_features, target_node_features):
    source_conn, source_edge = pygm.utils.dense_to_sparse(source_Adj_M)
    target_conn, target_edge = pygm.utils.dense_to_sparse(target_Adj_M)
    aff = functools.partial(pygm.utils.gaussian_aff_fn, sigma=.0001)
    # aff = functools.partial(pygm.utils.inner_prod_aff_fn)

    n_source = np.array([len(source_node_features)])
    n_target = np.array([len(target_node_features)])
    norm_source_node_features = normalize(source_node_features)
    norm_target_node_features = normalize(target_node_features)
    K = pygm.utils.build_aff_mat(norm_source_node_features, source_edge, source_conn,
                                 norm_target_node_features, target_edge, target_conn,
                                 n_source, None,
                                 n_target, None,
                                 edge_aff_fn=aff)

    np.fill_diagonal(K, 0)
    # X = pygm.rrwm(K, n_source, n_target, max_iter=5, alpha=0, beta=30)
    X = pygm.ipfp(K, n_source, n_target)
    X = pygm.sinkhorn(X)
    # X = pygm.hungarian(X)
    return X


def get_matching_pairs(source_match_idx, X, target_nodes_loc, MxIF_selection):
    match_idx_target = []
    for i in source_match_idx:
        j = np.argmax(X[i]).item()
        match_idx_target.append(j)
    sorted_MxIF_nodes_loc = []
    sorted_MxIF_nodes_idx = []
    for idx, sn in enumerate(match_idx_target):
        sorted_MxIF_nodes_loc.append(target_nodes_loc[sn, :])
        sorted_MxIF_nodes_idx.append(MxIF_selection[sn])
    return np.array(sorted_MxIF_nodes_loc), np.array(sorted_MxIF_nodes_idx), match_idx_target


def draw_graph_with_matching_links(source_g, target_g,
                                   source_pos, target_pos, source_color, target_color,
                                   source_match_idx, target_match_idx, filter=None):
    plt.figure(figsize=(8, 4))
    ax1 = plt.subplot(1, 2, 1)
    plt.title('Source')
    plt.gca().margins(0.4)
    nx.draw_networkx(source_g, pos=source_pos, node_color=source_color, with_labels=True, node_size=5)
    ax2 = plt.subplot(1, 2, 2, sharex=ax1, sharey=ax1)
    plt.title('Target')
    nx.draw_networkx(target_g, pos=target_pos, node_color=target_color, with_labels=True, node_size=5)

    for i in source_match_idx:
        if filter is None:
            j = target_match_idx[i]
            print([i, j, source_pos[i], target_pos[j]])
            con = ConnectionPatch(xyA=source_pos[i], xyB=target_pos[j], coordsA="data", coordsB="data",
                                  axesA=ax1, axesB=ax2, color="green")
            plt.gca().add_artist(con)
        else:
            if filter[i]:
                j = target_match_idx[i]
                print([i, j, source_pos[i], target_pos[j]])
                con = ConnectionPatch(xyA=source_pos[i], xyB=target_pos[j], coordsA="data", coordsB="data",
                                      axesA=ax1, axesB=ax2, color="green")
                plt.gca().add_artist(con)
    plt.show()

def check_CPD_GM(CPD_M, GM_M, rot_threshold=1, dis_threshold=20):
    CPD_theta = math.asin(CPD_M[0, 1]) * (180.0 / math.pi)
    CPD_delta = math.sqrt(CPD_M[0, 2]**2 + CPD_M[1, 2]**2)
    GM_theta = math.asin(GM_M[0, 1]) * (180.0 / math.pi)
    GM_delta = math.sqrt(GM_M[0, 2]**2 + GM_M[1, 2]**2)
    if abs(CPD_theta-GM_theta) < rot_threshold and abs(CPD_delta-GM_delta) < dis_threshold:
        return True
    else:
        return False
def get_M_from_cv2_affine(M, source_pix_size, target_pix_size):
    real_S = source_pix_size / target_pix_size
    M[0, 0] = M[0, 0] * real_S
    M[1, 1] = M[1, 1] * real_S
    a = M[0, 2] / target_pix_size
    b = M[1, 2] / target_pix_size  # please note, the x,y coordinate were swapped to w,h
    M[0, 2] = b
    M[1, 2] = a
    return M
