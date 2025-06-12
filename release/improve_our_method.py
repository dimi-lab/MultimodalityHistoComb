
import numpy as np
import open3d as o3
import transforms3d as t3d
from probreg import cpd
from probreg import callbacks
import os
import pandas as pd
import logging
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from utils import plot_centroids_trans, get_cell_loc
from gm_utils import get_points_in_window, create_graph_for_point_set, create_graph_use_query, GM_match, get_matching_pairs
from gm_utils import draw_graph_with_matching_links, create_delaunay_graph, filter_by_slope
from LPM import *
from sklearn.metrics.pairwise import euclidean_distances
# 


log = logging.getLogger('probreg')
log.setLevel(logging.DEBUG)


# 
data_dir = "/Users/jjiang10/Data/OV_TMA/AlignedCellQuant"
output_dir = "/Users/jjiang10/Data/OV_TMA/output"


morp_features = ["Nucleus: Area µm^2", "Nucleus: Length µm", "Nucleus: Circularity", "Nucleus: Solidity",
                 "Nucleus: Max diameter µm", "Nucleus: Min diameter µm", "Cell: Area µm^2", "Cell: Length µm",
                 "Cell: Circularity", "Cell: Solidity", "Cell: Max diameter µm", "Cell: Min diameter µm",
                 "Nucleus/Cell area ratio"]

HE_H_stain_features = ["Hematoxylin: Nucleus: Mean", "Hematoxylin: Nucleus: Median", "Hematoxylin: Nucleus: Min",
                     "Hematoxylin: Nucleus: Max", "Hematoxylin: Nucleus: Std.Dev."]

MxIF_DAPI_stain_features = ["DAPI_AF_R01: Nucleus: Mean", "DAPI_AF_R01: Nucleus: Median", "DAPI_AF_R01: Nucleus: Min",
                            "DAPI_AF_R01: Nucleus: Max", "DAPI_AF_R01: Nucleus: Std.Dev."]

case_list = ["H-13", "I-15", "I-17"]


for c in case_list:
  HE_fn = os.path.join(data_dir, c+ "_1on1_HE_quant.tsv")
  MxIF_fn = os.path.join(data_dir, c+ "_1on1_MxIF_quant.tsv")

  MxIF_df = pd.read_table(MxIF_fn, sep='\t', header=0)
  HE_df = pd.read_table(HE_fn, sep='\t', header=0)

  HE_centroids_um = get_cell_loc(HE_fn)
  MxIF_centroids_um = get_cell_loc(MxIF_fn)

  # cbs = [callbacks.Open3dVisualizerCallback(source, target, save=True)]
  tf_param, _, _ = cpd.registration_cpd(HE_centroids_um, MxIF_centroids_um,
                                      #   callbacks=cbs,
                                        update_scale=False)

  trans_HE_centroids_um = tf_param.transform(copy.deepcopy(HE_centroids_um))

  title = "Cell Centroids Before and after alignment"
  legend = ["HE cells", "MxIF cells", "trans_HE"]
  fn = f"log_{c}_centroids_alignment.png"
  plot_centroids_trans(HE_centroids_um, MxIF_centroids_um, trans_HE_centroids_um, legend, title, output_dir, fn)

  sample_loc = (max(HE_centroids_um[:, 0]/2), max(HE_centroids_um[:, 1]/2)) # sample the graph from the center of the core
  sample_wind = 150 # sample window size 100 um

  HE_selection, HE_nodes_loc = get_points_in_window(trans_HE_centroids_um, location=sample_loc, window_size=sample_wind)
  HE_nodes_features = HE_df.loc[HE_selection, morp_features+HE_H_stain_features].values
  HE_nodes_features = np.hstack((HE_nodes_features, HE_nodes_loc))
  HE_nodes_features = (HE_nodes_features - HE_nodes_features.min(axis=0)) / (HE_nodes_features.max(axis=0) - HE_nodes_features.min(axis=0))
  # HE_graph, HE_Adj_M = create_graph_for_point_set(HE_nodes_loc)
  HE_graph, HE_Adj_M = create_delaunay_graph(HE_nodes_loc)
  search_range = 1
  q_MxIF_nodes_dis, q_MxIF_idx = create_graph_use_query(HE_nodes_loc, MxIF_centroids_um, K=search_range)

  MxIF_selection = list(set(np.array(q_MxIF_idx).flatten()))
  MxIF_nodes_loc = MxIF_centroids_um[MxIF_selection]
  MxIF_nodes_features = MxIF_df.loc[MxIF_selection, morp_features+MxIF_DAPI_stain_features].values
  MxIF_nodes_features = np.hstack((MxIF_nodes_features, MxIF_nodes_loc))
  MxIF_nodes_features = (MxIF_nodes_features - MxIF_nodes_features.min(axis=0)) / (MxIF_nodes_features.max(axis=0) - MxIF_nodes_features.min(axis=0))

  # MxIF_graph, MxIF_Adj_M = create_graph_for_point_set(MxIF_nodes_loc)
  MxIF_graph, MxIF_Adj_M = create_delaunay_graph(MxIF_nodes_loc)


  # Calculate the similarity matrix using L2 distance (Euclidean distance)
  feature_similarity_matrix = euclidean_distances(HE_nodes_features, MxIF_nodes_features)

  # Calculate the distance matrix between the points in HE_nodes_loc and MxIF_nodes_loc
  distance_matrix = euclidean_distances(HE_nodes_loc, MxIF_nodes_loc)

  # Combine the feature similarity matrix and the distance matrix
  alpha = 0.5  # Weight for combining the matrices
  combined_similarity_matrix = alpha * feature_similarity_matrix + (1 - alpha) * distance_matrix

  # Plot the combined similarity matrix
  # plt.figure()
  # plt.imshow(combined_similarity_matrix)
  # plt.colorbar()
  # plt.title("Combined Similarity Matrix")
  # plt.show()

  # Find the most similar point in MxIF_nodes_loc for each point in HE_nodes_loc
  most_similar_indices = np.argmin(combined_similarity_matrix, axis=1)

  # Plot the two point sets and add links between the most similar points
  # plt.figure()
  # plt.scatter(HE_nodes_loc[:, 0], HE_nodes_loc[:, 1], c='r', label='HE nodes', alpha=0.2)
  # plt.scatter(MxIF_nodes_loc[:, 0], MxIF_nodes_loc[:, 1], c='b', label='MxIF nodes')

  # for i, idx in enumerate(most_similar_indices):
  #     plt.plot([HE_nodes_loc[i, 0], MxIF_nodes_loc[idx, 0]], [HE_nodes_loc[i, 1], MxIF_nodes_loc[idx, 1]], 'k--', alpha=0.5)

  # plt.legend()
  # plt.title("Point Sets with Links Between Most Similar Points")
  # plt.gca().set_aspect('equal', adjustable='box')
  # plt.show()


  # Perform graph matching using the similarity matrix
  X = GM_match(HE_Adj_M, MxIF_Adj_M, HE_nodes_features, MxIF_nodes_features)

  #
  pos1 = {_: HE_nodes_loc[_, :] for _ in range(len(HE_nodes_loc))}
  pos2 = {_: MxIF_nodes_loc[_, :] for _ in range(len(MxIF_nodes_loc))}
  color1 = ['#FF5733' for _ in range(len(HE_nodes_loc))]
  color2 = ['#1f78b4' for _ in range(len(MxIF_nodes_loc))]

  ### get matching pairs
  source_match_idx = range(len(HE_nodes_loc))
  sorted_MxIF_nodes_loc, sorted_MxIF_nodes_idx, target_match_idx = get_matching_pairs(source_match_idx, X, MxIF_nodes_loc, MxIF_selection)
  # ###  filtering the matching pairs
  # s_mask = filter_by_slope(HE_nodes_loc, sorted_MxIF_nodes_loc)
  # src_points = HE_centroids_um[HE_selection][s_mask]
  # dst_points = sorted_MxIF_nodes_loc[s_mask]
  # # s_target_match_idx = target_match_idx[s_mask]
  # s_target_match_idx =[x for x, y in zip(target_match_idx, s_mask) if y]
  # s_source_match_idx =[x for x, y in zip(source_match_idx, s_mask) if y]
  # mask = LPM_filter(src_points, dst_points)
  # mask = mask & s_mask

  mask = LPM_filter(HE_nodes_loc, sorted_MxIF_nodes_loc)

  sv_fn = os.path.join(output_dir, c+ "_gm_before_filtering.png")
  ### Draw the matching links before filtering
  draw_graph_with_matching_links(HE_graph, MxIF_graph,
                                  pos1, pos2, color1, color2,
                                  source_match_idx, target_match_idx, sv_fn, filter=None)

  ## draw the links between matching points after filtering
  sv_fn = os.path.join(output_dir, c+ "_gm_after_filtering.png")
  draw_graph_with_matching_links(HE_graph, MxIF_graph,
                                  pos1, pos2, color1, color2,
                                  source_match_idx, target_match_idx, sv_fn, filter=mask)

print("Done")

