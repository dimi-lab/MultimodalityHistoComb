import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
import os
import natsort
import seaborn as sns

data_root_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/AlignmentEval"
aa_root_dir = os.path.join(data_root_dir, "ApplyAlignment")
def get_TMA_core_list(img_path: str) -> list:
    img_fn_list = os.listdir(img_path)
    roi_list = [i.split("_")[0] for i in img_fn_list]
    return list(set(roi_list))

ROI_list = natsort.natsorted(get_TMA_core_list(os.path.join(data_root_dir, "Sec1GroundTruth")))

def read_metric(root_dir, roi_list, seg_method, metric, sec_id, align_method):
    '''
    :param seg_method: "StarDist" or "Watershed"
    :param metric: "RMSE", "theta" or "delta"
    :param sec_id: "Sec1" or "Sec2"
    :param align_method: "GT" or "CPD"
    :return:
    '''
    gt_metric_dict = {"theta": "gt_degrees", "delta": "gt_delta", "RMSE": "rmse_gt"}
    cpd_metric_dict = {"theta": "degrees", "delta": "delta", "RMSE": "rmse_cpd"}
    data_list = []
    for roi in roi_list:
        fn = os.path.join(root_dir, sec_id+"Output", roi + "_" + seg_method + "_metrics.csv")
        df = pd.read_csv(fn, sep=",")
        if align_method == "GT":
            data = df[gt_metric_dict[metric]]
        elif align_method == "CPD":
            data = df[cpd_metric_dict[metric]]
        else:
            raise Exception("metric not defined")
        data_list.append(data[0])
    return data_list




seg_method_list = ["StarDist", "Watershed"]
metric_list = ["RMSE", "theta", "delta"]
units = ["μm", "degree", "μm"]

for idx, metric in enumerate(metric_list):
    gt_sec1 = read_metric(aa_root_dir, ROI_list, seg_method_list[0], metric, "Sec1", "GT")
    cpd_sec1 = read_metric(aa_root_dir, ROI_list, seg_method_list[0], metric, "Sec1", "CPD")
    gt_sec2 = read_metric(aa_root_dir, ROI_list, seg_method_list[0], metric, "Sec2", "GT")
    cpd_sec2 = read_metric(aa_root_dir, ROI_list, seg_method_list[0], metric, "Sec2", "CPD")

    res_sec1 = np.array(gt_sec1) - np.array(cpd_sec1)
    res_sec2 = np.array(gt_sec2) - np.array(cpd_sec2)

    stardist_result_df = pd.DataFrame({
        'Sec1': res_sec1,
        'Sec2': res_sec2
    })

    gt_sec1 = read_metric(aa_root_dir, ROI_list, seg_method_list[1], metric, "Sec1", "GT")
    cpd_sec1 = read_metric(aa_root_dir, ROI_list, seg_method_list[1], metric, "Sec1", "CPD")
    gt_sec2 = read_metric(aa_root_dir, ROI_list, seg_method_list[1], metric, "Sec2", "GT")
    cpd_sec2 = read_metric(aa_root_dir, ROI_list, seg_method_list[1], metric, "Sec2", "CPD")

    res_sec1 = np.array(gt_sec1) - np.array(cpd_sec1)
    res_sec2 = np.array(gt_sec2) - np.array(cpd_sec2)

    watershed_result_df = pd.DataFrame({
        'Sec1': res_sec1,
        'Sec2': res_sec2
    })

    datasets = [stardist_result_df, watershed_result_df]

    # Define which colours you want to use
    colours = ['blue', 'red']
    # Define the groups
    seg_method = ["StarDist", "Watershed"]
    tissue_sec = ['Sec1', 'Sec2']

    # Get the max of the dataset
    all_maximums = [d.max(axis=1).values for d in datasets]
    dataset_maximums = [max(m) for m in all_maximums]
    y_max = max(dataset_maximums)
    # Get the min of the dataset
    all_minimums = [d.min(axis=1).values for d in datasets]
    dataset_minimums = [min(m) for m in all_minimums]
    y_min = min(dataset_minimums)
    # Calculate the y-axis range
    y_range = y_max - y_min

    # Set x-positions for boxes
    x_pos_range = np.arange(len(datasets)) / (len(datasets) - 1)
    x_pos = (x_pos_range * 0.25) + 0.75
    # Plot
    for i, data in enumerate(datasets):
        positions = [x_pos[i] + j * 1 for j in range(len(data.T))]
        bp = plt.boxplot(
            np.array(data), sym='', whis=[0, 100], widths=0.3 / len(datasets),
            labels=list(datasets[0]), patch_artist=True,
            positions=positions
        )
        # Fill the boxes with colours (requires patch_artist=True)
        k = i % len(colours)
        for box in bp['boxes']:
            box.set(facecolor=colours[k])
        # Make the median lines more visible
        plt.setp(bp['medians'], color='black')

        # Get the samples' medians
        medians = [bp['medians'][j].get_ydata()[0] for j in range(len(data.T))]
        medians = [str(round(s, 2)) for s in medians]
        # Increase the height of the plot by 5% to fit the labels
        plt.ylim([y_min - 0.1 * y_range, y_max + 0.05 * y_range])
        # Set the y-positions for the labels
        y_pos = y_min - 0.075 * y_range
        # for tick, label in zip(range(len(data.T)), plt.xticks()):
        #     k = tick % 2
        #     # plt.text(
        #     #     positions[tick], y_pos, r'$\tilde{x} =' + fr' {medians[tick]}$m',
        #     #     horizontalalignment='center', size='medium'
        #     # )

    # Titles
    plt.title('Alignment Error: %s' % metric)
    plt.ylabel(units[idx])
    # Axis ticks and labels
    plt.xticks(np.arange(len(list(datasets[0]))) + 0.85)
    plt.gca().xaxis.set_minor_locator(ticker.FixedLocator(
        np.array(range(len(list(datasets[0])) + 1)) + 0.5)
    )
    plt.gca().tick_params(axis='x', which='minor', length=4)
    plt.gca().tick_params(axis='x', which='major', length=0)
    # Change the limits of the x-axis
    plt.xlim([0.5, len(list(datasets[0])) + 0.25])
    # Legend
    legend_elements = []
    for i in range(len(datasets)):
        j = i % len(seg_method)
        k = i % len(colours)
        legend_elements.append(Patch(facecolor=colours[k], label=seg_method[j]))
    plt.legend(handles=legend_elements, fontsize=8)
    plt.tight_layout()
    plt.show()

    print("Done")


