import numpy as np
from PIL import Image
from natsort import natsorted
from sklearn.naive_bayes import GaussianNB
import sys,os
import metrikz
import multiprocessing
import matplotlib.pyplot as plt

def calculate_metrics(Img_a, Img_b):
    PSNR = metrikz.psnr(Img_a, Img_b)
    SSIM = metrikz.ssim(Img_a, Img_b)
    VIF = metrikz.pbvif(Img_a, Img_b)
    # print("A: PSNR: %f, SSIM: %f, VIF: %f" % (PSNR, SSIM, VIF))
    return PSNR, SSIM, VIF
    # PSNR = psnr(Img_a, Img_b)
    # SSIM = ssim(Img_a, Img_b)
    # VIF = vifp(Img_a, Img_b)
    # print("B: PSNR: %f, SSIM: %f, VIF: %f" % (PSNR, SSIM[0], VIF))
    # return PSNR, SSIM[0], VIF

def char_to_index(char):
    """
    Converts an alphabet character to its corresponding number index (0-based).

    Args:
        char (str): The alphabet character (case-insensitive).

    Returns:
        int: The number index corresponding to the character (0-based).
    """
    char = char.lower()  # Convert to lowercase
    if 'a' <= char <= 'z':
        return ord(char) - ord('a')
    else:
        raise ValueError(f"Invalid input: {char} is not an alphabet character.")

def read_CellDive_vHE(vHE_dir, roi, loc_x, loc_y, patch_size=256):
    '''

    :param vHE_dir:
    :param roi:
    :param loc_x:
    :param loc_y:
    :param patch_size:
    :param core_per_row: number of cores per row
    :return:
    '''
    core_vHE_fn = os.path.join(vHE_dir, roi + ".tif")
    img_core_arr = np.array(Image.open(core_vHE_fn))
    img_arr = img_core_arr[loc_x:loc_x+patch_size, loc_y:loc_y+patch_size, :]
    # plt.imshow(img_core_arr)
    # plt.show()
    return img_arr

vHE_dir = "/data_folder/Ovarian_TMA/CellDive_virtual_HE"
def cal_eval(eval_files):
    mxif_fn, my_vHE_fn, real_HE_fn = eval_files
    mxif_patch = Image.open(mxif_fn)
    my_vHE_patch = Image.open(my_vHE_fn)
    real_HE_patch = Image.open(real_HE_fn)
    ele = os.path.split(mxif_fn)[1].split("_")
    roi, loc_x, loc_y = [ele[0], int(ele[1]), int(ele[2])]
    patch_size = real_HE_patch.width
    CellDive_vHE_patch = read_CellDive_vHE(vHE_dir, roi, loc_x, loc_y, patch_size=patch_size)

    PSNR, SSIM, VIF = calculate_metrics(np.array(real_HE_patch), np.array(mxif_patch))
    PSNR_r, SSIM_r, VIF_r = calculate_metrics(np.array(real_HE_patch), np.array(my_vHE_patch))
    PSNR_rr, SSIM_rr, VIF_rr = calculate_metrics(np.array(real_HE_patch), np.array(CellDive_vHE_patch))

    # fig, axes = plt.subplots(1, 4, dpi=250)
    # axes[0].imshow(mxif_patch)
    # axes[0].set_title('MxIF')
    # axes[0].set_xlabel("PSNR:%.2f\nSSIM:%.2f\nVIF:%.2f" % (PSNR, SSIM, VIF), fontsize=10)
    # axes[0].set_xticks([])  # Hide the axis ticks and labels
    # axes[0].set_yticks([])
    # axes[1].imshow(real_HE_patch)
    # axes[1].set_title('Real H&E')
    # axes[1].set_xlabel("PSNR:%.2f\nSSIM:%.2f\nVIF:%.2f" % (48, 1.0, 1.0), fontsize=10)
    # axes[1].set_xticks([])  # Hide the axis ticks and labels
    # axes[1].set_yticks([])
    # axes[2].imshow(my_vHE_patch)
    # axes[2].set_title('Our vH&E')
    # axes[2].set_xlabel("PSNR:%.2f\nSSIM:%.2f\nVIF:%.2f" % (PSNR_r, SSIM_r, VIF_r), fontsize=10)
    # axes[2].set_xticks([])  # Hide the axis ticks and labels
    # axes[2].set_yticks([])
    # axes[3].imshow(CellDive_vHE_patch)
    # axes[3].set_title('Vendor vH&E')
    # axes[3].set_xlabel("PSNR:%.2f\nSSIM:%.2f\nVIF:%.2f" % (PSNR_rr, SSIM_rr, VIF_rr), fontsize=10)
    # axes[3].set_xticks([])  # Hide the axis ticks and labels
    # axes[3].set_yticks([])
    # plt.tight_layout()
    # plt.suptitle("_".join([roi, str(loc_x), str(loc_y)]))
    # plt.show()
    #

    wrt_str = mxif_fn + "," + str(PSNR) + "," + str(PSNR_r) + "," + str(PSNR_rr) + "," \
              + str(SSIM) + "," + str(SSIM_r) + "," + str(SSIM_rr) + "," \
              + str(VIF) + "," + str(VIF_r) + "," + str(VIF_rr) + "\n"
    metric_arr = [PSNR, PSNR_r, PSNR_rr, SSIM, SSIM_r, SSIM_rr, VIF, VIF_r, VIF_rr]
    print(wrt_str)
    print(metric_arr)
    return wrt_str, metric_arr



def parse_all_fn(img_file_name):
    ele = img_file_name.split("_", maxsplit=5)
    TMA_ID, loc_x, loc_y, real_fake, AB_png = ele
    AB = AB_png.split(".")[0]
    return TMA_ID, loc_x, loc_y, real_fake, AB


def get_pairwise_fn(input_path_name):
    path_ele = os.path.split(input_path_name)
    ele = path_ele[1].rsplit("_", maxsplit=2)
    pre_fix, real_fake, others = ele
    if real_fake == "fake":
        target_fn = os.path.join(path_ele[0], pre_fix + "_real_B." + others.split(".")[1])
        input_fn = os.path.join(path_ele[0], pre_fix + "_real_A." + others.split(".")[1])
        if os.path.exists(target_fn) and os.path.exists(input_fn):
            return target_fn, input_fn
        else:
            print(target_fn)
            print(input_fn)
            raise Exception("Target or output file missing")
    else:
        raise Exception("Input parameter should be inputs")


def get_all_val_fns(data_dir):
    val_file_names = []
    img_fns = os.listdir(data_dir)
    print("There are %d samples to be calculated" % (len(img_fns) / 3))
    for fn in img_fns:
        TMA_ID, loc_x, loc_y, real_fake, AB = parse_all_fn(fn)
        if "fake" == real_fake:
            output_fn = os.path.join(data_dir, fn)
            target_fn, input_fn = get_pairwise_fn(output_fn)
            val_file_names.append([input_fn, output_fn, target_fn])
    return val_file_names


def init_batch_eval():
    global cnt
    cnt = 0

if __name__ == "__main__":
    data_dir = "/data_folder/Ovarian_TMA/virtual_HE/pytorch_model_testing/vHE_pix2pix/test_latest/images"
    metrics_csv = "/data_folder/Ovarian_TMA/virtual_HE/pytorch_model_testing/vHE_pix2pix/metrics/psnr_ssim_vif_all.csv"
    bufsize = 1  # 0 means unbuffered, 1 means line buffered,
    metrics_csv_fp = open(metrics_csv, 'a', buffering=bufsize)
    wrt_str = "image,PSNR,PSNR_r,PSNR_rr,SSIM,SSIM_r,SSIM_rr,VIF,VIF_r,VIF_rr\n"
    metrics_csv_fp.write(wrt_str)

    val_file_names = get_all_val_fns(data_dir)
    metric_arr_stack = np.empty([1, 9])

    # for eval_files in val_file_names:
    #     wrt_str, metric_arr = cal_eval(eval_files)
    #     metrics_csv_fp.write(wrt_str)

    multiprocessing.set_start_method('spawn')
    pool = multiprocessing.Pool(processes=5, initializer=init_batch_eval)
    for wrt_str, metric_arr in pool.imap(cal_eval, val_file_names):
        metrics_csv_fp.write(wrt_str)
        metric_arr_stack = np.vstack([metric_arr_stack, metric_arr])

    metrics_csv_fp.close()
    print("MxIF vs virtual H&E, mean_PSNR, mean_SSIM, mean_VIF")
    print(np.mean(metric_arr_stack, axis=0))







