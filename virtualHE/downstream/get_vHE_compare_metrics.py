import numpy as np
from PIL import Image
from natsort import natsorted
from sklearn.naive_bayes import GaussianNB
import sys,os
import metrikz
import multiprocessing

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

def cal_eval(eval_files):
    input_fn, output_fn, target_fn = eval_files
    input_patch = Image.open(input_fn)
    restored_patch = Image.open(output_fn)
    clean_patch = Image.open(target_fn)
    PSNR, SSIM, VIF = calculate_metrics(np.array(clean_patch), np.array(input_patch))
    PSNR_r, SSIM_r, VIF_r = calculate_metrics(np.array(clean_patch), np.array(restored_patch))
    # TMA_ID, loc_x, loc_y, real_fake, AB = parse_all_fn(input_fn)
    # wrt_str = marked_wsi_uuid + "," + str(marked_loc_x) + "," + str(marked_loc_y) + "," + clean_wsi_uuid + "," + str(clean_loc_x) + "," + str(clean_loc_y)
    # wrt_str += "," + str(PSNR) + "," + str(PSNR_r) + "," + str(SSIM) + "," + str(SSIM_r) + "," + str(VIF) + "," + str(VIF_r) + "\n"
    wrt_str = input_fn + "," + str(PSNR) + "," + str(PSNR_r) + "," + str(SSIM) + "," + str(SSIM_r) + "," + str(VIF) + "," + str(VIF_r) + "\n"
    metric_arr = [PSNR, PSNR_r, SSIM, SSIM_r, VIF, VIF_r]
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
    metrics_csv = "/data_folder/Ovarian_TMA/virtual_HE/pytorch_model_testing/vHE_pix2pix/metrics/psnr_ssim_vif.csv"
    bufsize = 1  # 0 means unbuffered, 1 means line buffered,
    metrics_csv_fp = open(metrics_csv, 'a', buffering=bufsize)
    wrt_str = "image,PSNR,PSNR_r,SSIM,SSIM_r,VIF,VIF_r\n"
    metrics_csv_fp.write(wrt_str)

    val_file_names = get_all_val_fns(data_dir)
    metric_arr_stack = np.empty([1, 6])

    # for eval_files in val_file_names:
    #     wrt_str, metric_arr = cal_eval(eval_files)

    multiprocessing.set_start_method('spawn')
    pool = multiprocessing.Pool(processes=5, initializer=init_batch_eval)
    for wrt_str, metric_arr in pool.imap(cal_eval, val_file_names):
        metrics_csv_fp.write(wrt_str)
        metric_arr_stack = np.vstack([metric_arr_stack, metric_arr])

    metrics_csv_fp.close()
    print("MxIF vs virtual H&E, mean_PSNR, mean_SSIM, mean_VIF")
    print(np.mean(metric_arr_stack, axis=0))







