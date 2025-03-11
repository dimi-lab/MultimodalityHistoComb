import os, io
import glob
import numpy as np
import multiprocessing
import tifffile as tf
import matplotlib.pyplot as plt
from natsort import natsorted
from PIL import Image
from skimage.color import rgb2lab
import xml.etree.ElementTree

# Define some run parameters
num_processors = 20  # Number of processes that can be running at once
he_wsi_fn_dir = "/data_folder/Ovarian_TMA/AlignmentEval/ApplyAlignment/Sec1Output"
he_wsi_fn_list = natsorted(glob.glob(os.path.join(he_wsi_fn_dir, "*_gt_affined.tif")))
print(len(he_wsi_fn_list))

mxif_wsi_fn_dir = "/data_folder/Ovarian_TMA/MxIF_FOVs/Slide2050_24Plex"


output_dir = "/data_folder/Ovarian_TMA/virtual_HE/dataset/ImgPair"  # Define an output directory

def get_tiff_channel_names(fn):
    ome_tiff = tf.TiffFile(fn)
    # omexml_string = ome_tiff.pages[0].description
    # root = xml.etree.ElementTree.parse(io.StringIO(omexml_string))
    # namespaces = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
    # channels = root.findall('ome:Image[1]/ome:Pixels/ome:Channel', namespaces)
    # channel_names = [c.attrib['Name'] for c in channels]
    channel_names = ome_tiff.imagej_metadata.get("Labels")

    return channel_names

def get_channel_idx(channel_name_list, channel_query_list):
    channel_index = []
    for c in channel_query_list:
        for idx, i in enumerate(channel_name_list):
            if i.lower() == c.lower():
                channel_index.append(idx)
    return channel_index

def tissue_detection(rgb_img_arr, threshold=85):
    lab_img = rgb2lab(rgb_img_arr)
    l_img = lab_img[:, :, 0]
    # tissue is darker than background, recommend threshold value: 85
    binary_img_array_1 = np.array(0 < l_img)
    binary_img_array_2 = np.array(l_img < threshold)
    binary_img_array = np.logical_and(binary_img_array_1, binary_img_array_2) 
    return binary_img_array

def get_case_id_from_fn(abs_fn):
    fn = os.path.split(abs_fn)[1]
    return fn.replace("_gt_affined.tif", "")

def get_HE_patch_locations(img_arr, stride=128, patch_size=256):
    h, w, c = img_arr.shape
    assert h == w
    valid_start_loc_list = []
    for loc in range(0, h, stride):
        if loc + patch_size < h:
            valid_start_loc_list.append(loc)
    return valid_start_loc_list

def read_img_arr(img_fn):
    ext = os.path.splitext(img_fn)[1]
    if ext == ".tif" or ext == ".tiff":
        img_arr = tf.imread(img_fn)
        return img_arr
    else:
        print("Currently, doesn't support this format")
        return None

# %
if __name__ == '__main__':
    selected_channels = ["DAPI_AF_R01", "NAK", "PANCK"]
    patch_sz = 256
    pixel_cnt = 256*256
    for he_wsi_fn in he_wsi_fn_list:
        print(he_wsi_fn)
        case_id = get_case_id_from_fn(he_wsi_fn)
        mxif_wsi_fn = os.path.join(mxif_wsi_fn_dir, case_id+".tif")
        print(mxif_wsi_fn)
        he_img_arr = read_img_arr(he_wsi_fn)
        he_img_arr[np.all(he_img_arr==0, axis=-1)] = [240, 242, 240] # replace black background with gray
        
        channel_names = get_tiff_channel_names(mxif_wsi_fn)
        selected_channel_idx = get_channel_idx(channel_names, selected_channels)
        assert len(selected_channel_idx) == 3
        mxif_img_arr = read_img_arr(mxif_wsi_fn)
        mxif_img_arr = np.swapaxes(mxif_img_arr, 0, 2)
        mxif_img_arr = np.swapaxes(mxif_img_arr, 0, 1)[:, :, selected_channel_idx]
        # mxif_img_arr[np.all(mxif_img_arr == 0, axis=-1)] = [255, 255, 255]
        mxif_img_arr[:, :, 1] = mxif_img_arr[:, :, 0] + mxif_img_arr[:, :, 1]
        mxif_img_arr[:, :, 2] = mxif_img_arr[:, :, 0] + mxif_img_arr[:, :, 2]
        start_loc_list = get_HE_patch_locations(he_img_arr, stride=128, patch_size=256)

        if not os.path.exists(os.path.join(output_dir, case_id)):
            os.makedirs(os.path.join(output_dir, case_id))

        for idx_r, s_loc_r in enumerate(start_loc_list):
            for idx_c, s_loc_c in enumerate(start_loc_list):
                he_img_patch_arr = he_img_arr[s_loc_r: s_loc_r+patch_sz, s_loc_c:s_loc_c+patch_sz]
                foreground_pixs = np.sum(tissue_detection(he_img_patch_arr))
                if foreground_pixs/pixel_cnt < 0.3:
                    continue
                mxif_img_patch_arr = mxif_img_arr[s_loc_r: s_loc_r+patch_sz, s_loc_c:s_loc_c+patch_sz]
                # plt.imshow(he_img_patch_arr)
                # plt.savefig("he_img.png")
                # plt.imshow(mxif_img_patch_arr)
                # plt.savefig("mxif_img.png")

                save_img_fn = os.path.join(output_dir, case_id, "_".join([case_id,str(s_loc_r),str(s_loc_c)])+".jpg")
                sv_img_arr = np.zeros((patch_sz, patch_sz*2, 3), dtype=np.uint8)
                sv_img_arr[:, 0:patch_sz] = mxif_img_patch_arr
                sv_img_arr[:, patch_sz:patch_sz*2] = he_img_patch_arr
                Image.fromarray(sv_img_arr).save(save_img_fn)

                print(save_img_fn)

    # # Run the extraction process
    # multiprocessing.set_start_method('spawn')
    # pool = multiprocessing.Pool(processes=num_processors)
    # pool.map(patch_extractor.extract, [wsi_fn_list])
