import os
from skimage.morphology import binary_dilation, disk, binary_erosion
import numpy as np

'''
steps:
1. create synthetic volume image
2. run matlab code: https://github.com/smujiang/my_matlab_lib/blob/master/CONVERT_voxels_to_stl/CONVERT_voxels_to_stl.m
3. load the generated STL file to viewer.
'''

def create_synthetic_volume_img(cell_mask, radius=1, layers=5):
    '''
    create synthetic volume image based on cell segmentation, each time for a cell
    :param cell_mask:
    :param layers:
    :return:
    '''
    img_volume_arr = np.zeros((cell_mask.shape[0], cell_mask.shape[1], layers*2 + 1), dtype=np.uint8)
    img_volume_arr[:, :, layers] = cell_mask
    for i in range(layers):
        tmp_mask = binary_erosion(img_volume_arr[:, :, layers - i], disk(radius, dtype=bool))
        img_volume_arr[:, :, layers - i-1] = tmp_mask
    for i in range(layers):
        tmp_mask = binary_erosion(img_volume_arr[:, :, layers + i], disk(radius, dtype=bool))
        img_volume_arr[:, :, layers + i+1] = tmp_mask
    return img_volume_arr


if __name__ == "__main__":
    from skimage.data import human_mitosis
    import matplotlib.pyplot as plt
    from skimage import (
        color, feature, filters, measure, morphology, segmentation, util
    )

    image = human_mitosis()

    fig, ax = plt.subplots()
    ax.imshow(image, cmap='gray')
    ax.set_title('Microscopy image of human cells stained for nuclear DNA')
    plt.show()

    higher_threshold = 125
    dividing = image > higher_threshold

    smoother_dividing = filters.rank.mean(util.img_as_ubyte(dividing),
                                          morphology.disk(4))

    binary_smoother_dividing = smoother_dividing > 10

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(binary_smoother_dividing)
    ax.set_title('Dividing nuclei')
    ax.axis('off')
    plt.show()

    cleaned_dividing = measure.label(binary_smoother_dividing)
    print(cleaned_dividing.max())

    for i in range(1, cleaned_dividing.max()):
        cell_mask = np.zeros(image.shape)
        cell_mask[cleaned_dividing==i] = 1
        plt.imshow(cell_mask)
        plt.show()

        cell_volume_arr = create_synthetic_volume_img(cell_mask, radius=1, layers=5)
        np.save("./temp.npy", cell_volume_arr)
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.imshow(cell_volume_arr[:, :, 3])
        ax.set_title('Dividing nuclei')
        ax.axis('off')
        plt.show()


        fig, ax_list = plt.subplots(2, 5)
        for idx, ax in enumerate(ax_list.flatten()):
            ax.imshow(cell_volume_arr[:, :, idx], 'gray')
            ax.axis('off')
        plt.tight_layout()
        plt.show()














