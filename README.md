# Multimodality Histochemistry Combination
Multiplexed immunofluorescence (MxIF) signals are often compromised due to the complexities of cyclic staining protocols. We present a novel framework incorporating clinical-grade H&E staining as a reference to enhance the reliability and interpretability of MxIF analyses, focusing on cell spatial co-localization, cell feature validation, and virtual H&E generation.




![Framework](./imgs/framework.png)
Figure 1. Overview of our framework. A) MxIF and H&E images to be aligned; B) Cell segmentations and quantifications; C) Cell-level alignment. 1) coherent point drift. 2) graph matching; D) Aligned MxIF and H&E images with combined cell-level representations.

## Installation
```Shell
conda create --prefix /path_to_your_envs/alignment python=3.8.2
conda activate alignment
pip install scikit-learn
pip install probreg
pip install pygmtools
pip install cupy
```

## Modules
### Preprocessing
#### TMA de-array
De-array is a process that gets the tissue core out of the whole slide image, and then save the cores into FOV images.
#### Cell segmentation
StarDist was used as our baseline for cell segmentation. However, our framework is portable to different cell segmentation methods.
### Cell spatial co-localization 
Image Alignment
![Alignment](./imgs/alignment.png)

### Cell feature validation


Staining and morphology feature concordance for both re-stained and serial sections.
![Feature validation](./imgs/feature_validation.png)

### Virtual H&E generation

Generating virtual H&E if re-stained image is not available.

![Virtual H&E](./imgs/virtualHE.png)


## References
Preprint: [Jiang, Jun, Raymond Moore, Brenna Novotny, Leo Liu, Zachary Fogarty, Ray Guo, Markovic Svetomir, and Chen Wang. "Multimodal Alignment of Histopathological Images Using Cell Segmentation and Point Set Matching for Integrative Cancer Analysis." arXiv preprint arXiv:2410.00152 (2024).](https://arxiv.org/abs/2410.00152)
