# Multimodality Histochemistry Combination
Multiplexed immunofluorescence (MxIF) signals are often compromised due to the complexities of cyclic staining protocols. We present a novel framework incorporating clinical-grade H&E staining as a reference to enhance the reliability and interpretability of MxIF analyses, focusing on cell spatial co-localization, cell feature validation, and virtual H&E generation.

![Framework](./imgs/framework.png)
Figure 1. Overview of our framework. A) MxIF and H&E images to be aligned; B) Cell segmentations and quantifications; C) Cell-level alignment. 1) coherent point drift. 2) graph matching; D) Aligned MxIF and H&E images with combined cell-level representations.
## Cell spatial co-localization 
Image Alignment
![Alignment](./imgs/alignment.png)

## Cell feature validation


Staining and morphology feature concordance for both re-stained and serial sections.
![Feature validation](./imgs/feature_validation.png)

## Virtual H&E generation

Generating virtual H&E if re-stained image is not available.

![Virtual H&E](./imgs/virtualHE.png)
