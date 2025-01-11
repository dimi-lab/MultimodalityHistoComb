To run our alignment, cells within MxIF and H&E images need to be detected respectively. We used StarDist as our baseline for both H&E and MxIF, as there are pretrained model available. 

Cell segmentation methods we tried.
## H&E cell segmentation
StarDist: https://qupath.readthedocs.io/en/0.4/docs/deep/stardist.html
Watershed in QuPath
HoverNet: https://github.com/vqdang/hover_net

## MxIF cell segmentation
DeepCell: https://gist.github.com/petebankhead/db8548a0112bad089492061bf8046430
StarDist: https://qupath.readthedocs.io/en/0.4/docs/deep/stardist.html
Watershed in QuPath


## Other options and references:
DeepLIIF: An Online Platform for Quantification of Clinical ...
There are currently three tools that are available as web apps, 
namely DeepCell (https://deepcell.org/) [2, 6], 
Cellpose (https://www.cellpose.org/) [7], and 
NucleAIzer (https://www.nucleaizer.org/) [3]. 
DeepCell and Cellpose focus on nuclei and cytoplasm segmentation in multiplex images 
whereas Nucle-Alzer can segment nuclei across H&E, IHC, and multiplex images (via style transfer). 
DeepCell and Cellpose can also perform H&E and IHC nuclei segmentation if these are 
provided as grayscale images. DeepCell is the only tool among the above three that allows 
interaction with the segmented results via the DeepCell-Label web module (https://label.deepcell.org/). 
The Cellpose web app allows 512 × 512 input patches whereas DeepCell restricts input to 2048 × 2048 pixels 
in order to output results in a reasonable amount of time (less than a minute).

papers:
DeepLIIF: An Online Platform for Quantification of Clinical Pathology Slides
DeepSMILE: Contrastive self-supervised pre-training benefits MSI and HRD classification directly from H&E whole-slide images in colorectal and breast cancer

Question:
1. which tool for which modality?
2. how to evaluate the segmentation accuracy?














