To run our alignment, cells within MxIF and H&E images need to be detected respectively. We used StarDist as our baseline for both H&E and MxIF, as there are pretrained models available for both modalities. 

Cell segmentation methods we tried.
## H&E cell segmentation options
1. StarDist: https://qupath.readthedocs.io/en/0.4/docs/deep/stardist.html
2. Watershed in QuPath
3. HoverNet: https://github.com/vqdang/hover_net

## MxIF cell segmentation
1. DeepCell: https://gist.github.com/petebankhead/db8548a0112bad089492061bf8046430
2. StarDist: https://qupath.readthedocs.io/en/0.4/docs/deep/stardist.html
3. Watershed in QuPath


## TODO:
1. Evaluate the segmentation accuracy
2. which tool works best for

### Other options and references:
1. DeepLIIF: An Online Platform for Quantification of Clinical Pathology Slides
2. Cellpose (https://www.cellpose.org/)  
3. NucleAIzer (https://www.nucleaizer.org/)

### Notes:
DeepCell and Cellpose focus on nuclei and cytoplasm segmentation in multiplex images 
Nucle-Alzer can segment nuclei across H&E, IHC, and multiplex images (via style transfer). 
DeepCell and Cellpose can also perform H&E and IHC nuclei segmentation if these are 
provided as grayscale images. DeepCell is the only tool among the above three that allows 
interaction with the segmented results via the DeepCell-Label web module (https://label.deepcell.org/). 
The Cellpose web app allows 512 × 512 input patches whereas DeepCell restricts input to 2048 × 2048 pixels 
in order to output results in a reasonable amount of time (less than a minute).















