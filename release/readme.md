# MxIF and H&E image registration tool
This tool is developed to establishing cell-level correspondence for cells from MxIF and H&E.

## Install dependencies
```shell
pip install stardist # if detect cells on the fly
pip install tifffile
pip install probreg==0.3.7
pip install opencv-python
pip install numpy pandas matplotlib 

```
## Apply
### Start from MxIF and H&E cell segmentation results 
Require the tsv file exported from QuPath. The tsv file saves the coordinates and staining features of the cells from cell segmentation. Our workflow does not rely on a specific cell segmentation model. You can use any model as long as the segmentation results looks good to you.  
The original MxIF and H&E images (must be tiff file) were also needed to provide extra information (for example, the size of the target image) or read in the original image and save the transformed image.
You will also need specify a folder to save the transformed image.
Please refer to code [cli.py](cli.py)
```shell
python cli.py -s /opt/HE.tiff.cell_quant.tsv -t /opt/MxIF.tiff.cell_quant.tsv -is /opt/HE.tiff  -t /opt/MxIF.tiff -o /opt/temp
```
### Directly start from MxIF and H&E image
Different from the previous one. This client code includes the cell segmentation part in the workflow, which means the 
client takes less input, but it takes longer time to run. 
Please refer to code [cli_seg.py](cli_seg.py)
```shell
python cli_seg.py -is /opt/HE.tiff  -t /opt/MxIF.tiff -o /opt/temp
```



