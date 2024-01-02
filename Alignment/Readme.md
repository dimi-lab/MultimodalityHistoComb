Align MxIF, H&E and IHC images to the same space, and 
then apply the affine transformation to the cell detections to 
establish cellular spatial correspondence.

Invert the MxIF image to match the foreground range of H&E images? 
Why the alignment fails? The cells are too dense? How about sample the cells to sparser?


Question:
1. which tool to use?
2. how to evaluate the segmentation accuracy?
    a. use manual annotation as ground truth
    b. calculate cell distances?



# For papers
## robustness
1. mask a region to simulate tissue damage (probreg_points_mask.py)
2. different segmentation methods

## quantitative evalution
1. compare the metrics for annotations and same and serial section
















