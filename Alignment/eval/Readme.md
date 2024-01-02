## Evaluate alignment accuracy

### Get ground truth affine transformations for selected samples
1. Create QuPath projects, add selected images into the project
2. Annotate corresponding landmark pairs
3. Export the annotated points to csv files ([export_anno_points.groovy](export_anno_points.groovy))
4. Read source and target points from the csv files, and calculate transformation based on annotation ([ground_truth_from_anno.ipynb](ground_truth_from_anno.ipynb))
5. For each FOV, save affine transformation into a file

**Note**: For both section 1 and section 2, both H&E stains were matched to MxIF from section 1

### Get affine transformations for selected samples
1. Create QuPath projects, add selected images into the project.
2. Run cell segmentation (StarDist [StarDist_cell_seg.groovy](StarDist_cell_seg.groovy), watershed [Watershed_cell_seg.groovy](Watershed_cell_seg.groovy)), and export cell measurements[ExportMeasurements.groovy](ExportMeasurements.groovy).
3. Read source and target points from the csv files, and use CPD to calculate transformation.
4. For each FOV, save affine transformation into a file. [cpd_alignment.ipynb](cpd_alignment.ipynb)

### Quantitatively evaluate alignment
1. Apply the "ground truth" affine transformation to the points and image [apply_alignment.ipynb](apply_alignment.ipynb)
2. Apply the "CPD" affine transformation to the points and image [apply_alignment.ipynb](apply_alignment.ipynb)
3. Quantitatively compare the alignment accuracy, using three metrics (rmse, theta, delta)


### Extended work 
1. Create better visualization https://github.com/BIOP/qupath-extension-biop




