De-array is a process that gets the tissue core out of the whole slide image, and then save the cores into FOV images.
Since the QuPath Dearrayer doesn't work well in MxIF image, we propose to de-array in H&E, 
and the affine the core locations to MxIF, so that we can minimize human interaction in directly 
annotating cores within MxIF. To achieve this, we only need to annotate 3 pairs of points in H&E and MxIF 
to calculate the transformation between the two images.

## 1. H&E de-array
    a. Use Dearrayer within QuPath. May need to manually adjust some parameters.
    b. Run groovy script to save the FOV images. [tma_dearray_to_individual.groovy](tma_dearray_to_individual.groovy)
   ![QuPath_dearrayer.png](img%2FQuPath_dearrayer.png)
## 2. Establish spatial correspondence between MxIF and H&E
    a. Manually annotate three pairs of landmarks. Within QuPath, choose the similar cells within the same tissue core.
    b. Export H&E TMA core locations out of QuPath, denote as xy. [export_dearray_FOV_locations.groovy](export_dearray_FOV_locations.groovy)
    c. Calculate affine transformation matrix (denote as M), and then apply affine transormation (M) to the core locations (xy). Denote the transformed locations as trans_xy. [calculate_affined_loc.py](calculate_affined_loc.py)
## 3. MxIF de-array
    a. Import trans_xy into QuPath (with openning MxIF image).   [import_dearray_FOV_locations.groovy](import_dearray_FOV_locations.groovy)
    b. Manually adjust locations and set core missing/valid. 
    c. Run save FOV locations. [export_dearray_FOV_locations.groovy](export_dearray_FOV_locations.groovy)
    d. Run groovy script to save the MxIF FOV images. [tma_dearray_to_individual.groovy](tma_dearray_to_individual.groovy)





