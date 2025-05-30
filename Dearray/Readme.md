De-array is a process that gets the tissue cores out of the whole slide image, and then save the cores into FOV(field of view) images.
Since the QuPath Dearrayer doesn't work well in MxIF image, we propose to de-array in H&E, 
and the affine the core locations to MxIF, so that we can minimize human interaction in directly 
annotating cores within MxIF. To achieve this, we only need to annotate 3 pairs of points in H&E and MxIF 
to calculate the transformation between the two images. The detected TMA core grids within H&E can be transferred to MxIF using
our groovy script. Using the transferred core detection grids from H&E, most tissue cores can be covered within MxIF. Some could not be encapsulated well within the whole slide MxIF image, maybe because of regional tissue warp. Despite this, this method has minimized human efforts in locating TMA cores and extracting FOV images.

## 1. H&E de-array
a. Use Dearrayer within QuPath. May need to manually adjust some parameters.

b. Run groovy script to save the FOV images. [Check code](./tma_dearray_to_individual.groovy)

![QuPath_dearrayer.png](img%2FQuPath_dearrayer.png)
## 2. Establish spatial correspondence between MxIF and H&E
a. Manually annotate three pairs of landmarks. Within QuPath, choose the similar cells within the same tissue core.

b. Export H&E TMA core locations out of QuPath, denote as xy. [Check code](./export_dearray_FOV_locations.groovy)

c. Calculate affine transformation matrix (denote as M), and then apply affine transormation (M) to the core locations (xy). Denote the transformed locations as trans_xy. [Check code](./calculate_affined_loc.py)

## 3. MxIF de-array
a. Import trans_xy into QuPath (with openning MxIF image).   [Check code](./import_dearray_FOV_locations.groovy)

b. Manually adjust locations and set core missing/valid. 

c. Run save FOV locations. [Check code](./export_dearray_FOV_locations.groovy)

d. Run groovy script to save the MxIF FOV images. [Check code](./tma_dearray_to_individual.groovy)





