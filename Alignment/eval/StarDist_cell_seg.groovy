import qupath.ext.stardist.StarDist2D
import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.EllipseROI;
import qupath.lib.objects.PathDetectionObject
import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathCellObject

String seg_method = "StarDist"

def entry = getProjectEntry()
def imageData = getCurrentImageData()

String imgName = entry.getImageName().split("\\.")[0]
print(imgName)

def detObj = getDetectionObjects()
removeObjects(detObj, true)

createFullImageAnnotation(true)

if (imageData.isFluorescence()){
    //setPixelSizeMicrons(0.32499898, 0.32499898)
    //setImageType('FLUORESCENCE')
    def pathModel = '\\\\mfad\\researchmn\\HCPR\\HCPR-GYNECOLOGICALTUMORMICROENVIRONMENT\\Archive\\WSIs\\Ovarian_TMA\\QuPath_cell_segmentation_models\\dsb2018_heavy_augment.pb'

    def stardist = StarDist2D.builder(pathModel)
            .threshold(0.55)              // Probability (detection) threshold
            .channels('DAPI_AF_R01')            // Specify detection channel
            .normalizePercentiles(1, 99) // Percentile normalization
            .pixelSize(0.32499898)              // Resolution for detection
            .cellExpansion(5.0)
            .cellConstrainScale(1.5)
            .measureShape()
            .measureIntensity()
            .includeProbability(true)
            .build()

    // Run detection for the selected objects
    def pathObjects = getSelectedObjects()
    if (pathObjects.isEmpty()) {
        Dialogs.showErrorMessage("StarDist", "Please select a parent object!")
        return
    }
    stardist.detectObjects(imageData, pathObjects)
}
else if (imageData.getImageType() == ImageData.ImageType.BRIGHTFIELD_H_E){
    //setPixelSizeMicrons(0.2201, 0.2201)

    def pathModel = '\\\\mfad\\researchmn\\HCPR\\HCPR-GYNECOLOGICALTUMORMICROENVIRONMENT\\Archive\\WSIs\\Ovarian_TMA\\QuPath_cell_segmentation_models\\he_heavy_augment.pb'

    def stardist = StarDist2D.builder(pathModel)
          .threshold(0.5)              // Prediction threshold
          .normalizePercentiles(0, 95) // Percentile normalization
          .pixelSize(0.2201)              // Resolution for detection
          .cellExpansion(5.0)
          .measureShape()
          .measureIntensity()
          .includeProbability(true)
          .build()

    // Run detection for the selected objects
    def pathObjects = getSelectedObjects()
    if (pathObjects.isEmpty()) {
        Dialogs.showErrorMessage("StarDist", "Please select a parent object!")
        return
    }
    stardist.detectObjects(imageData, pathObjects)
}
else{
    throw new Exception("Doesn't support this modality")
}

println 'Done!'

