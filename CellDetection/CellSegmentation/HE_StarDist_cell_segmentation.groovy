import qupath.ext.stardist.StarDist2D
import qupath.lib.objects.PathObjects
import static qupath.lib.gui.scripting.QPEx.*

// Specify the model file (you will need to change this!)
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

setPixelSizeMicrons(0.2201, 0.2201)
setImageType('BRIGHTFIELD_H_E')

def cellObj = getAnnotationObjects() //getDetectionObjects() // or getAnnotationObjects() if you converted your cells to annotations
removeObjects(cellObj, true)
def detObj = getDetectionObjects()
removeObjects(detObj, true)

createFullImageAnnotation(true)
// Run detection for the selected objects
def imageData = getCurrentImageData()
def pathObjects = getSelectedObjects()
if (pathObjects.isEmpty()) {
    Dialogs.showErrorMessage("StarDist", "Please select a parent object!")
    return
}
stardist.detectObjects(imageData, pathObjects)
println 'Done!'