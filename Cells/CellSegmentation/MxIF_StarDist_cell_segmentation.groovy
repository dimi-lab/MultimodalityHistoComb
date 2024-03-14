import qupath.ext.stardist.StarDist2D

// Specify the model file (you will need to change this!)
def pathModel = '\\\\mfad\\researchmn\\HCPR\\HCPR-GYNECOLOGICALTUMORMICROENVIRONMENT\\Archive\\WSIs\\Ovarian_TMA\\QuPath_cell_segmentation_models\\dsb2018_heavy_augment.pb'

def stardist = StarDist2D.builder(pathModel)
        .threshold(0.5)              // Probability (detection) threshold
        .channels('DAPI_AF_R01')            // Specify detection channel
        .normalizePercentiles(1, 99) // Percentile normalization
        .pixelSize(0.32499898)              // Resolution for detection
        .cellExpansion(5.0)
        .measureShape()
        .measureIntensity()
        .includeProbability(true)
        .build()

setPixelSizeMicrons(0.32499898, 0.32499898)
setImageType('FLUORESCENCE')

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