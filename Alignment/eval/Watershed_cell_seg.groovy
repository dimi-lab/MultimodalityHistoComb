import qupath.ext.stardist.StarDist2D
import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.EllipseROI;
import qupath.lib.objects.PathDetectionObject
import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathCellObject

String seg_method = "Watershed"

def entry = getProjectEntry()
def imageData = getCurrentImageData()

String imgName = entry.getImageName().split("\\.")[0]
print(imgName)

def detObj = getDetectionObjects()
removeObjects(detObj, true)

createFullImageAnnotation(true)

if (imageData.isFluorescence()){
    //setPixelSizeMicrons(0.32499898, 0.32499898)
    runPlugin('qupath.imagej.detect.cells.WatershedCellDetection',
    '{"DetectionChannel": "DAPI_AF_R01", requestedPixelSizeMicrons": 0.32499898, "backgroundRadiusMicrons": 8, "medianRadiusMicrons": 0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 10,  "watershedPostProcess": true,  "cellExpansionMicrons": 5,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

}
else if (imageData.getImageType() == ImageData.ImageType.BRIGHTFIELD_H_E){
    //setPixelSizeMicrons(0.2201, 0.2201)
    runPlugin('qupath.imagej.detect.cells.WatershedCellDetection',
    '{"requestedPixelSizeMicrons": 0.2201, "backgroundRadiusMicrons": 8, "medianRadiusMicrons": 0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 0.1,  "watershedPostProcess": true,  "cellExpansionMicrons": 5,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

}
else{
    throw new Exception("Doesn't support this modality")
}


println 'Done!'

