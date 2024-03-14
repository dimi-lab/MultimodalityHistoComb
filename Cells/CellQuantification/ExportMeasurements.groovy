import qupath.lib.regions.ImagePlane
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.EllipseROI;
import qupath.lib.objects.PathDetectionObject
import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathCellObject

def plane = ImagePlane.getPlane(0, 0)
def myProject = getProject()
//myProject.getImageList().subList(0,10).each{
myProject.getImageList().each{
    String imgName = it.getImageName().split("\\.")[0]
    print(imgName)
    def separator = "\t"
    def exportType = PathCellObject.class
    //def exportType = PathDetectionObject.class
    def outputFile = new File(buildFilePath(PROJECT_BASE_DIR, "export", imgName, imgName+"_QUANT.tsv" ))
    File newDir = new File(buildFilePath(PROJECT_BASE_DIR, "export", imgName))
    if (!newDir.exists()) {
        newDir.mkdirs()
    }
    // Create the measurementExporter and start the export
    def exporter  = new MeasurementExporter()
	.imageList([it]) // Images from which measurements will be exported
	.separator(separator) // Character that separates values
	.exportType(exportType) // Type of objects to export
	.exportMeasurements(outputFile) // Start the export process
}


println("")
println("Done.")