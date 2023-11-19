import qupath.lib.regions.*
import qupath.lib.objects.PathObjects
import qupath.lib.objects.hierarchy.DefaultTMAGrid
import ij.*
import java.awt.*
import java.awt.image.BufferedImage
import javax.imageio.ImageIO


def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()

csv_fn = "\\\\mfad\\researchmn\\HCPR\\HCPR-GYNECOLOGICALTUMORMICROENVIRONMENT\\Archive\\WSIs\\Ovarian_TMA\\Slide2050_24Plex.ome.tif - SLIDE-2050_c1.tif_FOV_loc.csv"

fh = new File(csv_fn)

int numHorizontal = 27
int numVertical = 13
def cores = []

fh.withReader('UTF-8') { reader ->
    def line = reader.readLine()
    while ((line = reader.readLine()) != null) {
        ele = line.tokenize( ',' )
        def FOV_ID = ele[0]
        def c_x = ele[1] as Double
        def c_y = ele[2] as Double
        def W = ele[3] as Double
        def H = ele[3] as Double
        def Missing = ele[5]
        print(FOV_ID)
//        print(c_x)
//        print(c_y)
//        print(Missing)
        if (Missing == "TRUE"){
//            TMACore = new TMACoreObject()
            TMACore = PathObjects.createTMACoreObject(c_x, c_y, W*0.73, H*0.73, true)
        }
        else{
//            TMACore = new TMACoreObject()
            TMACore = PathObjects.createTMACoreObject(c_x, c_y, W*0.73, H*0.73, false)
        }
        TMACore.setName(FOV_ID)
        cores << TMACore
    }
}


// Create & set the grid
def tmaGrid = new DefaultTMAGrid(cores, numHorizontal)
getCurrentHierarchy().setTMAGrid(tmaGrid)
relabelTMAGrid("1-27","A-M", true)
