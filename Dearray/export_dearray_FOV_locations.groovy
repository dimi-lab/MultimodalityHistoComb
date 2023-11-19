import qupath.lib.regions.*
import ij.*
import java.awt.*
import java.awt.image.BufferedImage
import javax.imageio.ImageIO


def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()

// Request all objects from the hierarchy & filter only the annotations
def TMACores = hierarchy.getFlattenedObjectList(null).findAll{it.isTMACore()}

print(TMACores)
data_dir = "\\\\mfad\\researchmn\\HCPR\\HCPR-GYNECOLOGICALTUMORMICROENVIRONMENT\\Archive\\WSIs\\Ovarian_TMA"
case_id = new File(server.getMetadata().getName())
csv_fn = data_dir + File.separator + case_id +"_FOV_loc.csv"
fh = new File(csv_fn)
print(csv_fn)
String content = 'FOV,X,Y,W,H,Missing\n'

for (c in TMACores) {
    roi = c.getROI()
    line = c.toString() + "," + roi.getBoundsX().toString()\
    + "," + roi.getBoundsY().toString() + "," + roi.getBoundsWidth().toString()\
    + "," + roi.getBoundsHeight().toString()+ "," + c.isMissing().toString() + "\n"
    print(line)
    content += line
}
fh.write(content)


