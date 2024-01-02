
PathObjectHierarchy hierarchy = getCurrentHierarchy()
def annotations = hierarchy.getAnnotationObjects()
def entry = getProjectEntry()
img_fn = entry.getImageName()
print(img_fn)
img_id = img_fn.substring(0, img_fn.length()-4)
print(img_id)
String wrt_str = "x,y\n"

def PROJECT_BASE_DIR = QuPathGUI.getInstance().getProject().getBaseDirectory().toString()
String fn = buildFilePath(PROJECT_BASE_DIR, "export", img_id +"_align_anno.csv")

File fp = new File(fn)

for (anno in annotations) {
    points = anno.getROI().getAllPoints()
    for (p in points){
        points_str = String.format("%f,%f\n", p.x, p.y)
        wrt_str += points_str
    }
}
fp.write(wrt_str)
print("Done")