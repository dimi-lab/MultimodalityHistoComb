
// export contour of the cell boundaries (polygons).

String root_dir = "H:\\temp"


def entry = getProjectEntry()

String imgName = entry.getImageName().split("\\.")[0]
print(imgName)
def outputFile = new File(buildFilePath(PROJECT_BASE_DIR, "export", imgName+ "_cell_poly.csv" ))
def detObj = getDetectionObjects()
wrt_str = "cell_id,point_str\n"
cell_id = 0
for (c in detObj) {
    cell_id_str = cell_id.toString()
   
    points = c.getROI().getAllPoints()
    point_str = points.toString()
    
    wrt_str += cell_id_str + "," + point_str + "\n"
    cell_id +=1
}

outputFile.write(wrt_str)

print("Done")
