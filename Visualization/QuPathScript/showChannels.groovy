/**
 * Simple script for QuPath v0.3.2 that loops through all available channels,
 * displaying them in turn for 1 second.
 *
 * Written for https://forum.image.sc/t/how-do-i-script-qupath-to-change-the-displayed-channel/67216
 *
 * @author Pete Bankhead
 */


def viewer = getCurrentViewer()
def display = viewer.getImageDisplay()
def available = display.availableChannels()
for (int i = 0; i < available.size()-1; i++) {
//    display.selectedChannels.setAll(available[i])
    print(available[i])
    display.selectedChannels.setAll([available[i],available[i+1]])
    viewer.repaintEntireImage()
    Thread.sleep(1000)
}