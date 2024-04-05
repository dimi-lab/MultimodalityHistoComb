import os

core_list_fn = "core_list.txt"
HE_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/HE_FOVs/same_section"
MxIF_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/MxIF_FOVs/Slide2050_24Plex"
output_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/test_align_wsi/Sec1"

core_list = open(core_list_fn, 'r').readlines()
wrt_str = ""
for i in core_list:
    i = i.strip()
    # command = "python cli_seg_plus.py -s %s/%s -t %s/%s -o %s -v\n" % (HE_dir, i, MxIF_dir, i, output_dir)
    command = "python cli_seg_plus.py -s %s/%s -t %s/%s -o %s \n" % (HE_dir, i, MxIF_dir, i, output_dir)
    wrt_str += command
fn = "./run_sec1.sh"
fp = open(fn, 'w')
fp.write(wrt_str)
fp.close()


output_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/test_align_wsi/Sec2"
wrt_str = ""
for i in core_list:
    i = i.strip()
    # command = "python cli_seg_plus.py -s %s/%s -t %s/%s -o %s -v\n" % (HE_dir, i, MxIF_dir, i, output_dir)
    command = "python cli_seg_plus.py -s %s/%s -t %s/%s -o %s \n" % (HE_dir, i, MxIF_dir, i, output_dir)
    wrt_str += command
fn = "./run_sec2.sh"
fp = open(fn, 'w')
fp.write(wrt_str)
fp.close()
