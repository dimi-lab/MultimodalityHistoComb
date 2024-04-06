import os

core_list_fn = "core_list.txt"
HE_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/HE_FOVs/same_section"
MxIF_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/MxIF_FOVs/Slide2050_24Plex"
MxIF_quant_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/AlignmentEval/QuPathAnnoProj_MxIF/export"
HE_quant_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/AlignmentEval/QuPathAnnoProj_HE_Sec1/export"

output_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/test_align_seg_qupath/Sec1"
core_list = open(core_list_fn, 'r').readlines()
wrt_str = ""
for i in core_list:
    i = i.strip()
    _s = os.path.join(HE_dir, i)
    _t = os.path.join(MxIF_dir, i)
    _is = os.path.join(HE_quant_dir, i.replace(".tif", "_StarDist_QUANT.tsv"))
    _it = os.path.join(MxIF_quant_dir, i.replace(".tif", "_StarDist_QUANT.tsv"))
    command = "python cli_plus.py -s %s -t %s -is %s -it %s -o %s -v false\n" % (_s, _t, _is, _it, output_dir)
    wrt_str += command
fn = "./q_run_sec1.sh"
fp = open(fn, 'w')
fp.write(wrt_str)
fp.close()


output_dir = "/infodev1/non-phi-data/junjiang/Ovarian_TMA/test_align_seg_qupath/Sec2"
wrt_str = ""
for i in core_list:
    i = i.strip()
    _s = os.path.join(HE_dir, i)
    _t = os.path.join(MxIF_dir, i)
    _is = os.path.join(HE_quant_dir, i.replace(".tif", "_StarDist_QUANT.tsv"))
    _it = os.path.join(MxIF_quant_dir, i.replace(".tif", "_StarDist_QUANT.tsv"))
    command = "python cli_plus.py -s %s -t %s -is %s -it %s -o %s -v false\n" % (_s, _t, _is, _it, output_dir)
    wrt_str += command
fn = "./q_run_sec2.sh"
fp = open(fn, 'w')
fp.write(wrt_str)
fp.close()
