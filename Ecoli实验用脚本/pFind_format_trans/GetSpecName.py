# -*- coding: utf-8 -*-
#从.mgf文件中提取title行并存储

def GetTitle(MGF_file, out_file):
    f = open(out_file, "w")
    with open(MGF_file, 'r') as ff:
        lines = ff.readlines()
    for line in lines:
        if line[0:5] == "TITLE":
            f.write(line.strip() + "\n")
    f.close()

if __name__ == "__main__":
    MGF_file = "/Users/kaifeiwang/Desktop/mgf/Ecoli-1to1to1-un-C13-N15-150mM-20150823_HCDFT.mgf"
    out_file = "/Users/kaifeiwang/Desktop/mgf/Ecoli-1to1to1-un-C13-N15-150mM-20150823_HCDFT.title"
    GetTitle(MGF_file, out_file)
