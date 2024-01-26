# -*- coding: utf-8 -*-
"""只提取seq"""

def GetSeq(input_path, output_path):
    with open(input_path) as f1:
        f11 = f1.readlines()
    res = set()
    for i in range(0, len(f11)):
        if f11[i][0] == '>':
            continue
        res.add(f11[i])
    with open(output_path, "w") as f:
        for s in res:
            f.write(s)


if __name__ == "__main__":
    input_path = "Z:\\kfwang\\pAnno\\pt1\\All_GSSP.txt"
    output_path = "Z:\\kfwang\\pAnno\\pt1\\All_GSSP_only.txt"
    #input_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\open\\pt1_pFind_open_mgf_GSSP.txt"
    #output_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\open\\pt1_pFind_open_mgf_GSSP_only.txt"
    GetSeq(input_path, output_path)