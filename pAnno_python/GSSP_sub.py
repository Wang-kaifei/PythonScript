# -*- coding: utf-8 -*-
"""判断两个GSSP的提取结果是否对齐"""
def ReadSequence(GSSP_path):
    #提取蛋白库中的seq
    res = set()
    with open(GSSP_path) as f1:
        f11 = f1.readlines()
    i = 0
    while i < len(f11):
        if(f11[i][0] == '\n'):
            i += 1
            continue
        res.add(f11[i])
        i += 1
    return res

def TestSame(seq_seta, seq_setb):
    #判断两个set是否相同
    print(len(seq_seta))
    print(len(seq_setb))
    return len(seq_seta) == len(seq_setb) and len(seq_seta - seq_setb) == 0

if __name__ == "__main__":
    GSSP_patha = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\trypsin\\pFind_r_gssp.txt"
    GSSP_pathb = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\trypsin\\GSSP_restricted_cpp.txt"
    print(TestSame(ReadSequence(GSSP_patha), ReadSequence(GSSP_pathb)))