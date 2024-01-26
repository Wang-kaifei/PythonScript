# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os

"""脚本功能：考察.spectra文件中的每一行，判定该psm是否为非法，若不含ecoli相关蛋白，则判定为非法"""
def test_psm(proteins, pro_dict):
    """函数功能：考察该psm是否非法，输入为对应的蛋白质和提前存好的蛋白质的名称-信息字典对"""
    for i in range(0,len(proteins) - 1):
        if "REV" == proteins[i][0:3]:
            continue
        if "gi" == proteins[i][0:2] or "CON" == proteins[i][0:3] or "coli" in pro_dict[proteins[i]]:
            return True
    return False

def test(path, pro_dict):
    wrong = 0
    right = 0
    with open(path) as f:
        f = f.readlines()
    for i in range(1, len(f)):
        seg_12 = f[i].split('\t')[12]
        proteins = seg_12.split('/')
        if test_psm(proteins, pro_dict):
            right += 1
        else:
            wrong += 1
    print wrong // (wrong + right)
    return right, wrong

def get_dict(fasta):
    pro_dict = {}
    with open(fasta) as f:
        f = f.readlines()
    for i in range(0, len(f)):
        if f[i][0] == '>':
            segs = f[i].split(' ')
            pro_dict[segs[0][1:]] = f[i]
    return pro_dict
        

if __name__ == "__main__":
    path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\big\\VS_C_SVM\\OR_0.01\\pFind-Filtered.spectra"
    fasta = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\big\\fasta\\ecoli_plus_all.fasta"
    pro_dict = get_dict(fasta)
    right, wrong = test(path, pro_dict)
    print "right: ", right
    print "wrong: ", wrong
    
