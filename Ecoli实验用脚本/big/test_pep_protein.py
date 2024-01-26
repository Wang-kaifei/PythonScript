# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os

"""脚本功能：肽段层级陷阱库考察，若不含ecoli相关蛋白，则判定为非法
思路：每考察一行需要判定一下是否已经考察过，考察过的肽段放在set"""
def test_psm(proteins, pro_dict):
    """函数功能：考察该pep是否非法，输入为对应的蛋白质和提前存好的蛋白质的名称-信息字典对"""
    for i in range(0,len(proteins) - 1):
        if "REV" == proteins[i][0:3]:
            continue
        if "gi" == proteins[i][0:2] or "CON" == proteins[i][0:3] or "coli" in pro_dict[proteins[i]]:
            return True
    return False

def test(path, pro_dict):
    wrong = 0
    right = 0
    set_pep = set()
    with open(path) as f:
        f = f.readlines()
    for i in range(1, len(f)):
        segs = f[i].split('\t')
        try:
            if float(segs[4]) > 0.001:
                break
            if segs[15] == "decoy":
                continue
        except IndexError:
            continue
        try:
            pep = segs[5] + segs[10]
        except IndexError:
            continue
        if pep in set_pep:
            continue
        else:
            set_pep.add(pep)
            proteins = segs[12].split('/')
            if test_psm(proteins, pro_dict):
                right += 1
            else:
                wrong += 1
    print float(wrong) / (wrong + right)
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
    path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\big\\O_mod_100_6\\OR_0.01\\pFind.spectra"
    fasta = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\big\\fasta\\ecoli_plus_all.fasta"
    pro_dict = get_dict(fasta)
    right, wrong = test(path, pro_dict)
    print "right: ", right
    print "wrong: ", wrong
    
