# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os

"""脚本功能：MSFragger肽段层级陷阱库考察，若不含ecoli相关蛋白，则判定为非法
思路：每考察一行需要判定一下是否已经考察过，考察过的肽段放在set   8 17"""
def test_pep(proteins, pro_dict):
    """函数功能：考察该pep是否非法，输入为对应的蛋白质和提前存好的蛋白质的名称-信息字典对"""
    for i in range(0, len(proteins)):
        try:
            if proteins[i] == '' or "REV" == proteins[i][0:3]:
                continue
            if "gi" == proteins[i][0:2] or "CON" == proteins[i][0:3] or "coli" in pro_dict[proteins[i]]:
                return True
        except:
            continue
    return False

def test(path, pro_dict):
    wrong_pro = 0
    right_pro = 0
    wrong_lead = 0
    right_lead = 0
    set_pep = set()
    with open(path) as f:
        f = f.readlines()
    names = f[0].split('\t')
    index_pro = 0
    index_lead = 0
    for i in range(1, len(f)):
        segs = f[i].split('\t')
        proteins = segs[len(segs) - 1].split(', ')
        if test_pep(proteins, pro_dict):
            right_pro += 1
        else:
            wrong_pro += 1
        if test_pep([segs[8]], pro_dict):
            right_lead += 1
        else:
            wrong_lead += 1
    print "right_pro: ", right_pro
    print "wrong_pro: ", wrong_pro
    print "right_lead: ", right_lead
    print "wrong_lead: ", wrong_lead
    print "pro ratio:", float(wrong_pro) / (wrong_pro + right_pro)
    print "lead ratio:", float(wrong_lead) / (wrong_lead + right_lead)

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
    path = "C:\\Users\\kfwang\\test_wkf\\Meta_Frag\\peptide.tsv"
    fasta = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\big\\fasta\\ecoli_plus_all.fasta_td.fasta"
    pro_dict = get_dict(fasta)
    test(path, pro_dict)
