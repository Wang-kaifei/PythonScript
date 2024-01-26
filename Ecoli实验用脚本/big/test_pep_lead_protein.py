# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os

"""脚本功能：使用pfind.protein文件计算lead protein陷阱库比例"""
def test_pep(proteins, pro_dict):
    """函数功能：考察该protein是否非法，输入为对应的蛋白质和提前存好的蛋白质的名称-信息字典对"""
    for protein in proteins:
        if "gi" == protein[0:2] or "CON" == protein[0:3] or "REV" == protein[0:3] or "coli" in pro_dict[protein]:
            return True
    return False

def test(path, pro_dict):
    wrong_lead = 0
    right_lead  = 0
    i = 2
    pep_set = set()
    with open(path) as f:
        f = f.readlines()
    while(i < len(f)):
        segs = f[i].split('\t', 5)
        if "---" in segs[0]:
            break
        if segs[0].isdigit():
            proteins = []
            proteins.append(segs[1])
            flag = False
            i += 1
            segs_pep = f[i].split('\t')
            while(not segs_pep[2].isdigit()):
                if segs_pep[1] == "SameSet":
                    proteins.append(segs_pep[2])
                i += 1
                segs_pep = f[i].split('\t')
            flag = test_pep(proteins, pro_dict)
            print proteins
            try:
                while(segs_pep[2].isdigit()):
                    pep = segs_pep[3] + segs_pep[8]
                    if(pep not in pep_set):
                        if flag:
                            right_lead += 1
                        else:
                            wrong_lead += 1
                        pep_set.add(pep)
                    i += 1
                    segs_pep = f[i].split('\t')
            except IndexError:
                print segs_pep
    print "right_lead: ", right_lead
    print "wrong_lead: ", wrong_lead
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
    path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\big\\C\\pFind.protein"
    fasta = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\big\\fasta\\ecoli_plus_all.fasta"
    pro_dict = get_dict(fasta)
    test(path, pro_dict)

    
