# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os

"""脚本功能：数据库中CON开头蛋白质名称命名格式不符合Max_Quant的要求，故去掉"""
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
    fasta = "D:\\users\\kfwang\\two-step\\Ecoli\\standard\\database\\uniprot_unreviewed.fasta"
    out = "D:\\users\\kfwang\\two-step\\Ecoli\\Ecoli-unreviewed\\database\\uniprot_unreviewed(Max_Quant).fasta"
    res = ""
    with open(fasta) as f:
        f = f.readlines()
    i = 0
    count1 = 0
    count2 = 0
    while(i < len(f)):
        if f[i][0] != '>':
            res += f[i]
        elif f[i][1 : 4] == 'CON':
            count2 += 1
            i += 1
            while(f[i][0] != '>'):
                i += 1
            i -= 1
        else:
            count1 += 1
            res += f[i]
        i += 1
    with open(out, 'w') as f:
        f.write(res)
        f.close()
    print count1,count2

            
    
