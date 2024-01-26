#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 17:05:34 2020

@author: kaifeiwang
"""
"""脚本功能：.spectra文件中的差集提取，求其中nan的比例"""
import os
    


def get_ratio(path, psms):
    """读取nan文件，判断是否为psm，并返回psm中在nan文件中出现的比例"""
    print(path)
    for dirpath, dirnames, filenames in os.walk(path):
        for file in filenames:
            if file == "pQuant_spectra.list_rnan":
                nanfile = os.path.join(dirpath, file)
    with open(nanfile) as f:
        f = f.readlines()
    num = 0
    for i in range(0, len(f)):
        segs = f[i].split('\t')
        try:
            mods = segs[2].split('|', 1)[1].split('|')
        except:
            print i, segs
        this_psm = segs[0] + segs[1]
        for j in range(0, len(mods) - 1):
            this_psm += mods[j]
        for psm in psms:
            if psm == this_psm:
                num += 1
    
    print len(f), float(num), len(psms), float(num) / len(psms)
    return float(num) / len(psms)


def get_psm(filename, merge = False):
    """函数功能：从.spectra文件中按照nan文件的格式提取psm"""
    if merge:
        #filename = filename.rsplit('\\', 1)[0] + "\\result5.spectra"
        filename += "\\temp.spectra"
    else:
        filename += "\\pFind-Filtered.spectra"
    with open(filename) as f:
        f = f.readlines()
    psms = []
    for i in range(1, len(f)):
        segs = f[i].split('\t')
        psm = segs[0] + segs[5]
        mods = segs[10].split(';')
        for j in range(0, len(mods) -  1):
            psm += mods[j]
        psms.append(psm)
    return psms



def get_sub(psmA, psmB):
    psmA = set(psmA)
    psmB = set(psmB)
    return psmA.difference(psmB)
  
        


if __name__ == "__main__":
    """将所有Ecoli定量结果提取nan值生成新文件"""

    pathlist = [
        "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\C",
       "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\VS_C_SVM\\OR_0.01",
        "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\VS_C_SVM\\OR_0.05",
        "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\VS_C_SVM\\OR_1",
                   ]
    #flag = [False, False, False, False, False, False, False, False]
    flag = [False, True, True, True, False, True, True, True]
    outpath = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\quant_sub.txt"
    res = "Test_Name\tNan_Ratio\n"
    psm_list = []
    name_list = []
    for i in range(0, len(pathlist)):
        psm_list.append(get_psm(pathlist[i], flag[i]))    #存储各路径下psm
        if flag[i]:    
            name_list.append(pathlist[i].rsplit('\\', 1)[1])  #存储各路径的实验名称
        else:
            name_list.append(pathlist[i].rsplit('\\', 1)[1])
    for i in range(0, len(psm_list) - 1):
        for j in range(i + 1, len(psm_list)):
            res += name_list[i] + ' vs ' + name_list[j] + '\t'
            sub = get_sub(psm_list[i], psm_list[j])
            res += str(get_ratio(pathlist[i], sub)) + '\n'
            res += name_list[j] + ' vs ' + name_list[i] + '\t'
            sub = get_sub(psm_list[j], psm_list[i])
            res += str(get_ratio(pathlist[j], sub)) + '\n'
    with open(outpath, 'w') as f:
        f.write(res)
        f.close()        
                         
                         

