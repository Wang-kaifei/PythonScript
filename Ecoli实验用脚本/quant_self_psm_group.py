#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 17:05:34 2020

@author: kaifeiwang
"""
"""脚本功能：求nan文件中psm占总体psm的比例，特殊适用于求各个分组的比例"""
import os
    
def get_ratio(path):
    """读取nan文件，判断是否为psm，并返回psm中在nan文件中出现的比例"""
    for dirpath, dirnames, filenames in os.walk(path):
        for file in filenames:
            if file == "pQuant_spectra.list_rnan":
                nanfile = os.path.join(dirpath, file)
    with open(nanfile) as f_nan:
        f_nan = f_nan.readlines()
    total = 1
    name = path.rsplit('\\', 1)[1]
    spectra = path + "\\" + name + ".spectra"
    with open(spectra) as f:
        f = f.readlines()
    print len(f_nan), len(f) - 1, float(len(f_nan)) / (len(f) - 1)
    return float(len(f_nan)) / (len(f) - 1)

if __name__ == "__main__":
    pathlist = [
        "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_0.05\\group\\group0",
        "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_0.05\group\\group1",
        "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_0.05\\group\\group2",
        "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_0.05\\group\\group3"
                ]
    flag = [True, True, True, True, True]
    outpath = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\quant_self.txt"
    res = "Test_Name\tNan_Ratio\n"
    for i in range(0, len(pathlist)):
        res += pathlist[i].rsplit('\\', 1)[1] + '\t' + str(get_ratio(pathlist[i])) + '\n'
    with open(outpath, 'w') as f:
        f.write(res)
        f.close()        
                         
                         
            
