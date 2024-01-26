#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 17:05:34 2020

@author: kaifeiwang
"""
"""脚本功能：求nan文件中psm占总体psm的比例"""
import os
    
def get_ratio(path, merge):
    """读取nan文件，判断是否为psm，并返回psm中在nan文件中出现的比例"""
    '''for dirpath, dirnames, filenames in os.walk(path):
        for file in filenames:
            if file == "pQuant_spectra.list_rnan":
                nanfile = os.path.join(dirpath, file)'''
    nanfile = path + "\\pQuant_spectra.list_rnan"
    with open(nanfile) as f_nan:
        f_nan = f_nan.readlines()
    total = 1
    if not merge:
        f_spec = path +  "\\pFind.summary"
        with open(f_spec) as f_sum:
            f_sum = f_sum.readlines()
        total = float(f_sum[1].split()[1])
    print len(f_nan), total, len(f_nan) / total 
    return len(f_nan) / total

if __name__ == "__main__":
    pathlist = [
        "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\database2\\O_mod_100\\VS_R\\OR2_0.001",
        "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\database2\\O_mod_100\\VS_R\\OR2_0.01",

                ]
    #flag = [False, True, True, True, True, True, True, True, True]
    flag = [False, False, False, False, False, False, False, False]
    outpath = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\database2\\O_mod_100\\VS_R\\OR2_0.001\\quant_self.txt"
    res = "Test_Name\tNan_Ratio\n"
    for i in range(0, len(pathlist)):
       # res += pathlist[i].rsplit('\\', 2)[1] +  pathlist[i].rsplit('\\', 2)[2] + '\t' + str(get_ratio(pathlist[i], flag[i])) + '\n'
       res += pathlist[i].rsplit('\\', 2)[2] + '\t' + str(get_ratio(pathlist[i], flag[i])) + '\n'
    with open(outpath, 'w') as f:
        f.write(res)
        f.close()        
                         
                         
            
