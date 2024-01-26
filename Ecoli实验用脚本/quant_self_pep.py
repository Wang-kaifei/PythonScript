#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 17:05:34 2020

@author: kaifeiwang
"""
"""脚本功能：将Nan文件中的肽段提取求其占整体肽段数量的比例"""
import os

def get_pep(path):
    """将path路径下的.spectra.list文件中的nan值提取重新生成文件，新文件只包括第0列"""
    for dirpath, dirnames, filenames in os.walk(path):
        for file in filenames:
            if file == "pQuant_spectra.list_nan":
                nanfile = os.path.join(dirpath, file)
    with open(nanfile) as f:
        f = f.readlines()
    peps = set()
    for i in range(0, len(f)):
        segs = f[i].rstrip().split('\t')
        peps.add(segs[1] + segs[2])
    summary_file = path + "\\pFind.summary"
    with open(summary_file) as f:
        f = f.readlines()
    all_pep = float(f[3].split()[-1])
    print len(peps)
    print all_pep
    return len(peps) / all_pep  #返回Nan比例
        
        
        
if __name__ == "__main__":
    """将所有Ecoli定量结果提取nan值生成新文件"""
    pathlist = ["C:\\test_wkf\\two_step_test\\Ecoli\\database1\\R1",
                "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\R2",
                "C:\\test_wkf\\two_step_test\\Ecoli\\vs_R_mod\\OR1_0.01",
                "C:\\test_wkf\\two_step_test\\Ecoli\\vs_R_mod\\OR1_0.05",
                "C:\\test_wkf\\two_step_test\\Ecoli\\vs_R_mod\\OR1_1",
                "C:\\test_wkf\\two_step_test\\Ecoli\\vs_R_mod\\OR2_0.01",
                "C:\\test_wkf\\two_step_test\\Ecoli\\vs_R_mod\\OR2_0.05",
                "C:\\test_wkf\\two_step_test\\Ecoli\\vs_R_mod\\OR2_1"
                ]
    outpath = "C:\\test_wkf\\two_step_test\\Ecoli\\quant_self.txt"
    res = "Test_Name\tNan_Ratio\n"
    for path in pathlist:
        res += path.rsplit('\\', 1)[1] + '\t' + str(get_pep(path)) + '\n'
    with open(outpath, 'w') as f:
        f.write(res)
        f.close()        
                         
                         
            
