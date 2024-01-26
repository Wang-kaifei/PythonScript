#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 17:05:34 2020

@author: kaifeiwang
"""
"""把quant文件里的nan值提取出来生成新文件"""
import os

def get_nan(path):
    """将path路径下的.spectra.list文件中的nan值提取重新生成文件，新文件只包括第0列"""
    for dirpath, dirnames, filenames in os.walk(path):
        for file in filenames:
            if file == "pQuant.spectra.list":
                quantfile = os.path.join(dirpath, file)
    with open(quantfile) as f1:
        f11 = f1.readlines()
    res = ""  #记录结果
    cnt = 0   #记录结果数量
    title = 1 #记录干扰行数量
    for i in range(1, len(f11)):
        f11_sp = f11[i].split('\t')
        try: 
            r1 = float(f11_sp[9])
            r2 = float(f11_sp[18])
        except Exception:
            print("the next 30000")
            title += 1
            continue
        if (r1 == 1024 or r1 - 0.0009765625 < abs(1e-10)) and (r2 == 1024 or r2 - 0.0009765625 < abs(1e-10)):
            res += f11_sp[0] + '\t' + f11_sp[1] + '\t' + f11_sp[2] + '\n'
            cnt += 1
    newfile = path + "\\pQuant_spectra.list_rnan"
    with open(newfile, 'w') as f:
        f.write(res)
        f.close()
        
if __name__ == "__main__":
    """将所有Ecoli定量结果提取nan值生成新文件"""
    pathlist = [
        "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\database2\\O_mod_100\\VS_R\\OR2_0.001"

            ]
    for path in pathlist:
        print path
        get_nan(path)
                         
                         
# "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database2\\VS_C\\OR2_0.01\\result7",
 #       "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database2\\VS_C\\OR2_0.05\\result7",
  #      "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database2\\VS_C\\OR2_1\\result7",
   #     "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_0.01\\result7",
    #    "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_0.05\\result7",
     #   "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_1\\result7"       
