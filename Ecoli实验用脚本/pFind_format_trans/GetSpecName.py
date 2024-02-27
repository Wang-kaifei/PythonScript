'''
Descripttion: 
version: 
Author: Kaifei
Date: 2024-01-25 17:14:52
LastEditors: Kaifei
LastEditTime: 2024-02-27 18:07:51
'''
# -*- coding: utf-8 -*-
#从.mgf文件中提取title行并存储

import os
import glob

def GetTitle(MGF_file, out_file):
    f = open(out_file, "w")
    with open(MGF_file, 'r') as ff:
        lines = ff.readlines()
    for line in lines:
        if line[0:5] == "TITLE":
            f.write(line.strip() + "\n")
    f.close()

if __name__ == "__main__":
    mgf_folder = "/Users/kaifeiwang/Desktop/mgf/"
    files_list = glob.glob(os.path.join(mgf_folder, '*.mgf'))
    for path in files_list:
        out_path = path.rsplit('.', 1)[0] + ".title"
        GetTitle(path, out_path)
