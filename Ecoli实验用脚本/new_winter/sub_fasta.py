#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""构建两个数据库的差集输出"""


import pandas as pd
import csv
from Tkinter import _flatten
import math
def get_protein(fasta_path):
    """从fasta文件读取蛋白质名"""
    res = set()
    with open(fasta_path) as file:
        file = file.readlines()
        for i in range(0, len(file)):  #遍历每一行
            if file[i][0] == ">":  #如果该行是蛋白质名
                protein_name = file[i].split(' ', 2)[0] #提取与get_protein函数匹配的形式
                res.add(protein_name)
    return res

def get_sub (protein_names1, protein_names2, fasta_path, outfile):
    """求protein_names1 - protein_names2并写出"""
    res = ""
    sub1 = protein_names1 - protein_names2
    sub2 = protein_names2 - protein_names1
    print(len(sub1))
    print(len(sub2))
    with open(fasta_path) as file:
        file = file.readlines()
    print("len  ",len(file))
    for i in range(0, len(file)):  #遍历每一行
        if i % 100000 == 0:
            print(i)
        if file[i][0 : 1] == ">":  #如果该行是蛋白质名
            protein_name = file[i].split(' ', 2)[0] #提取与get_protein函数匹配的形式
            if protein_name in sub1:  #如果该蛋白质名在候选集合中
                print(protein_name)
                res += file[i]   #存储蛋白质信息
                i += 1
                while i < len(file) and file[i][0 : 1] != ">":
                    res += file[i]
                    i += 1
                i -= 1
    with open(outfile, "w") as f:
            f.write(res)


    


if __name__ == "__main__":
    fasta_path1 = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_comp\\new.fasta"
    fasta_path2 = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_comp\\nn.fasta"
    fasta_path3 = "Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"
    outfile = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_comp\\sub.fasta"
    protein_names1 = get_protein(fasta_path1)
    protein_names2 = get_protein(fasta_path2)
    get_sub(protein_names1, protein_names2, fasta_path3, outfile)





