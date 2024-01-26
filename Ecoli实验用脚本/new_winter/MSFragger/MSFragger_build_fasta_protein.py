#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""pFind实现二次搜索需要的建库脚本，直接从protein文件中提取"""


import pandas as pd
import csv
from Tkinter import _flatten
import math

def get_proteins(filepath):
    """读取peptide.tsv文件中的sequence"""
    train=pd.read_csv(filepath, sep='\t', usecols=[2], header = None, skiprows=1)
    return list(_flatten(train.values.tolist()))

def get_sub (protein_names, fasta_path, outfile):
    """函数功能：从路径为fasta_path的.fasta文件中提取与protein_names集合中相同的
    蛋白质信息，写入新的fasta文件。
    """
    res = ""
    with open(fasta_path) as file:
        file = file.readlines()
    for i in range(0, len(file)):  #遍历每一行
        if file[i][0 : 1] == ">":  #如果该行是蛋白质名
            protein_name = file[i].split(' ')[0][1:] #提取与get_protein函数匹配的形式
            if protein_name in protein_names:  #如果该蛋白质名在候选集合中
                #print("right")
                res += file[i]   #存储蛋白质信息
                i += 1
                while i < len(file) and file[i][0] != ">":
                    res += file[i]
                    i += 1
                i -= 1
    with open(outfile, "w") as f:
        f.write(res)


if __name__ == "__main__":
    protein_file = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_comp\\protein.tsv"
    fasta_path = "Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"
    outfile = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_comp\\new.fasta"
    proteins = get_proteins(protein_file)
    print(len(proteins))
    get_sub(proteins, fasta_path, outfile)





