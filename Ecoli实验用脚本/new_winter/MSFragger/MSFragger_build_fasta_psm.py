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

def get_seq(filepath):
    """读取psm.tsv文件中的sequence"""
    seq_cow = pd.read_csv(filepath, sep='\t', usecols=[2], header = None, skiprows=1)
    seqs = set()
    seq_list = list(_flatten(seq_cow.values.tolist()))
    print(len(seq_list))
    for seq in seq_list:
        seqs.add(seq)
    print(len(seqs))
    return seqs

def get_seq_uni(filepath):
    """读取psm.tsv文件中的sequence,同时标记为unique肽段"""
    seq_cow = pd.read_csv(filepath, sep='\t', usecols=[2], header = None, skiprows=1)
    flag_cow = pd.read_csv(filepath, sep='\t', usecols=[24], header = None, skiprows=1)
    seqs = set()
    seq_list = list(_flatten(seq_cow.values.tolist()))
    flag_list = list(_flatten(flag_cow.values.tolist()))
    for i in range(0, len(seq_list)):
        if flag_list[i] == True:
            seqs.add(seq_list[i])
    print(len(seqs))
    return seqs

def get_protein(filepath):
    """读取psm.tsv文件中的proteins"""
    pro_cow = pd.read_csv(filepath, sep='\t', usecols=[34], header = None, skiprows=1)
    proteins = set()
    pro_cow_list = list(_flatten(pro_cow.values.tolist()))
    for line in pro_cow_list:
        try:
            pros = line.split(',')
            for pro in pros:
                pro = pro.strip()
                if pro[0:3] == "REV":
                    continue
                else:
                    proteins.add(pro.strip())
        except:
            continue
    print(len(proteins))
    return proteins

def get_sub (protein_names, fasta_path, outfile):
    """函数功能：从路径为fasta_path的.fasta文件中提取与protein_names集合中相同的
    蛋白质信息，写入新的fasta文件。
    """
    res = ""
    with open(fasta_path) as file:
        file = file.readlines()
    for i in range(0, len(file)):  #遍历每一行
        if file[i][0] == ">":  #如果该行是蛋白质名
            protein_name = file[i].split(' ')[0][1:] #提取与get_protein函数匹配的形式
            if protein_name in protein_names:  #如果该蛋白质名在候选集合中
                res += file[i]   #存储蛋白质信息
                i += 1
                while i < len(file) and file[i][0] != ">":
                    res += file[i]
                    i += 1
                i -= 1
    with open(outfile, "w") as f:
        f.write(res)


if __name__ == "__main__":
    psm_file = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\test\\psm.tsv"
    fasta_path = "Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"
    outfile = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\test\\new_psm.fasta"
    get_seq(psm_file)
    get_sub(get_protein(psm_file), fasta_path, outfile)





