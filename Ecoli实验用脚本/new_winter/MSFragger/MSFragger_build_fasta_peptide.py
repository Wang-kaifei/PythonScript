#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""MSFragger根据peptide.tsv文件建库的脚本
问题：proteins不在一列，是在每行的最后一列"""


import pandas as pd
import csv
from Tkinter import _flatten
import math

def get_proteins(filepath):
    """读取peptide.tsv文件中的proteins"""
    pro = set()
    with open(filepath) as f:
        f = f.readlines()
    for i in range(1, len(f)):  #遍历每一行
        segs = f[i].strip().split('\t')[-1]
        proteins = segs.split(',')
        if len(proteins) == 1:
            proteins = [segs]
        for protein in proteins:
            protein = protein.strip()
            if protein == "":
                print("???")
                continue
            if protein[0:3] == "REV":
                continue
            pro.add(protein)
    print(len(pro))
    return pro

def test(protein_names, protein):
    for pro in protein_names:
        if pro in protein:
            return True, pro
    return False, pro

def get_sub (protein_names, fasta_path, outfile):
    """函数功能：从路径为fasta_path的.fasta文件中提取与protein_names集合中相同的
    蛋白质信息，写入新的fasta文件。
    """
    used = set()
    res = ""
    count = 0
    print(len(protein_names))
    with open(fasta_path) as file:
        file = file.readlines()
    for i in range(0, len(file)):  #遍历每一行
        if file[i][0 : 1] == ">":  #如果该行是蛋白质名
            #protein_name = file[i].split(' ')[0][1:] #提取与get_protein函数匹配的形式
            protein_name = file[i]
            #if protein_name in protein_names or protein_name.split('|')[0] in protein_names:  #如果该蛋白质名在候选集合中
            a, b = test(protein_names, protein_name)
            if a:
                #print("right")
                if b in used:
                    print(b)
                used.add(b)
                res += file[i]   #存储蛋白质信息
                count += 1
                i += 1
                while i < len(file) and file[i][0] != ">":
                    res += file[i]
                    i += 1
                i -= 1
    print(count)
    print(len(used))
    with open(outfile, "w") as f:
        f.write(res)


if __name__ == "__main__":
    protein_file = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\MSFragger\\peptide.tsv"
    fasta_path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\ecoli_plus_human_con.fasta"
    outfile = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\database\\MSFragger\\0.01.fasta"
    proteins = get_proteins(protein_file)
    get_sub(proteins, fasta_path, outfile)





