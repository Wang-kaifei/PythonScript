#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""pFind实现二次搜索需要的建库脚本，直接字符串匹配"""


import pandas as pd
import csv
from Tkinter import _flatten


import math

def same_mass(peptide):
    site = []
    pep = []
    for i in range(0, len(peptide)):
        if peptide[i] == 'L' or peptide[i] == 'I':
            site.append(i)
    for i in range(0, int(math.pow(2, len(site)))):
        temp = list(peptide)
        b = bin(i)[2:]  #得到二进制排列
        for j in range(0, len(site) - len(b)):
            b = "0" + b
        for j in range(0, len(b)):
            if b[j] == '1':
                temp[site[j]] = 'L'
            else:
                temp[site[j]] = 'I'
        pep.append("".join(temp))
   
    return pep


def get_uniq_seq(seqs):
    """对传入的seq存储为set，且处理同构体问题"""
    res = set()
    print(len(seqs))
    flag = 0
    for seq in seqs:
        #if 'I' in seq or 'L' in seq:
            #flag += 1
           # temps = same_mass(seq)
          #  for temp in temps:
         #       res.add(temp)
        #else:
        res.add(seq)
    #print(flag)
    print(len(res))
    return res

def get_peptides(filepath):
    """读取pfind.spectra文件中的sequence"""
    train=pd.read_csv(filepath, sep='\t', usecols=[5], header = None, skiprows=1)
    list_seq = list(_flatten(train.values.tolist()))
    set_seq = set()
    for seq in list_seq:
        set_seq.add(seq)
    return get_uniq_seq(set_seq)

def build_dic(fastapath):
    """将原始数据库建立词典"""
    fasta = {}
    with open(fastapath) as file:
        file = file.readlines()  #以行形式读取.fasta文件
        protein = file[0].split(' ')[0]
        sequence = ""
        start = 1
        while file[start][0] != '>':
            sequence += file[start].strip()
            start += 1
        for i in range(start, len(file)):  #遍历每一行
            if file[i][0] == '>':  #如果是蛋白质开头
                fasta[protein] = sequence  #存储上个蛋白
                protein = file[i].split(' ')[0]
                sequence = ""
            else:
                sequence += file[i].strip()
    print("build end")
    return fasta

def get_pro_name(fasta, peptides):
    """根据输入的fasta字典和序列信息提取被找到的蛋白质放入set"""
    proteins = set()
    for peptide in peptides:
        for key in list(fasta.keys()):
            if peptide in fasta[key]:
                proteins.add(key)
                del fasta[key] #剪枝，删除已存入set的蛋白质
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
        if file[i][0 : 1] == ">":  #如果该行是蛋白质名
            protein_name = file[i].split(' ')[0] #提取与get_protein函数匹配的形式
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
    fastapath = "Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"
    idpath = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind-Filtered.spectra"
    outpath = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\new.fasta"
    peptides = get_peptides(idpath)
    fasta = build_dic(fastapath)
    protein_names = get_pro_name(fasta, peptides)
    get_sub(protein_names, fastapath, outpath)




