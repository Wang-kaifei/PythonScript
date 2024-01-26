#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""pFind与MSFragger鉴定结果的差集，主要得到pFind多鉴定到的那些PSM(sequence)"""


import pandas as pd
import csv
from Tkinter import _flatten
import math

def MSF_get_seq(filepath):
    """读取peptide.tsv文件中的sequence"""
    seq_cow = pd.read_csv(filepath, sep='\t', usecols=[0], header = None, skiprows=1)
    seqs = set()
    seq_list = list(_flatten(seq_cow.values.tolist()))
    for seq in seq_list:
        seqs.add(seq)
    print(len(seqs))
    return seqs

def pFind_sub(filepath, outpath1, outpath2, MSF_seq):
    """考察pFind-filterd.spectra文件中的PSM，如果不在MSF的set中，则写出"""
    res1 = ""
    res2 = ""
    with open(filepath) as file:
        file = file.readlines()
    res1 = file[0]
    res2 = file[0]
    for i in range(1, len(file)):  #遍历每一行
        seq = file[i].split('\t')[5]
        if seq not in MSF_seq:
            res1 += file[i]
        else:
            res2 += file[i]
    with open(outfile1, "w") as f:
        f.write(res1)
    with open(outfile2, "w") as f:
        f.write(res2)


if __name__ == "__main__":
    pfind = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind-Filtered.spectra"
    MSFragger = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_comp\\peptide.tsv"
    outfile1 = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind-MSF_sub.spectra"
    outfile2 = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind-MSF_inter.spectra"
    pFind_sub(pfind, outfile1, outfile2, MSF_get_seq(MSFragger))
    #get_sub(proteins, fasta_path, outfile)





