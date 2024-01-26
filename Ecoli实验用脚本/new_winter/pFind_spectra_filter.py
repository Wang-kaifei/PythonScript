#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""spectra文件过滤后建库
1. 删除对应5个蛋白质以上的sequence"""

import math



def get_dict(threshold, spec_path):
    """构建字典，
    key为sequence，value为该序列出现的次数（多少谱图鉴定到这个sequence）
    key为sequence，value为该序列对应的蛋白质list
    """
    with open(spec_path) as file:
        file = file.readlines()  #以行形式读取spectra文件
    seq_count = {}
    seq_pro = {}
    tpro_count = 0.0
    psm = 0.0
    for i in range(1, len(file)):  #遍历每一行
        segs = file[i].split('\t')
        psm += 1
        if float(segs[4]) > threshold:
            break
        key = segs[5]
        if key in seq_count:  #如果已经存储该sequence
            seq_count[key] += 1  #计数加一
            tpro_count += len(seq_pro[key])
        else:
            seq_count[key] = 1
            proteins = segs[12].split('/')[0 : -1]
            seq_pro[key] = proteins
            tpro_count += len(proteins)

    tspec_count = psm / len(seq_count)
    tpro_count /= len(seq_count)
    print(tpro_count)
    print(tspec_count)
    print(psm)
    return seq_count, seq_pro, math.floor(tspec_count), math.ceil(tpro_count)

def get_protein (seq_count, seq_pro, tspec_count = 1, tpro_count = 4):
    """根据过滤规则得到建库蛋白质
    支持谱图数不得少于tspec_count，
    sequence支持蛋白质数不得多于tpro_count
    """
    proteins = set()
    filters = 0
    for key in seq_count:
        if seq_count[key] >= tspec_count or len(seq_pro[key]) <= tpro_count:
            for pro in seq_pro[key]:
                if pro[0:3] != 'REV':
                    proteins.add(pro)
        else:
            filters += 1
    print(filters)
                
    return proteins

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
    threshold = 0.05
    spec_path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\pFind\\pFind.spectra"
    fasta_path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\ecoli_plus_human_con.fasta"
    outfile = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\database\\pFind_2\\0.05.fasta"
    seq_count, seq_pro, tspec_count, tpro_count = get_dict(threshold, spec_path)
    print(tpro_count)
    get_sub(get_protein(seq_count, seq_pro, tspec_count = tspec_count, tpro_count = tpro_count), fasta_path, outfile)





