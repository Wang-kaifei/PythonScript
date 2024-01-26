#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""comet实现二次搜索需要的建库脚本"""


def get_pro_name(spectrafile):
    """第10列是protein"""
    proteins = set()
    with open(spectrafile) as file:
        file = file.readlines()
    for i in range(1, len(file)):
        segs = file[i].split('\t')
        pros = segs[15].strip().split(',')
        for pro in pros:
            if pro == "":
                continue
            proteins.add(pro.strip())
    print(len(proteins))
    #test, print 10 proteins.
    i = 0
    for protein in proteins:
        if i  == 10:
            break
        print (protein)
        i += 1
    return proteins

def test(protein_names, protein):
    for pro in protein_names:
        if pro in protein:
            return True, pro
    return False, ""

def get_sub (protein_names, fasta_path, outfile):
    """函数功能：从路径为fasta_path的.fasta文件中提取与protein_names集合中相同的
    蛋白质信息，写入新的fasta文件。
    """
    res = ""
    count = 0
    with open(fasta_path) as file:
        file = file.readlines()
    i = 0
    while i < len(file):
        if file[i][0] == ">":  #如果该行是蛋白质名
            protein_name = file[i]
            a, b = test(protein_names, protein_name)
            if a:
                res += file[i]   #存储蛋白质信息
                count += 1
                i += 1
                while i < len(file) and file[i][0] != ">":
                    res += file[i]
                    i += 1
                protein_names.remove(b)
                continue
        i += 1
          
    print(count)
    with open(outfile, "w") as f:
        f.write(res)




if __name__ == "__main__":
    fastapath = "Z:\\kfwang\\meta_griffin\\database\\Galaxy6-[Merged_database_for_input.fasta].fasta"
    idpath = "C:\\Users\\kfwang\\results\\griffin\\one-step\\comet\\result\\result_filter.txt"
    outpath = "C:\\Users\\kfwang\\results\\griffin\\one-step\\comet\\result\\new.fasta"
    get_sub(get_pro_name(idpath), fastapath, outpath)




