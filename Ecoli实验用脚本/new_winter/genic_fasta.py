#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""生成只包含蛋白质名称和sequence的通用fasta文件"""



def rename(name, num):
    res = ""
    if ">gi" == name.lower()[0:3]:
        res = '>G' + str(num)
    elif ">con" == name.lower()[0:4]:
        res = '>C' + str(num)
    else:
        res = '>P' + str(num)
    return res + '\n'

def rewrite(fasta_path, outfile):
    """
    将fasta文件的header只保留name输出
    """
    name_dict = {}
    res = ""
    count = 0
    with open(fasta_path) as file:
        file = file.readlines()
    for i in range(0, len(file)):  #遍历每一行
        if count % 500000 == 0:
            print(count)
        if file[i][0] == ">":  #如果该行是蛋白质名
            name = rename(file[i], count)
            res += name
            count += 1
        else:
            res += file[i]
    with open(outfile, "w") as f:
            f.write(res)
    print(count)
    return name_dict

def write_dict(dictionary, filepath):
    f = open(filepath, 'w')
    f.write(str(dictionary))
    f.close()

def read_dict(filepath):
    f = open(filepath, 'r')
    a = f.read()
    dictionary = eval(a)
    f.close()

    
if __name__ == "__main__":
    fasta = "D:\\users\\kfwang\\two-step\\Ecoli\\Ecoli-unreviewed\\database\\uniprot_unreviewed.fasta"
    out = "D:\\users\\kfwang\\two-step\\Ecoli\\Ecoli-unreviewed\\database\\uniprot_unreviewed_MaxQuant.fasta"
    name_dict = rewrite(fasta, out)
    dict_path = "D:\\users\\kfwang\\two-step\\Ecoli\\Ecoli-unreviewed\\database\\name_convert.txt"
    write_dict(name_dict, dict_path)





