#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""翻转数据库的decoy和target，（名字互换）"""



def rename(name, num):
    res = ""
    if ">rev" == name.lower()[0:4]:
        res = '>' + ''.join(name.strip().split()[0].split('_')[1:])
    else:
        res = '>REV_' + name.strip().split()[0][1:] + '\t' + '_REV'+ str(num)
    return res

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
        if file[i][0 : 7] == ">REV_sp":  #如果该行是decoy蛋白质名
            name = rename(file[i], count)
            name_dict[name] = file[i]
            res = res + name  + '\n'
            i += 1
            res += file[i]
    with open(outfile, "w") as f:
            f.write(res)
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
    fasta = "Y:\\chihao_dong111\\dong_111_raw\\ecoli+1_1human.fasta_td.fasta"
    out = "Y:\\chihao_dong111\\dong_111_raw\\ecoli+1_1human.fasta_td_reverse.fasta"
    name_dict = rewrite(fasta, out)
    dict_path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\uniprot\\name_convert.txt"
    write_dict(name_dict, dict_path)





