#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""从pFind.protein文件提取蛋白质建库！且只选Sameset"""

import math

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    return False

def get_protein (filepath, threshold):
    """从pFind.protein文件提取蛋白质建库
    每一个group行，第3列是q_value，第0列是数字。一旦第3列不是数字，则break
    """
    proteins = set()
    with open(filepath) as f:
        f = f.readlines()
    for i in range(0, len(f)):
        segs = f[i].split('\t')
        try:
            if is_number(segs[0]):
                if float(segs[3]) <= threshold:
                    if 'REV' != segs[1].strip()[0:3]:
                        proteins.add(segs[1].strip())
                    else:
                        print ("filter", segs[3], segs[0])
                else:
                    break
            elif is_number(segs[3]):
                if 'sameset' == segs[1].strip().lower() and 'REV' != segs[2].strip()[0:3]:
                    ##print("??")
                    proteins.add(segs[2].strip())
        except:
            print(f[i])
            continue
        
    print(len(proteins))
    return proteins

def test(protein_names, protein):
    for pro in protein_names:
        if pro == "":
            continue
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
    threshold = 0.05
    protein_path = "C:\\Users\\kfwang\\results\\griffin\\one-step\\pFind-no-coelute\\pFind.protein"    #"C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\pFind\\pFind.protein" 
    fasta_path = "Z:\\kfwang\\meta_griffin\\database\\Galaxy6-[Merged_database_for_input.fasta].fasta"         #"C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\ecoli_plus_human_con.fasta"
    outfile = "C:\\Users\\kfwang\\results\\griffin\\one-step\\pFind-no-coelute\\0.05.fasta"
    get_sub(get_protein(protein_path, threshold), fasta_path, outfile)