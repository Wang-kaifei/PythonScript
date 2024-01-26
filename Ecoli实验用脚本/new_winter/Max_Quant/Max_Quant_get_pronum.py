#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""得到Max_Quant的搜索结果进一步建库
原料：proteinGroups.txt文件，第一列是proteinnames，‘；’分割，末尾无额外加入"""



def get_protein(file_path):
    """从proteinGroups.txt文件读取蛋白质名"""
    res = set()
    with open(file_path) as file:
        file = file.readlines()
        for i in range(1, len(file)):  #遍历每一行
            proteins = file[i].split('\t')[0].split(';')
            for protein in proteins:
                if 'rev' != protein[0:3].lower():
                    res.add(protein)
    print(len(res))
    return res


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
    for i in range(0, len(file)):  #遍历每一行
        if file[i][0 : 1] == ">":  #如果该行是蛋白质名
            protein_name = file[i]
            a, b = test(protein_names, protein_name)
            if a:
                res += file[i]   #存储蛋白质信息
                count += 1
                i += 1
                while i < len(file) and file[i][0] != ">":
                    res += file[i]
                    i += 1
                i -= 1
                protein_names.remove(b)
    print(count)
    with open(outfile, "w") as f:
        f.write(res)


    


if __name__ == "__main__":
    path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\MaxQuant\\combined\\txt\\proteinGroups.txt"
    fasta_path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\uniprot\\ecoli_plus_all.fasta"
    out_path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\database\\Max_Quant\\Max_Quant.fasta"
    protein = get_protein(path)
    get_sub(protein, fasta_path, out_path)



