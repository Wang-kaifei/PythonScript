#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""MaxQuant 建库结果概览"""



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


def get_pronum(proteins):
    """函数功能：获取fasta文件中蛋白质的个数"""
    gi = 0
    con = 0

    for protein in proteins:  #遍历每一行
        if protein[0:2].lower() == 'gi':
            gi += 1
        if protein[0:3].lower() == 'con':
            con += 1
    print(gi, con)


if __name__ == "__main__":
    path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\two_step\\pFind\\MaxQuant\\0.01\\combined\\txt\\proteinGroups.txt"
    protein = get_protein(path)
    get_pronum(protein)


