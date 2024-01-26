'''
Descripttion: 以.spectra文件中的蛋白质为准，将鉴定到的蛋白质归group【没有蛋白质q_value这一步了】
注意，不能让某个蛋白出现在多个protein group中，所以subset应该更严格一些，只有当肽段全部被包含时，才算subset
另外，此文件应该仅包含target100%文件中被覆盖的蛋白，所以最后写出还需要再做一次过滤
version: 
Author: kfwang
Date: 2022-04-19 10:20:01
LastEditors: sueRimn
LastEditTime: 2022-04-27 11:51:31
'''
# -*- coding: utf-8 -*-
# @author Wang kaifei


import pandas as pd
import sys
from tqdm import tqdm
import ahocorasick

class Protein:
    def __init__(self, name = "", peps: list = []):
        self.name= name
        self.peps = peps
    def __eq__(self, other):
        return self.name == other.name
    def __str__(self) -> str:
        return self.name + "\t" + ";".join(self.peps) + "\n"     
    def __le__(self, other):
        return len(self.peps) < len(other.peps)
    def __gt__(self, other):
        return len(self.peps) > len(other.peps)
    def __hash__(self) -> int:
        return hash(self.name)
    
def ReadProtein(path):
    """从spectra文件中读取，将相同蛋白的行归为一个protein
    即变成protein--peptides

    Args:
        path (_type_): 输入文件路径

    Returns:
        _type_: 存放proteins的list
    """
    temp_dic = dict() # key = protein_id value = set(peptides)
    with open(path, "r") as f:
        lines = f.readlines()
    del lines[0]
    for line in lines: # 遍历所有PSM
        segs = line.strip().split("\t")
        seq = segs[5].strip()
        pros = segs[12].strip().split('/')
        for pro in pros:
            if pro == "" or pro[0:3] == "REV":
                continue
            if pro in temp_dic: # 如果蛋白名已经存在
                if seq not in temp_dic[pro]:
                    temp_dic[pro].append(seq) # 添加对应肽段
            else: # 如果蛋白名不存在
                temp_dic[pro] = [seq] # 添加该蛋白
    proteins = []
    print(len(temp_dic))
    for key, value in temp_dic.items():
        proteins.append(Protein(key, value))
    return proteins


def IsContain(new_protein:Protein, lead_protein:Protein)->bool:
    """判断蛋白质是否符合同组的关系

    Args:
        new_protein (Protein): 待判断蛋白
        lead_protein (Protein): 主蛋白

    Returns:
        bool: 是否能分到此group
    """
    for new_pep in new_protein.peps:
        if new_pep not in lead_protein.peps: # 有肽段不被包含，则不在该组
            return False
    return True
    
    
def Classify(group_res : dict, protein : Protein)->bool:
    """将新蛋白质进组

    Args:
        group_res (dict): 目前的分组
        protein (Protein): 新蛋白

    Returns:
        bool: 新蛋白是否能分入现有组
    """
    flag = False
    for key, value in group_res.items():
        if IsContain(protein, key): # 如果新基因是subset，则存入
            value.append(protein)
            flag = True
    return flag


def PartGroup(proteins:list):
    """将蛋白质分组

    Args:
        proteins (_type_): 存放蛋白质的list
    """
    proteins.sort(reverse = True) # 根据包含肽段的数量从大到小排序
    group_res = dict() # key = Protein value = [subset or sameset]
    for protein in tqdm(proteins, ascii=True):
        if (not Classify(group_res, protein)):  # 如果没有成功存入现有Group
            group_res[protein] = [] # 自立一个group
    return group_res
                
def IsSameset(new_protein:Protein, lead_protein:Protein)->bool:
    if len(new_protein.peps) != len(lead_protein.peps):
        return False
    for new_pep in new_protein.peps:
        if new_pep not in lead_protein.peps:
            return False
    return True

def BuildACAutomaton(keywords):
    A = ahocorasick.Automaton()
    for idx, keyword in enumerate(keywords):
        A.add_word(keyword, (idx, keyword))
    A.make_automaton()
    return A

def CheckSubStr(pro_name, automaton):
    matches = set()
    for end_index, (idx, keyword) in automaton.iter(pro_name):
        matches.add(idx)
    return matches

def ReadName (fasta_file):
    """读取蛋白质名"""
    names = set()
    for line in open(fasta_file):
        if line[0] == '>':
            names.add(line.split()[0][1:])
    return names

def WriteGroup(out_path:str, group_res:dict, fasta_file:str):
    pros_name = ReadName (fasta_file)
    trie = BuildACAutomaton(pros_name)
    cnt_lead = 0
    write_name = set() # 临时存储被写出的蛋白名
    with open(out_path, "w") as f:
        f.write("ID\tProtein id\tPeptides\n")
        for lead_pro, sub_pros in group_res.items():
            if not CheckSubStr(lead_pro.name, trie):
                continue
            cnt_lead += 1
            f.write(str(cnt_lead) + "\t" + str(lead_pro)) # 写出lead protein
            cnt_sub = 0
            write_name.add(lead_pro.name)
            for sub_pro in sub_pros: # 写出其sameset和subset
                if not CheckSubStr(sub_pro.name, trie):
                    continue
                cnt_sub += 1
                write_name.add(sub_pro.name)
                if IsSameset(sub_pro, lead_pro):
                    f.write("\t" + str(cnt_sub) + "\tSameSet\t")
                else:
                    f.write("\t" + str(cnt_sub) + "\tSubSet\t")
                f.write(str(sub_pro))
    print(f"#Total pro: {len(write_name)}")

if __name__ == "__main__":
    test_data = r"D:\users\kfwang\pAnno_2\human\database\pFind_Full\pFind-Filtered.spectra" # 待处理文件
    out_data = r"D:\users\kfwang\pAnno_2\human\database\pFind_Full\mygroup.protein" # 输出文件
    target_fasta_path = r"D:\users\kfwang\pAnno_2\human\database\target100%.fasta" # 有质谱数据覆盖的蛋白质DB
    group_res = PartGroup(ReadProtein(test_data)) # 划分group
    WriteGroup(out_data, group_res, target_fasta_path) # 写出划分结果（只写靶子库覆盖的部分）