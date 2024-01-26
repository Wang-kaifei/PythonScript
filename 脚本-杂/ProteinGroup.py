'''
Descripttion: 
version: 
Author: kfwang
Date: 2022-04-19 10:20:01
LastEditors: sueRimn
LastEditTime: 2022-04-21 09:35:07
'''
# -*- coding: utf-8 -*-
# @author Wang kaifei

import pandas as pd
import sys
from tqdm import tqdm

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
    """从提供的protein--peptide文件中读取，将相同蛋白的行归为一个protein
    即变成protein--peptides

    Args:
        path (_type_): 输入文件路径

    Returns:
        _type_: 存放proteins的list
    """
    temp_dic = dict() # key = protein_id value = list[peptides]
    data = pd.read_excel(path) # 只支持一个sheet
    col_name = list(data)
    print(col_name)
    for row in data.itertuples():
        if getattr(row, col_name[0]) in temp_dic: # 如果蛋白名已经存在
            temp_dic[getattr(row, col_name[0])].append(getattr(row, col_name[1])) # 添加对应肽段
        else: # 如果蛋白名不存在
            temp_dic[getattr(row, col_name[0])] = [getattr(row, col_name[1])] # 添加该蛋白
    proteins = []
    print(len(temp_dic))
    for key, value in temp_dic.items():
        proteins.append(Protein(key, value))
    return proteins

def IsSub(protein:Protein, stored_peps:list)->bool:
    for new_pep in protein.peps:
        if new_pep not in stored_peps:
            return False
    return True

def IsContain(new_protein:Protein, lead_protein:Protein)->bool:
    for new_pep in new_protein.peps:
        if new_pep in lead_protein.peps:
            return True
    return False
    
    
def Classify(group_res : dict, protein : Protein, stored_peps : list)->bool:
    if IsSub(protein, stored_peps): # 如果是subset，需要将其分配到合适的group
        for key, value in group_res.items():
            if IsContain(protein, key): # 如果新基因是subset，则存入
                value.append(protein)
        return True
    return False


def PartGroup(proteins:list):
    """将蛋白质分组

    Args:
        proteins (_type_): 存放蛋白质的list
    """
    proteins.sort(reverse = True) # 根据包含肽段的数量从大到小排序
    group_res = dict() # key = Protein value = [subset or sameset]
    stored_peps = [] # 存储已经出现过的肽段
    i = 0
    for protein in tqdm(proteins):
        if (not Classify(group_res, protein, stored_peps)):  # 如果没有成功存入现有Group
            group_res[protein] = [] # 自立一个group
            for pep in protein.peps:
                if pep not in stored_peps:
                    stored_peps.append(pep)
        i += 1
    return group_res
                
                
def IsSameset(new_protein:Protein, lead_protein:Protein)->bool:
    if len(new_protein.peps) != len(lead_protein.peps):
        return False
    for new_pep in new_protein.peps:
        if new_pep not in lead_protein.peps:
            return False
    return True


def WriteGroup(out_path:str, group_res:dict):
    with open(out_path, "w") as f:
        f.write("Protein id\tPeptides\n")
        for lead_pro, sub_pros in group_res.items():
            print(len(sub_pros))
            f.write(str(lead_pro)) # 写出lead protein
            for sub_pro in sub_pros: # 写出其sameset和subset
                if IsSameset(sub_pro, lead_pro):
                    f.write("\tSameSet\t")
                else:
                    f.write("\tSubSet\t")
                f.write(str(sub_pro))


if __name__ == "__main__":
    if (len(sys.argv) != 3):
        print("Error: Wrong param number!")
        exit(0)
    test_data = sys.argv[1] # 待处理文件
    out_data = sys.argv[2] # 输出文件
    group_res = PartGroup(ReadProtein(test_data))
    WriteGroup(out_data, group_res)