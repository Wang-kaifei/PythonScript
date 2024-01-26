# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""脚本功能：从pFind.spectra文件中提取蛋白，然后再到原fasta数据库中查找，并放入新的数据库
   输入：卡的FDR值，pFind.spectra、原fasta文件路径、输出新fasta的路径
"""
"""脚本功能:将.protein文件中的主蛋白名称提取出来，需要提取的是在.fasta文件中有的名字"""
import os

def protein_group(Filename):
    """处理后缀为.protein的文件，提取符合要求的protein_group名，放入set中返回"""
    with open(Filename) as f1:
        f11 = f1.readlines()
    listA = []
    for i in range(0, len(f11)):
        if f11[i].split('\t', 3)[0].isdigit():
            if str.lower(f11[i].split('\t', 3)[1][0 : 3]) != "rev":
                listA.append(f11[i].split('\t', 3)[1])
                # print(f11[i].split('\t', 3)[1])
            else:
                print f11[i]
        elif "--" in f11[i].split('\t', 3)[0]:
            break
    return listA

def get_protein_name(Filename):
    """处理后缀为.fasta的文件，将蛋白质名称存储，以第一列为键，整行为值存成字典"""
    pro_dict = {}
    with open(Filename) as f1:
        f11 = f1.readlines()
    for i in range(0, len(f11)): #遍历整个文件
        if f11[i][0] == '>':   #如果该行是蛋白质名称
            index = f11[i].split(' ')[0]
            index = index[1 : len(index)] #去除>只保留名称
            value = f11[i]
            pro_dict[index] = value  #添加进字典
    return pro_dict

if __name__ == "__main__":
    protein = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\C\\wrong3.protein"
    fasta = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\fasta\\ecoli_plus_all.fasta"
    out = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\C\\wrong3_name.txt"
    names = protein_group(protein)
    dic = get_protein_name(fasta)
    print len(dic)
    print len(names)
    res = ""
    KeyError_name = set()
    for index in names:
        if index[0] != 'g' or index[1] != 'i':
            try:
                res += dic[index]
            except KeyError:
                KeyError_name.add(index)
    print  KeyError_name
    with open(out, "w") as f:
        f.write(res)
        f.close()
    



            
