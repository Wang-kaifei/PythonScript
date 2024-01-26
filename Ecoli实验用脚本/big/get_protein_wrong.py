# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
"""脚本功能:考察filtered3.protein文件中的所有蛋白质group的SameSet，输出开头为gi或名称中
带有coli的总数量"""
import os

def test_coli(index, dic):
    name = ""
    try:   
        name = dic[index]
    except KeyError:
        print (index)
    if 'coli' in name:
        return True
    return False


def test(group, dic):
    """处理蛋白质group，判断该蛋白质group是否是错误的"""
    main_name = group[0].split('\t')[1]  #主蛋白名称
    if main_name[0 : 2] == 'gi' or main_name[0 : 3] == 'CON' or test_coli(main_name, dic):  #主蛋白是对的，返回False
        return False
    for i in range(1, len(group)):
        segs = group[i].split('\t')
        if segs[1] == 'SameSet': #考察该SameSet
            if segs[2][0 : 2] == 'gi' or segs[2][0 : 3] == 'CON' or test_coli(segs[2], dic): #SameSet出现正确蛋白质，返回False
                return False
        else:
            break
    return True

def get_string(group, num):
    res = ""
    segs = group[0].split('\t')
    segs[0] = str(num)
    for i in range(0, len(segs) - 1):
        res += segs[i] + '\t'
    res += '\n'
    for i in range(1, len(group)):
        res += group[i]
    return res

def get_protein_filter(Filename, dic):
    """处理后缀为.protein的文件，将每个protein group存储依次判断是否被过滤
    返回最终输出字符串"""
    res = ""
    num = 0
    with open(Filename) as f:
        f = f.readlines()
    group_now = []
    for i in range(0, len(f)):
        segs = f[i].split('\t')
        if segs[0].isdigit():  #group头
            if len(group_now) != 0 and test(group_now, dic):
                num += 1
                res += get_string(group_now, num)
            group_now = [] #开始记录新的group 
            group_now.append(f[i])
        elif "--" in f[i].split('\t', 3)[0]:
            break
        else:
            group_now.append(f[i])
    return res

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
    protein = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\C\\pFind_filtered3.protein"
    fasta = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\fasta\\ecoli_plus_all.fasta"
    out = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\C\\wrong3.protein"
    dic = get_protein_name(fasta)
    res = get_protein_filter(protein, dic)
    with open(out, 'w') as f:
        f.write(res)
        f.close()
