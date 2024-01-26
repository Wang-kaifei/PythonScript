'''
Descripttion: 
version: 
Author: sueRimn
Date: 2022-08-30 14:04:50
LastEditors: Kaifei
LastEditTime: 2023-08-18 18:50:19
'''
# -*- coding: utf-8 -*-

"""脚本功能：处理pAnno翻译成的数据库
将定制数据库中的已注释蛋白质切换

"""

"""构造name2seq dict"""
def BuildDic(fasta_path):
    name_seq = {}
    with open(fasta_path) as file:
        lines = file.readlines()
    name = ""
    seq = ""
    for line in lines:
        if line[0] == ">":
            if seq != "":
                name_seq[name] = seq
            name = line.strip()
            seq = ""
        else:
            seq += line.strip()
    if seq != "" and name != "":
        name_seq[name] = seq
    print("Build end, all:", len(name_seq), " proteins")
    return name_seq

def BuildNewTransDB(out_path, name2seq, in_path):
    f_out = open(out_path, 'w')
    with open(in_path) as file:
        lines = file.readlines()
    for line in lines:
        if line[1] == "$":
            break
        f_out.write(line)
    for name, seq in name2seq.items():
        new_name = ">$$" + name[1:] + "\n"
        f_out.write(new_name)
        f_out.write(seq + "\n")

# 删除已注释蛋白质
def DelAnnoDB(out_path, in_path):
    f_out = open(out_path, 'w')
    with open(in_path) as file:
        lines = file.readlines()
    for line in lines:
        if line[1] == "$":
            break
        f_out.write(line)

def GerTransBase(out_path, in_path):
    f_out = open(out_path, 'w')
    with open(in_path) as file:
        lines = file.readlines()
    for line in lines:
        if line[1] == "$":
            break
        f_out.write(line)

#为数据库拼接已注释蛋白质
def CatAnnoDB(in_path, name2seq):
    """为数据库拼接"""
    f_out = open(in_path, 'a') # 追加形式打开文件
    for name, seq in name2seq.items():
        new_name = ">$$" + name[1:] + "\n"
        f_out.write(new_name)
        f_out.write(seq + "\n")

if __name__ == "__main__":
    ref_db = r"G:\kfwang\human\database\target_group_pace2.fasta" # 要拼接的已注释库
    old_db = r"G:\kfwang\human\PaceEva\2\trans_res\database.fasta" # 旧定制数据库
    out_path = r"G:\kfwang\human\PaceEva\2\trans_res\database_new.fasta" # 新结果
    name2seq = BuildDic(ref_db)
    # BuildNewTransDB(out_path, name2seq, old_db)
    #DelAnnoDB(out_path, old_db)

    #trans_out = r"D:\users\kfwang\pAnno\yeast-dbr\3\trans_res\trans.fasta"
    DBR_path = r"G:\kfwang\human\PaceEva\2\output\pFind_res2\DBReducer\DBReducer.fasta"
    # GerTransBase(trans_out, old_db)

    CatAnnoDB(DBR_path, name2seq)