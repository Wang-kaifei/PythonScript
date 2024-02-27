#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""肽段鉴定层次的评测脚本
1. 构造待召回肽段集。遍历可信肽段集，判断肽段是否能对应到保留蛋白质，如果不能则该肽段判定为待召回。
2. 构造可信谱图集。每个引擎各写一个函数，如果该PSM对应到的肽段是可信肽段，则谱图被判定为可信谱图。 注意，这里的可信谱图集是全局的，全部引擎通用
3. 构造标注集。引擎的鉴定结果，只有可信谱图对应的PSM，才被视为标注集成员，参与到评测
4. 评测。对于每一个引擎，计算其召回率、精确率
    precision = 标注集肽段∩待召回肽段 / 标注集肽段数
    recall = 标注集肽段∩待召回肽段 / 待召回肽段数"""

import re
import os
import sys
from tqdm import tqdm
import time
from collections import Counter
from functools import partial
import concurrent.futures
import multiprocessing
import ahocorasick
import math

def ReadCrePep(input_file):
    """读取可信肽段集"""
    with open(input_file, 'r') as f:
        lines = f.readlines()
    return set([line.strip() for line in lines])

def ReadProSeq(fasta_path):
    """从fasta文件中读取蛋白质序列"""
    seqs = []
    with open(fasta_path) as file:
        lines = file.readlines()
    seq = ""
    for line in lines:
        if line[0] == ">":
            if seq != "":
                seqs.append(seq)
            seq = ""
        else:
            seq += line.strip()
    if seq != "":
        seqs.append(seq)
    tqdm.write(f"Build end, all: {len(seqs)} protein sequences.")
    time.sleep(0.1)
    return seqs

def GetPep(protein:str):
    """生成蛋白质酶切结果，最大遗漏酶切为2"""
    res = []
    pep = ""
    for aa in protein:
        if aa == 'K' or aa == 'R':
            pep += aa
            res.append(pep)
            pep = ""
        else:
            pep += aa
    if pep != "":
        res.append(pep)
    length = len(res)
    for i in range(1, length): # 遗漏一个酶切位点
        res.append(res[i - 1] + res[i])
    for i in range(2, length): # 遗漏两个酶切位点
        res.append(res[i - 2] + res[i - 1] + res[i])
    return res
    
def Digestion(proteins):
    """将proteins模拟酶切，返回生成的肽段列表"""
    peps = []
    for protein in proteins:
        peps.extend(GetPep(protein))
    return set(peps)

def DigestionUnique(proteins):
    """将proteins模拟酶切，返回生成的Unique肽段列表"""
    peps = []
    for protein in proteins:
        peps.extend(GetPep(protein))
    # 统计每个元素出现的次数
    counts = Counter(peps)
    # 筛选只出现一次的元素
    unique_elements = [element for element, count in counts.items() if count == 1]
    return unique_elements

def IsContainRes(pep_test, peps):
    """按照特异性酶切规则判定的结果"""
    return pep_test in peps

def IsContainOpen(pep_test, pros):
    """按照非特异酶切判定"""
    for pro in pros:
        if pep_test in pro:
            return True
    return False

def TestOnePep (pep, pros):
    """判断肽段是否以子串形式存在于pros list"""
    for pro in pros:
        if pep in pro:
            return True
    return False

def ExpectPep(cre_pep_path, anno_pro_path):
    """构造期望被召回的肽段集合，如果可信肽段不被蛋白质覆盖，则为期望召回肽段"""
    exp_pep = set()
    cre_peps = ReadCrePep(cre_pep_path)
    anno_pros = ReadProSeq(anno_pro_path)
    anno_peps = Digestion(anno_pros) # 已注释库酶切后的肽段
    for cre_pep in cre_peps:
        if cre_pep in anno_peps:
            continue
        elif not TestOnePep (cre_pep, anno_pros):
            exp_pep.add(cre_pep)
    return exp_pep

def ReadPSMpFind(pfind_res_path):
    """读取pFind的PSM  考虑到代码复用，其他引擎的结果都转化成pFind的格式来处理"""
    res = {} # key = spectra, value = sequence
    with open(pfind_res_path, 'r') as f:
        lines = f.readlines()
    del lines[0]
    for line in tqdm(lines, ascii=True):
        segs = line.strip().split('\t')
        if len(segs) < 6:
            break
        res[segs[0].strip()] = segs[5].strip()
    return res

def BuildCreSpec(res_path_list, cre_peps):
    """构造可信谱图集"""
    cre_specs = set()
    for path in res_path_list: # 读取所有引擎的结果
        psms = ReadPSMpFind(path)
        for spec, pep in psms.items(): # 存储可信肽段对应到的谱图
            if pep in cre_peps:
                cre_specs.add(spec)
    return cre_specs

def PepLevleTest(cre_peps, cre_specs, engine_res_path):
    """肽段层次评测"""
    psms = ReadPSMpFind(engine_res_path)
    label_peps = set() # 标注集对应到的所有肽段
    for spec, pep in psms.items():
        if spec in cre_specs:
            label_peps.add(pep)
    tp_pep = label_peps & cre_peps # 被召回的肽段
    precision = len(tp_pep) / len(label_peps)
    recall = len(tp_pep) / len(cre_peps)
    print(f"the res of {engine_res_path}")
    print(f"tp_pep: {len(tp_pep)}, label_peps: {len(label_peps)}, cre_peps: {len(cre_peps)}")
    print(f"Precision: {precision}, Recall: {recall}")

            
if __name__ == "__main__":
    pass