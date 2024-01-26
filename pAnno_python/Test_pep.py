#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""统计输出结果中，被已注释库覆盖到的肽段数量"""

import re
import sys
from tqdm import tqdm
import time
from collections import Counter

class pFind_PSM:
    def __init__(self, File_Name="", Scan_No="", Exp_MHplus="", Charge="", Q_value="", Sequence="", Calc_MHplus="",
                 Mass_Shift="", Raw_Score="", Final_Score="", Modification="", Specificity="", Proteins="",
                 Positions="", Label="", Target="", Miss_Clv_Sites="", Avg_Frag_Mass_Shift="", Others=""):
        self.File_Name = File_Name
        self.Scan_No = Scan_No
        self.Exp_MHplus = Exp_MHplus
        self.Charge = Charge
        self.Q_value = Q_value
        self.Sequence = Sequence
        self.Calc_MHplus = Calc_MHplus
        self.Mass_Shift = Mass_Shift
        self.Raw_Score = Raw_Score
        self.Final_Score = Final_Score
        self.Modification = Modification
        self.Specificity = Specificity
        self.Proteins = Proteins
        self.Positions = Positions
        self.Label = Label
        self.Target = Target
        self.Miss_Clv_Sites = Miss_Clv_Sites
        self.Avg_Frag_Mass_Shift = Avg_Frag_Mass_Shift
        self.Others = Others

    def __str__(self):
        return str(self.File_Name) + '\t' + str(self.Scan_No) + '\t' + str(self.Exp_MHplus) + '\t' + str(
            self.Charge) + '\t' + str(self.Q_value) + '\t' + str(self.Sequence) + '\t' + str(
            self.Calc_MHplus) + '\t' + str(self.Mass_Shift) + '\t' + str(self.Raw_Score) + '\t' + str(
            self.Final_Score) + '\t' + str(self.Modification) + '\t' + str(self.Specificity) + '\t' + str(
            self.Proteins) + '\t' + str(self.Positions) + '\t' + str(self.Label) + '\t' + str(self.Target) + '\t' + str(
            self.Miss_Clv_Sites) + '\t' + str(self.Avg_Frag_Mass_Shift) + '\t' + str(self.Others) + '\n'
        # return '\t'.join(('%s' % item for item in self.__dict__.values())) + '\n'

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

def ReadpFind(res_file):
    """读取pFind鉴定到的sequence"""
    peps = set()
    with open(res_file) as file:
        lines = file.readlines()
    del lines[0]
    for line in lines:
        pep = line.split('\t')[5].strip()
        peps.add(pep)
    return peps

def Venn(peps1, peps2):
    print(f"peps1: {len(peps1)}\npeps2: {len(peps2)}\npeps1 & peps2: {len(peps1 & peps2)}\npeps1 - peps2: {len(peps1 - peps2)}\npeps2 - peps1: {len(peps2 - peps1)}\n")

def TestPep(test_peps, peps, pros):
    """统计鉴定到的肽段有多少被已注释库包含，
    先按特异性酶切匹配，未匹配上的肽段按照非特异性酶切考察
    peps: 已注释库酶切后生成的肽段
    pros: 已注释库包含的蛋白"""
    cov_cnt = 0
    test_peps_2 = [] # 存储不能按照特异性酶切匹配的肽段
    for test_pep in tqdm(test_peps, ascii=True):
        if IsContainRes(test_pep, peps):
            if cov_cnt < 10:
                print(test_pep)
            cov_cnt += 1
        else:
            test_peps_2.append(test_pep)
    for test_pep_2 in tqdm(test_peps_2, ascii=True): # 按照非特异性酶切考察
        if IsContainOpen(test_pep_2, pros):
            cov_cnt += 1
    tqdm.write(f"Test end!\nTotal peptide number: {len(test_peps)}.\nCovered peptide number: {cov_cnt}.\nUncovered peptide number: {len(test_peps) - cov_cnt}.")

def TestPSM(in_file, out_file, peps, pros):
    """考察PSM，提取不能被已注释库覆盖的PSM，并写出文件"""
    cov_cnt = 0
    in_f = open(in_file)
    out_f = open(out_file, 'w')
    title = in_f.readline() # 读第一行，列头
    out_f.write(title)
    line = in_f.readline()
    while line:
        line = in_f.readline()
        segs = line.strip().split('\t')
        if len(segs) < 19:
            continue
            print(segs)
        psm = pFind_PSM(segs[0], segs[1], segs[2], segs[3], segs[4], segs[5], segs[6], segs[7], segs[8], segs[9], segs[10], segs[11], segs[12], segs[13], segs[14], segs[15], segs[16], segs[17], segs[18])
        if not IsContainRes(psm.Sequence, peps) and not IsContainOpen(psm.Sequence, pros):
            out_f.write(str(psm))
            cov_cnt += 1
    print(f"cov_cnt : {cov_cnt}")
    in_f.close()
    out_f.close()

if __name__ == "__main__":
    full_res = r"G:\kfwang\human\database\pFind_Full\pFind-Filtered.spectra"
    res_file1 = r"G:\kfwang\human\PaceEva\2\output\pFind_res1\pFind-Filtered.spectra"
    res_file2 = r"G:\kfwang\human\PaceEva\2\output\pFind_res2\pFind-Filtered.spectra"
    out_file = r"G:\kfwang\human\PaceEva\2\output\uncovered.spectra"
    target_fasta_path = r"G:\kfwang\human\PaceEva\2\output\pFind_res2\DBReducer\DBReducer.fasta"
    total_peps = ReadpFind(full_res) # 待召回的全部肽段
    print(f"#Benchmark peptides: {len(total_peps)}")
    fir_peps = ReadpFind(res_file1) # p第一轮鉴定到的肽段
    sec_peps = ReadpFind(res_file2) # pAnno第二轮鉴定到的肽段
    panno_peps = sec_peps | fir_peps
    print(f"pAnno peps: {len(panno_peps)}")
    first_rec = total_peps & fir_peps
    test_peps = total_peps - (sec_peps | fir_peps) # benchmark 未召回结果
    print(f"Total recall: {len(total_peps) - len(test_peps)}")
    print(f"First recall peps: {len(first_rec)}")
    print(f"Second recall peps: {len(total_peps) - len(test_peps) - len(first_rec)}")
    print(f"Second should recall peps: {len(total_peps) - len(first_rec)}")
    print(f"test_peps : {len(test_peps)}")
    print(f"pAnno Res 1 precision: {len(first_rec) / len(fir_peps)}")
    print(f"pAnno Res 2 precision: {len(total_peps & sec_peps) / len(sec_peps)}")
    print(f"pAnno Res precision: {len(total_peps & panno_peps) / len(panno_peps)}")
    Venn(total_peps, test_peps)

    pros_target = ReadProSeq(target_fasta_path)
    peps_target = Digestion(pros_target) # 数据库模拟酶切后的肽段
    # test_peps = ReadpFind(res_file1)
    TestPep(test_peps, peps_target, pros_target) # 判断肽段被已注释库的覆盖情况
    TestPSM(res_file1, out_file, peps_target, pros_target) # 提取不能被已注释库覆盖的肽段对应的PSM输出


    # peps_target = Digestion(pros_target) # 数据库模拟酶切后的肽段
    # tqdm.write(f"Read target pep size: {len(peps_target)}")
    # cnt_use_target = 0
    # for pep in peps_target:
    #     if len(pep) >= 6:
    #         cnt_use_target += 1
    # tqdm.write(f"Read target pep size >= 6: {cnt_use_target}")
    # pros_decoy = ReadProSeq(decoy_fasta_path)
    # peps_decoy = Digestion(pros_decoy) # 数据库模拟酶切后的肽段
    # tqdm.write(f"Read decoy pep size: {len(peps_decoy)}")
    # cnt_use_decoy = 0
    # for pep in peps_decoy:
    #     if len(pep) >= 6:
    #         cnt_use_decoy += 1
    # tqdm.write(f"Read decoy pep size >= 6: {cnt_use_decoy}")
    # tqdm.write(f"ratio >= 6: {cnt_use_decoy / cnt_use_target}")
    # tqdm.write(f"ratio: {len(peps_decoy) / len(peps_target)}")
