#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""肽段鉴定层次的评测脚本，各引擎内部结果独立评测
1. 各引擎分别提取可信肽段（0.1%FDR, 搜索已注释库），为待召回肽段
2. 各引擎构造可信谱图集（在提取可信肽段时就存储）
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

ALPHABET_SIZE = 26
PRIME_SIZE = 500
m_lfCode = [[0.0 for j in range(PRIME_SIZE)] for i in range(256)]
aacode = [0, 1, 2, 3, 4, 5, 6, 7, 11, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
prime = [
		2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,
		53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,
		127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
		199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,
		283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,
		383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,
		467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,
		577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,
		661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,
		769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,
		877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,
		983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,
		1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,
		1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,
		1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427,
		1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,
		1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,1601,1607,1609,1613,
		1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,1709,1721,1723,1733,
		1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,1831,1847,1861,1867,
		1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,1949,1951,1973,1979,1987,
		1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063,2069,2081,2083,2087,
		2089,2099,2111,2113,2129,2131,2137,2141,2143,2153,2161,2179,2203,2207,2213,
		2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,2293,2297,2309,2311,2333,
		2339,2341,2347,2351,2357,2371,2377,2381,2383,2389,2393,2399,2411,2417,2423,
		2437,2441,2447,2459,2467,2473,2477,2503,2521,2531,2539,2543,2549,2551,2557,
		2579,2591,2593,2609,2617,2621,2633,2647,2657,2659,2663,2671,2677,2683,2687,
		2689,2693,2699,2707,2711,2713,2719,2729,2731,2741,2749,2753,2767,2777,2789,
		2791,2797,2801,2803,2819,2833,2837,2843,2851,2857,2861,2879,2887,2897,2903,
		2909,2917,2927,2939,2953,2957,2963,2969,2971,2999,3001,3011,3019,3023,3037,
		3041,3049,3061,3067,3079,3083,3089,3109,3119,3121,3137,3163,3167,3169,3181,
		3187,3191,3203,3209,3217,3221,3229,3251,3253,3257,3259,3271,3299,3301,3307,
		3313,3319,3323,3329,3331,3343,3347,3359,3361,3371,3373,3389,3391,3407,3413,
		3433,3449,3457,3461,3463,3467,3469,3491,3499,3511,3517,3527,3529,3533,3539,
		3541,3547,3557,3559,3571
]

def GodelInitialize():
    for i in range(256):
        for j in range(PRIME_SIZE):
            m_lfCode[i][j] = (i + 1) * math.log(prime[j])
    return True

def GetGodel(peptide):
    lfCode = 0.0
    for i in range(len(peptide)):
        lfCode += m_lfCode[aacode[ord(peptide[i]) - ord('A')]][i]
    return lfCode
"""----------------------------------------------------------------------------------------------"""

def Pep2Codes(peps):
    """生成肽段的Godel编码"""
    pep2codes = {}
    for pep in peps:
        pep2codes[pep] = GetGodel(pep)
    return pep2codes

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

def ExpectPep(cre_peps, anno_pros, anno_peps_dig):
    """构造期望被召回的肽段集合，如果可信肽段不被蛋白质覆盖，则为期望召回肽段"""
    exp_pep = set()
    for cre_pep in cre_peps:
        if cre_pep in anno_peps_dig:
            continue
        elif not TestOnePep (cre_pep, anno_pros):
            exp_pep.add(cre_pep)
    print(f"exp_pep cnt: {len(exp_pep)}")
    return exp_pep

def ReadPSM(pfind_res_path, threshold = 0.001):
    spec2pep = {} # spec 2 pep
    pep2code = {} # pep to godel code
    with open(pfind_res_path, 'r') as f:
        lines = f.readlines()
    del lines[0]
    for line in tqdm(lines, ascii=True):
        segs = line.strip().split('\t')
        if float(segs[4]) >= threshold or len(segs) < 6:
            break
        spec2pep[segs[0].strip()] = segs[5].strip()
        pep2code[segs[5].strip()] = GetGodel(segs[5].strip())
    return spec2pep, pep2code
    
def PepLevleTestpFind(exp_peps, cre_specs, engine_res_path, anno_pros, anno_peps_dig):
    """肽段层次评测
    期望被召回的肽段, 可信谱图集, 引擎结果路径"""
    psms, pep2code = ReadPSM(engine_res_path)
    exp_codes = set(Pep2Codes(exp_peps).values()) # 期望被召回的肽段的Godel编码
    label_peps = set() # 标注集对应到的所有肽段
    for spec, pep in psms.items():
        if spec in cre_specs:
            label_peps.add(pep)
    pre_ud_pep = set() # precision的分母，标注集肽段中，无法对应到保留蛋白的部分
    for pep in label_peps:
        if pep in anno_peps_dig:
            continue
        elif not TestOnePep (pep, anno_pros):
            pre_ud_pep.add(pep)
    tp_pep_cnt = 0
    recall_code = set()
    for test_pep in label_peps:
        if pep2code[test_pep] in exp_codes:
            recall_code.add(pep2code[test_pep])
            tp_pep_cnt += 1
    # 错误的PSM写出文件
    with open(engine_res_path.rsplit('.', 1)[0] + "_wrong.txt", 'w') as f:
        for spec, pep in psms.items():
            if spec in cre_specs and pep in pre_ud_pep and pep2code[pep] not in exp_codes:
                f.write(f"{spec}\t{pep}\n")
    precision = tp_pep_cnt / len(pre_ud_pep)
    recall = len(recall_code) / len(exp_codes)
    print(f"the res of {engine_res_path}")
    print(f"tp_pep: {tp_pep_cnt}, pre_ud_pep: {len(pre_ud_pep)}, exp_peps: {len(exp_peps)}")
    print(f"Precision: {precision}, Recall: {recall}")
        
def PepLevleTest(cre_peps, cre_specs, engine_res_path, psms, pep2code):
    """肽段层次评测
    期望被召回的肽段, 可信谱图集, 引擎结果路径"""
    exp_codes = set(Pep2Codes(cre_peps).values()) # 期望被召回的肽段的Godel编码
    label_peps = set() # 标注集对应到的所有肽段
    for spec, pep in psms.items():
        if spec in cre_specs:
            label_peps.add(pep)
    label_codes = set(Pep2Codes(label_peps).values()) # 标注集对应到的所有godel code
    recall_codes = set() # 被召回的code
    recall_peps = set() # 被召回的肽段
    for test_pep in label_peps:
        if pep2code[test_pep] in exp_codes:
            recall_codes.add(pep2code[test_pep])
            recall_peps.add(test_pep)    
    # 错误的PSM写出文件
    with open(engine_res_path.rsplit('.', 1)[0] + "_wrong.txt", 'w') as f:
        for spec, pep in psms.items():
            if spec in cre_specs and pep in label_peps and pep2code[pep] not in exp_codes:
                f.write(f"{spec}\t{pep}\n")
    precision = len(recall_codes) / len(label_codes)
    recall_code_level = len(recall_codes) / len(exp_codes)
    recall_pep_level = len(label_peps & cre_peps) / len(cre_peps)
    print(f"recall_pep_level: {recall_pep_level}")
    print(f"the res of {engine_res_path}")
    print(f"#tp codes: {len(recall_codes)}, labeled codes: {len(label_codes)}, labeled_peps: {len(label_peps)}, cre_peps: {len(cre_peps)}")
    print(f"Precision: {precision}, Recall_code_level: {recall_code_level}")
    
    
    
#ECOLI
    # pfind_res_path = r"E:\wkf\Ecoli\Evaluation_pep\pFind\pFind-Filtered.spectra"
    # comet_res_path = r"E:\wkf\Ecoli\Evaluation_pep\Comet\comet.spectra"
    # msgf_res_path = r"E:\wkf\Ecoli\Evaluation_pep\MSGF+\msgf.spectra"
    # msfragger_res_path = r"E:\wkf\Ecoli\Evaluation_pep\MSFragger\msfragger.spectra"
    # maxquant_res_path = r"E:\wkf\Ecoli\Evaluation_pep\MaxQuant\maxquant.spectra"
    # pAnno_res_path1 = r"E:\wkf\Ecoli\Evaluation_pep\pAnno\output\pFind_res1\pFind-Filtered.spectra"
    # pAnno_res_path2 = r"E:\wkf\Ecoli\Evaluation_pep\pAnno\output\pFind_res2\pFind-Filtered.spectra"
    # anno_pro_path = r"E:\wkf\Ecoli\Evaluation_pep\pAnno\reference_db\target_group_pace2.fasta"
    # #############################################################################################
    # pfind_cre_path = r"E:\wkf\Ecoli\credible_pep\pFind\pFind-Filtered.spectra"
    # comet_cre_path = r"E:\wkf\Ecoli\credible_pep\Comet\res\comet.spectra"
    # msgf_cre_path = r"E:\wkf\Ecoli\credible_pep\MSGF+\msgf.spectra"
    # msfragger_cre_path = r"E:\wkf\Ecoli\credible_pep\MSFragger\msfragger.spectra"
    # maxquant_cre_path = r"E:\wkf\Ecoli\credible_pep\MaxQuant\maxquant.spectra"
    
if __name__ == "__main__":
    GodelInitialize()
    pfind_res_path = r"E:\wkf\Ecoli\Evaluation_pep\pFind\pFind-Filtered.spectra"
    comet_res_path = r"E:\wkf\Ecoli\Evaluation_pep\Comet\comet.spectra"
    msgf_res_path = r"E:\wkf\Ecoli\Evaluation_pep\MSGF+\msgf.spectra"
    msfragger_res_path = r"E:\wkf\Ecoli\Evaluation_pep\MSFragger\msfragger.spectra"
    maxquant_res_path = r"E:\wkf\Ecoli\Evaluation_pep\MaxQuant\maxquant.spectra"
    pAnno_res_path1 = r"E:\wkf\Ecoli\Evaluation_pep\pAnno\output\pFind_res1\pFind-Filtered.spectra"
    pAnno_res_path2 = r"E:\wkf\Ecoli\Evaluation_pep\pAnno\output\pFind_res2\pFind-Filtered.spectra"
    anno_pro_path = r"E:\wkf\Ecoli\Evaluation_pep\pAnno\reference_db\target_group_pace2.fasta"
    #############################################################################################
    pfind_cre_path = r"E:\wkf\Ecoli\credible_pep\pFind\pFind-Filtered.spectra"
    comet_cre_path = r"E:\wkf\Ecoli\credible_pep\Comet\res\comet.spectra"
    msgf_cre_path = r"E:\wkf\Ecoli\credible_pep\MSGF+\msgf.spectra"
    msfragger_cre_path = r"E:\wkf\Ecoli\credible_pep\MSFragger\msfragger.spectra"
    maxquant_cre_path = r"E:\wkf\Ecoli\credible_pep\MaxQuant\maxquant.spectra"
    cre_path_list = [pfind_cre_path, comet_cre_path, msgf_cre_path, msfragger_cre_path, maxquant_cre_path]
    res_path_list = [pfind_res_path, comet_res_path, msgf_res_path, msfragger_res_path, maxquant_res_path]
    anno_pros = ReadProSeq(anno_pro_path) # 已注释库蛋白质序列
    anno_peps_dig = Digestion(anno_pros) # 已注释库酶切后的肽段
    for i in range(len(cre_path_list)):
        cre_psm, _ = ReadPSM(cre_path_list[i], 0.001) # 标注集
        psms, pep2code = ReadPSM(res_path_list[i], 0.01) # 评测结果
        PepLevleTest(set(cre_psm.values()), set(cre_psm.keys()), res_path_list[i], psms, pep2code)
    # pAnno的结果，需要综合两轮搜索结果
    cre_psm_pfind, _ = ReadPSM(pfind_cre_path, 0.001) # 标注集
    psm1, pep2code1 = ReadPSM(pAnno_res_path1, 0.01) # 评测结果(第一轮)
    psm2, pep2code2 = ReadPSM(pAnno_res_path2, 0.01) # 评测结果(第二轮)
    psm2.update(psm1)
    pep2code2.update(pep2code1)
    PepLevleTest(set(cre_psm_pfind.values()), set(cre_psm_pfind.keys()), pAnno_res_path1, psm2, pep2code2)