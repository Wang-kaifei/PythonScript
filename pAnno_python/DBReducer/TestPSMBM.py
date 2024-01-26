#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""统计输出结果中，被已注释库覆盖到的肽段数量"""

import re
import os
import sys
from tqdm import tqdm
from tqdm_joblib import tqdm_joblib
import time
from collections import Counter
from joblib import Parallel, delayed
from functools import partial
import concurrent.futures
import multiprocessing
import ahocorasick
import math

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

"""从.spectra文件中读取PSM
返回PSM dict, key = spectra, value = sequence"""
def ReadPSM(pfind_res_path, threshold = 1):
    res = {} # spec 2 pep
    pep2code = {} # pep to godel code
    with open(pfind_res_path, 'r') as f:
        lines = f.readlines()
    del lines[0]
    for line in tqdm(lines, ascii=True):
        segs = line.strip().split('\t')
        if float(segs[4]) > threshold:
            break
        if len(segs) < 6:
            break
        res[segs[0].strip()] = segs[5].strip()
        pep2code[segs[5].strip()] = GetGodel(segs[5].strip())
    return res, pep2code

"""从PSMs中读取sequence"""
def GetSeqsbyPSM(psms):
    return set(psms.values())

"""从PSMs中读取MSMS"""
def GetMSMSbyPSM(psms):
    return set(psms.keys())

"""判断pep对应的psm是否为待考察psm"""
def ProcessData(fasta_peps, rec_peps, index_pep):
    index = index_pep[0]
    pep = index_pep[1]
    return index, (pep not in rec_peps) and (pep in fasta_peps)

"""判断pep对应的psm是否为待考察psm, 未召回的谱图"""
def ProcessDatatmp(fasta_peps, rec_peps, index_pep):
    index = index_pep[0]
    pep = index_pep[1]
    return index, (pep not in rec_peps)

"""获取benchmark结果中待考察PSM: 没有被召回的PSM
BMpfind_path: benchmark pfind鉴定结果
first_pfind_path: 删谱前第一轮搜索结果
DBRpfind_path: 测试实验 pfind鉴定结果(搜索缩减数据库)
fasta_path: DBReducer缩减后的库
"""
def BMTestPSM(BMpfind_path, DBRpfind_path, first_pfind_path, fasta_path):
    res_psm = {}
    BM_psms, BM_code = ReadPSM(BMpfind_path)
    test_psms, test_code = ReadPSM(DBRpfind_path)
    first_psm, first_code = ReadPSM(first_pfind_path)
    test_psms.update(first_psm)
    rec_peps = GetSeqsbyPSM(test_psms)
    print(f"# pAnno identified pep: {len(rec_peps)}")
    fasta_peps = Digestion(ReadProSeq(fasta_path))
    print(f"# DBR Database contain pep : {len(fasta_peps)}")
    #BM_peps = GetSeqsbyPSM(BM_psms)
    BM_peps = list(BM_psms.values())
    print(f"#Benchmark peptides: {len(BM_peps)}")
    #print(f"# pAnno recalled pep {len(rec_peps & BM_peps)}")
    ProcessDataWrapper = partial(ProcessDatatmp, fasta_peps, rec_peps)

    # TestAnnoRes_1 = partial(TestAnnoRes, name2seq, genes, anno_res)
    pool = multiprocessing.Pool(processes=50)
    results = list(tqdm(pool.imap(ProcessDataWrapper, [(i, d) for i, d in enumerate(BM_peps)], chunksize=5000), total=len(BM_peps), ascii=True))

    # results = list(tqdm(Parallel(n_jobs = 30)(delayed(ProcessData)(i, d) for i, d in enumerate(peps)), total=len(peps), ascii=True))
    print(results[0:10])
    # 根据任务的编号进行排序
    sorted_results = sorted(results, key=lambda x: x[0])
    # 提取排序后的结果
    sorted_results = [r[1] for r in sorted_results]
    i = 0
    for spec, pep in BM_psms.items():
        if sorted_results[i]:
            res_psm[spec] = pep
        i += 1
    return res_psm

"""读取.mgf文件中的谱图名"""
def ReadMGF(mgf_path):
    res = []
    with open(mgf_path, 'r') as f:
        lines = f.readlines()
    del lines[0]
    for line in tqdm(lines, ascii=True):
        if line[:5] != "TITLE":
            continue
        res.append(line.split('=')[1].strip())
    return res

def ContainSpec(test_specs, data_set_specs):
    for spec in test_specs:
        if spec in data_set_specs:
            return True
    return False

def PartBMPep(BMpfind_path, DBRpfind_path, first_pfind_path, fasta_path, mgf_path):
    data_set_specs = [] # pAnno删谱后剩余谱图
    cnt1 = 0 # 谱图被删除
    cnt2 = 0 # 谱图存在但匹配到的肽段发生变化
    cnt3 = 0 # PSM一致但未通过FDR过滤
    err_del_specs = [] # 被pAnno第一轮误删的谱图
    for filename in os.listdir(mgf_path):
        file_path = os.path.join(mgf_path, filename)
        if os.path.isfile(file_path):
            data_set_specs.extend(ReadMGF(file_path))
    print(f"# MGF Spectra: {len(data_set_specs)}")
    fasta_peps = Digestion(ReadProSeq(fasta_path))
    test_psm = BMTestPSM(BMpfind_path, DBRpfind_path, first_pfind_path, fasta_path) # benchmark待分类PSM
    test_pep_specs = {} # 将test_psm按照肽段整理，key = pep, value = spec list
    for spec, pep in test_psm.items():
        if pep not in test_pep_specs:
            test_pep_specs[pep] = [spec]
        else:
            test_pep_specs[pep].append(spec)
    print(f"# BM unrec PSM: {len(test_psm)}")
    cnt_no_spe_pep = 0 # 未召回的PSM且肽段谱图都不在
    print(f"# Fasta pepties: {len(fasta_peps)}")
    for spec, pep in tqdm(test_psm.items(), ascii=True):
        if (pep not in fasta_peps) and (spec not in data_set_specs):
            cnt_no_spe_pep += 1
    print(f"cnt_no_spe_pep : {cnt_no_spe_pep}")
    print(f"# BM unrec peptides: {len(GetSeqsbyPSM(test_psm))}")
    print(f"# pep dict length: {len(test_pep_specs)}")
    for spec, pep in test_pep_specs.items():
        if ContainSpec(specs, data_set_specs): # 有谱图存在
            cnt2 += 1
        else: # 对应的谱图都被删除
            cnt1 += 1
    print(f"In the benchmark results, there are {len(test_pep_specs)} peptides that are not recalled and which exist in the database.\n1. Spectrum was deleted: {cnt1}.\n2. Pep not identified: {cnt2}")

"""对于未被召回的分类
BMpfind_path: benchmark 搜索结果路径
DBRpfind_path: pAnno第二步搜索结果路径
first_pfind_path: pAnno第一步搜索结果路径
fasta_path: DBReducer缩库后的数据库
mgf_path: pAnno删谱后的数据集文件夹
""" 
def PartBMnRec(BMpfind_path, DBRpfind_path, first_pfind_path, fasta_path, mgf_path):
    specs = [] # pAnno删谱后剩余谱图
    cnt1 = 0 # 肽段和谱图都不存在
    cnt2 = 0 # 谱图存在肽段不存在
    cnt3 = 0 # 谱图、肽段都存在但发生变化
    cnt4 = 0 # PSM存在但未通过FDR过滤
    cnt5 = 0 # 肽段存在谱图不存在
    cnt6 = 0 # 谱图未解析 
    fasta_peps = Digestion(ReadProSeq(fasta_path)) # 第二次搜索的数据库
    err_del_specs = [] # 被pAnno第一轮误删的谱图
    part3spec = [] # part3 的谱图名
    for filename in os.listdir(mgf_path):
        file_path = os.path.join(mgf_path, filename)
        if os.path.isfile(file_path):
            specs.extend(ReadMGF(file_path))
    print(f"# MGF Spectra: {len(specs)}")
    test_psm = BMTestPSM(BMpfind_path, DBRpfind_path, first_pfind_path, fasta_path) # benchmark待分类PSM
    print(f"# BM unrec PSMs: {len(test_psm)}")
    print(f"# BM unrec peptides: {len(GetSeqsbyPSM(test_psm))}")
    pAnno_psm2, _ = ReadPSM(DBRpfind_path) # pAnno第二步搜索报告的PSM
    for spec, pep in test_psm.items():
        if spec not in specs: # 谱图没有在数据集中出现
            if pep not in fasta_peps: # part 1
                cnt1 += 1
            else: # part 5
                cnt5 += 1
                err_del_specs.append(spec)
        elif pep not in fasta_peps: # part 2
            cnt2 += 1
        elif spec not in pAnno_psm2:
            cnt6 += 1
        elif pep != pAnno_psm2.get(spec): # part 3 谱图对应的肽段发生变化
            part3spec.append(spec)
            cnt3 += 1
        else: # part 4
            cnt4 += 1
    print(f"In the benchmark results, there are {len(test_psm)} PSMs that are not recalled.\n1. Spectrum and peptide was deleted: {cnt1}.\n2. Datasets contain the spectrum but the peptied was deleted: {cnt2}.\n3. Pep changed: {cnt3}.\n4. Not filterd out: {cnt4}\n 5. Error deleted MS: {cnt5}\n 6. Spectrum was not analysed: {cnt6}")
    with open(os.path.join(os.path.dirname(BMpfind_path), "err_del.spec"), 'w') as file:
        for name in err_del_specs:
            file.write(name + "\n")
    return part3spec
    

def ReadSpecName(in_file):
    with open(in_file, 'r') as f:
        lines = f.readlines()
    return lines

"""将第一类，即pAnno第一次搜索结束后谱图被删除部分做细分
a. 在pAnno第一次搜索中被鉴定
b. 第一次中没有被鉴定，而是随着scan删除
spec_names: 第一类谱图名
psms: pAnno第一次搜索得到的psms
"""
def PartOne(spec_names, psms):
    cnt1 = 0
    cnt2 = 0
    print(len(psms))
    for name in spec_names:
        if name.strip() in psms:
            cnt1 += 1
        else:
            cnt2 += 1
    print(f"#Part1 spectra: {len(spec_names)}\n#Idetified in the first search of pAnno: {cnt1}\n#Delete by same scan {cnt2}")

def SubTitleMGF(mgf1, mgf2, out_path):
    title1 = set(ReadMGF(mgf1))
    title2 = set(ReadMGF(mgf2))
    print(f"#{mgf1}: {len(title1)}\n#{mgf2}: {len(title2)}\n#title1 - title2: {len(title1 - title2)}\n#title2 - title1: {len(title2 - title1)}\n")
    with open(out_path, 'w') as file:
        for name in (title1 - title2):
            file.write(name + "\n")

"""提取目标谱图的PSM，并写出
res_path1: pAnno第一轮搜索鉴定结果
res_path2: pAnno第二轮搜索鉴定结果
spec_name_path: 目标谱图名
out_path: 输出文件路径
"""
def GetPSMbySpecs(res_path1, res_path2, spec_name_path, out_path):
    specs = set()
    covered_specs = set()
    fout = open(os.path.join(out_path, "identiedpsms") , 'w') # 谱图被解析，输出PSM
    uncspec_file = open(os.path.join(out_path, "unany") , 'w') # 没有被解析的谱图，只输出谱图名
    with open(spec_name_path, 'r') as file:
        for line in file:
            specs.add(line.strip())
    print(f"#Specs: {len(specs)}")
    with open(res_path1, 'r') as file:
        lines = file.readlines()
    for line in tqdm(lines):
        segs = line.split('\t')
        if segs[0].strip() in specs:
            fout.write(line.strip() + "\t1\n")
            covered_specs.add(segs[0])
    ########################################
    with open(res_path2, 'r') as file:
        lines = file.readlines()
    for line in tqdm(lines):
        segs = line.split('\t')
        #if len(segs) < 7:
            #break
        if segs[0].strip() in specs:
            fout.write(line.strip() + "\t2\n")
            covered_specs.add(segs[0])

    for spec in specs:
        if spec not in covered_specs:
            uncspec_file.write(spec + '\n')
    uncspec_file.close()
    fout.close()    

def GetPepUncover(psm_r1, bm_code, sec_code):
    """求未被sec_code覆盖的肽段"""
    res = set() # 待召回的肽段序列
    print(f"#Second codes: {len(sec_code)}")
    for spec, pep in tqdm (psm_r1.items()):
        if bm_code[pep] not in sec_code:
            res.add(pep)
    return res

def GetInterPep(full_pep, full_codes, test_codes):
    res = set()
    for pep in full_pep:
        if full_codes[pep] in test_codes:
            res.add(pep)
    return res

def PepLevTest(BM_path, fir_path, sec_path, retainedDB_path):
    retain_pros = ReadProSeq(retainedDB_path) # 模拟已注释库蛋白质序列
    psm_r0, bm_code = ReadPSM(BM_path) # 全部已解析的谱图
    fir_psms, fir_code = ReadPSM(fir_path)
    sec_psms, sec_code = ReadPSM(sec_path)
    psm_r1 = FullPSMPart(retain_pros, psm_r0, set(fir_code.values())) # 全集中未被模拟已注释库覆盖到的PSM（待召回部分）
    part2bm = GetPepUncover(psm_r1, bm_code, set(fir_code.values())) # 待召回的肽段
    print(f"# Peptides to be recalled {len(part2bm)}\n")
    fir_seqs = GetSeqsbyPSM(fir_psms) # 第一轮鉴定到的肽段
    sec_seqs = GetSeqsbyPSM(sec_psms) - fir_seqs # 第二轮鉴定到的肽段
    print(f"# fir seqs {len(fir_seqs)}\n# sec seqs {len(sec_seqs)}")
    #test_seqs = fir_seqs | sec_seqs
    print(f"# Benchmark seq {len(part2bm)}\n# pAnno seq {len(sec_seqs)}")
    inter_seqs = GetInterPep(part2bm, bm_code, set(sec_code.values()))
    print(f"# inter seq {len(inter_seqs)}\n")
    print(f"precison={len(inter_seqs) / len(sec_seqs)}\nrecall={len(inter_seqs) / len(part2bm)}")
    return part2bm - inter_seqs

def TestOnePep (pep, pros):
    """判断肽段是否以子串形式存在于pros list"""
    for pro in pros:
        if pep in pro:
            return True
    return False

def TestPepConFast(peps, fasta_path):
    """判断肽段是否在fasta文件中存在，返回被fasta文件cover的肽段"""
    not_covered = set()
    db_pros = ReadProSeq(fasta_path)
    print(f"fasta size: {len(db_pros)}")
    print(f"#Total unrecalled peps: {len(peps)}")
    db_peps = Digestion(db_pros)
    for pep in peps:
        if pep in db_peps:
            continue
        elif not TestOnePep (pep, db_pros):
            not_covered.add(pep)
    print(f"uncovered pep : {len(not_covered)}")
    return peps - not_covered

def PSMFilterbyPro (trie, psms, pros):
    """过滤PSM，仅保留不被蛋白质序列包含的那些肽段对应的PSM"""
    matches = set() # 存储被覆盖到的肽段
    for pro in pros:
        for end_index, string in trie.iter(pro):
            matches.add(string)
    res = {}
    for spec, pep in psms.items():
        if pep not in matches:
            res[spec] = pep
    print(f"filtered psm: {len(res)}")
    return res

def FullPSMPart(retain_pros, total_psms, fir_codes):
    """对来自pFindFull的所有PSM，根据肽段是否对应到保留的模拟已注释库而分成两类"""
    peps = GetSeqsbyPSM(total_psms) # 提取所有的sequence
    trie = ahocorasick.Automaton() # Create a trie data structure
    for pep in peps:
        trie.add_word(pep, pep)
    trie.make_automaton()    # Build the trie data structure
    labeled_PSMs = PSMFilterbyPro (trie, total_psms, retain_pros) # 来自非模拟已注释库的PSM（recall分母）
    # 还需要加一个filter: 通过godel code 剔除第一轮被鉴定到的肽段，避免IL的影响
    for spec, pep in total_psms.items():
        if pep in fir_codes:
            del labeled_PSMs[spec]
    return labeled_PSMs

def GetPSMInter(psm_r1, bm_pep2code, sec_code):
    """求待召回部分和pAnno第二轮搜索结果的交集（recall、precision的分子）
    按照肽段的godel code确定来算"""
    res = {}
    for spec, pep in psm_r1.items():
        if bm_pep2code[pep] in sec_code:
            res[spec] = pep
    return res

def Spec2pFindRes(pfind_res_path):
    """key = spec, value = total line"""
    res = {}
    with open(pfind_res_path, 'r') as f:
        lines = f.readlines()
    del lines[0]
    for line in tqdm(lines, ascii=True):
        segs = line.strip().split('\t')
        if len(segs) < 6:
            break
        res[segs[0].strip()] = line.strip()
    return res
    
def PSMLevelEva(BMres_path, retainedDB_path, pFindres2, threshold, pFind_res1):
    bm_total_res = Spec2pFindRes(BMres_path) # 全部已解析的谱图对应的pFind结果
    not_rec = {} # 未被召回的部分
    """PSM层次的evaluation"""
    retain_pros = ReadProSeq(retainedDB_path) # 模拟已注释库蛋白质序列
    psm_r0, bm_code = ReadPSM(BMres_path) # 全部已解析的谱图
    annyed, fir_code = ReadPSM(pFind_res1) # 第一轮被解析的谱图
    psm_r1 = FullPSMPart(retain_pros, psm_r0, set(fir_code.values())) # 全集中未被模拟已注释库覆盖到的PSM（待召回部分）
    psm_r2, sec_code = ReadPSM(pFindres2, threshold) # pAnno第二轮搜索结果
    pep_r2 = GetSeqsbyPSM(psm_r2) # pAnno第二轮鉴定到的肽段
    print(f"Identified peptides by round 2: {len(pep_r2)}")
    msms_r2 = GetMSMSbyPSM(psm_r2) # pAnno第二轮搜索中被解析的谱图
    msms_r1 = GetMSMSbyPSM(psm_r1) # 待召回部分中被解析的谱图
    inter12 = GetPSMInter(psm_r1, bm_code, set(sec_code.values())) # 待召回部分和pAnno第二轮搜索结果的交集（recall、precision的分子）
    print(f"#Labeled PSM = {len(psm_r1)}\n#TP = {len(inter12)}\nRecall = {len(inter12)} / {len(psm_r1)} = {len(inter12) / len(psm_r1)}\nPrecision = {len(inter12)} / {len(msms_r2 & msms_r1)} = {len(inter12) / len(msms_r2 & msms_r1)}\n")
    for spec, line in bm_total_res.items():
        if spec in psm_r1 and spec not in inter12: #and spec not in annyed:
            not_rec[spec] = line
    print(f"#Not recalled PSM = {len(not_rec)}")
    WritePSM(not_rec, os.path.join(os.path.dirname(pFindres2), "notrec"))
    return not_rec

def WritePSM(psm, filename):
    with open(filename, 'w') as file:
        for spec, pep in psm.items():
            file.write(spec + "\t" + pep + "\n")
            
def SubPSM(res_file1, res_file2, outpath, ignor = None):
    psm1, code1 = ReadPSM(res_file1) # 全部已解析的谱图
    psm2, code2 = ReadPSM(res_file2) # 全部已解析的谱图
    print(code1["b1922_293T_proteinID_02A_QE3_122212.42936.42936.5.0.dta"])
    print(code2["b1922_293T_proteinID_02A_QE3_122212.42936.42936.5.0.dta"])
    sub_res1 = {}
    sub_res2 = {}
    for spec, code in code1.items():
        if spec in code2 and spec not in ignor:
            if code2[spec] != code:
                sub_res1[spec] = psm1[spec]
                sub_res2[spec] = psm2[spec]
    WritePSM(sub_res1, outpath + "\\psm1")
    WritePSM(sub_res1, outpath + "\\psm2")
            
if __name__ == "__main__":
    GodelInitialize()
    BMpfind_path = r"G:\kfwang\miniTest\database\pFind_Full\pFind-Filtered-full.spectra"
    retainedDB_path = r"G:\kfwang\miniTest\database\target_group_pace2.fasta"
    DBRpfind_path = r"G:\kfwang\miniTest\DBRPart\2\output\pFind_res7\pFind-Filtered.spectra"
    first_pfind_path = r"G:\kfwang\miniTest\DBRPart\2\output\pFind_res1\pFind-Filtered1.spectra"
    fasta_path = r"G:\kfwang\miniTest\DBRPart\2\output\pFind_res2\DBReducer\Reduced.fasta"
    mgf_path = "G:\\kfwang\\miniTest\\DBRPart\\2\\output\\pFind_res1\\MSdel\\"
    unrec_pep = PepLevTest(BMpfind_path, first_pfind_path, DBRpfind_path, retainedDB_path)
    recalled = PSMLevelEva(BMpfind_path, retainedDB_path, DBRpfind_path, 0.01, first_pfind_path) # 可以被召回的
    
    
    #WritePSM(ReadPSM(BMpfind_path)[0], os.path.join(os.path.dirname(BMpfind_path), "BMpsms"))
    # DBR_RNA = r"G:\kfwang\miniTest\DBRPart\2\output\pFind_res2\DBReducer\RNARed\pFind-Filtered.spectra"
    # fir_psm, _ = ReadPSM(first_pfind_path)
    # SubPSM(BMpfind_path, DBRpfind_path, os.path.dirname(DBRpfind_path), fir_psm.keys())
    
    db_coverd = TestPepConFast(unrec_pep, fasta_path)
    with open(os.path.join(os.path.dirname(DBRpfind_path), "db_cov"), 'w') as p3file:
        for item in db_coverd:
            p3file.write(item + "\n")
    # PartBMPep(BMpfind_path, DBRpfind_path, first_pfind_path, fasta_path, mgf_path)
    #part3spec = PartBMnRec(BMpfind_path, DBRpfind_path, first_pfind_path, fasta_path, mgf_path)
    #with open(os.path.join(os.path.dirname(DBRpfind_path), "part3sepc"), 'w') as p3file:
        #for item in part3spec:
            #p3file.write(item + "\n")
    #GetPSMbySpecs(first_pfind_path, DBRpfind_path, os.path.join(os.path.dirname(DBRpfind_path), "part3sepc"), os.path.dirname(DBRpfind_path))
    # spec_name_file = r"D:\users\kfwang\pAnno_2\human\database\pFind_Full\err_del.spec"
    # PartOne(ReadSpecName(spec_name_file), ReadPSM(first_pfind_path))

    # spec_name_file = r"D:\users\kfwang\pAnno_2\human\database\pFind_Full\err_del.spec"
    # PartOne(ReadSpecName(spec_name_file), ReadPSM(first_pfind_path))

    # mgf1 = r"D:\users\kfwang\pAnno_2\DelSpec\PaceEva\2\output\pFind_res2\MSdel\DEL_b1906_293T_proteinID_01A_QE3_122212_HCDFT.mgf"
    # mgf2 = r"D:\users\kfwang\pAnno_2\DelSpec\PaceEva\2\output\pFind_res2\MSdel_old\DEL_b1906_293T_proteinID_01A_QE3_122212_HCDFT.mgf"
    # out_path = r"D:\users\kfwang\pAnno_2\DelSpec\PaceEva\2\output\pFind_res2\MSdel\01A_sub.mgf"
    # SubTitleMGF(mgf1, mgf2, out_path)