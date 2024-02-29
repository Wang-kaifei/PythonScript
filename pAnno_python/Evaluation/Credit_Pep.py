'''
Descripttion: 
version: 
Author: Kaifei
Date: 2024-01-23 21:43:35
LastEditors: Kaifei
LastEditTime: 2024-02-29 17:58:58
'''
# -*- coding: utf-8 -*-

"""构造可信肽段以及可信谱图集，不包括画图的部分"""

import math
from tqdm import tqdm

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

def GetCredCodes(code1, code2, code3, code4, code5):
    inter12 = code1 & code2
    inter13 = code1 & code3
    inter14 = code1 & code4
    inter15 = code1 & code5
    inter23 = code2 & code3
    inter24 = code2 & code4
    inter25 = code2 & code5
    inter34 = code3 & code4
    inter35 = code3 & code5
    inter45 = code4 & code5
    #求并集
    union = inter12 | inter13 | inter14 | inter15 | inter23 | inter24 | inter25 | inter34 | inter35 | inter45
    print(f"Credible pep-code: {len(union)}")
    return union

def ReadPSM(pfind_res_path):
    spec2pep = {} # spec 2 pep
    pep2code = {} # pep to godel code
    with open(pfind_res_path, 'r') as f:
        lines = f.readlines()
    del lines[0]
    for line in tqdm(lines, ascii=True):
        segs = line.strip().split('\t')
        if float(segs[4]) >= 0.01 or len(segs) < 6:
            break
        spec2pep[segs[0].strip()] = segs[5].strip()
        pep2code[segs[5].strip()] = GetGodel(segs[5].strip())
    return spec2pep, pep2code

def BuildCreSpec(res_path_list, cre_pepcodes, outpath):
    """构造可信谱图集"""
    cre_specs = set()
    for path in res_path_list: # 读取所有引擎的结果
        psms, pep2code = ReadPSM(path)
        for spec, pep in psms.items(): # 存储可信肽段对应到的谱图
            if pep2code[pep] in cre_pepcodes:
                cre_specs.add(spec)
    fout = open(outpath, 'w')
    for spec in cre_specs:
        fout.write(spec + "\n")

def GetCredPep(pep1, pep2, pep3, pep4, pep5, outpath, cre_code):
    cre_peps = set()
    for pep in pep1:
        if GetGodel(pep) in cre_code:
            cre_peps.add(pep)
    for pep in pep2:
        if GetGodel(pep) in cre_code:
            cre_peps.add(pep)
    for pep in pep3:
        if GetGodel(pep) in cre_code:
            cre_peps.add(pep)
    for pep in pep4:
        if GetGodel(pep) in cre_code:
            cre_peps.add(pep)
    for pep in pep5:
        if GetGodel(pep) in cre_code:
            cre_peps.add(pep)
    print(f"Credible pep: {len(cre_peps)}")
    # 将pep写出
    with open(outpath, 'w') as f:
        for i in cre_peps:
            f.write(i + '\n')
    return cre_peps

if __name__ == "__main__":
    GodelInitialize()
    comet_path = r"E:\wkf\Ecoli\credible_pep\Comet\res\comet.spectra"
    msgf_path = r"E:\wkf\Ecoli\credible_pep\MSGF+\msgf.spectra"
    maxquant_path = r"E:\wkf\Ecoli\credible_pep\MaxQuant\maxquant.spectra"
    pFind_path = r"E:\wkf\Ecoli\credible_pep\pFind\pFind-Filtered.spectra"
    msfragger_path = r"E:\wkf\Ecoli\credible_pep\MSFragger\msfragger.spectra"
    cre_pep_path = r"E:\wkf\Ecoli\credible_pep\credible_pep.txt"
    cre_spec_path = r"E:\wkf\Ecoli\credible_pep\credible_spec.txt"
    comet_psm, comet_code = ReadPSM(comet_path)
    msgf_psm, msgf_code = ReadPSM(msgf_path)
    maxquant_psm, maxquant_code = ReadPSM(maxquant_path)
    pfind_psm, pfind_code = ReadPSM(pFind_path)
    msfragger_psm, msfragger_code = ReadPSM(msfragger_path)
    cre_code = GetCredCodes(set(comet_code.values()), set(msgf_code.values()), set(maxquant_code.values()), set(pfind_code.values()), set(msfragger_code.values())) # 可信肽段的编码
    GetCredPep(set(comet_code.keys()), set(msgf_code.keys()), set(maxquant_code.keys()), set(pfind_code.keys()), set(msfragger_code.keys()), cre_pep_path, cre_code) # 写出可信肽段
    BuildCreSpec([comet_path, msgf_path, maxquant_path, pFind_path, msfragger_path], cre_code, cre_spec_path) # 写出可信谱图集