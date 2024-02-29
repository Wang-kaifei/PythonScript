'''
Descripttion: 
version: 
Author: Kaifei
Date: 2024-01-23 21:43:35
LastEditors: Kaifei
LastEditTime: 2024-02-29 16:56:37
'''
# -*- coding: utf-8 -*-

from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import seaborn as sns
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

def CometPep(res_path):
    seqs = set()
    codes = set()
    with open(res_path, 'r') as f:
        ff = f.readlines()
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        seqs.add(segs[11])
        codes.add(GetGodel(segs[11]))
    print(f"Comet: {len(seqs)}")
    return seqs, "Comet", codes

def msgf_mod_seq(raw_pep):
    mod_dic = {"+15.995":"Oxidation[M]", "+57.021":"Carbamidomethyl[C]", "-17.027":"Gln->pyro-Glu[AnyN-termQ]", "+42.011":"Acetyl[ProteinN-term]"}
    # 如果最左端出现（N端）修饰，则index是0
    mods = ""
    seq = "" # pFind格式的序列
    aa_index = 0 # 目前有多少个氨基酸出现
    temp_mod = "" # 暂时存储的修饰
    mid = raw_pep.split('.', 1)[1].rsplit('.', 1)[0]
    for cc in mid:
        if cc.isalpha():
            if temp_mod != "": # 如果前面的氨基酸有修饰，需要存储
                mods += str(aa_index) + "," + mod_dic[temp_mod] + ";"
                temp_mod = ""
            seq += cc
            aa_index += 1
        elif cc == "+" or cc == "-":
            if temp_mod != "": # 如果前面的氨基酸有修饰，需要存储
                mods += str(aa_index) + "," + mod_dic[temp_mod] + ";"
            temp_mod = cc
        else:
            temp_mod += cc
    return seq, mods
  
def MSGFPep(res_path):
    seqs = set()
    codes = set()
    with open(res_path, 'r') as f:
        ff = f.readlines()
    for i in range(1, len(ff)):    
        segs = ff[i].strip().split('\t')
        seqence, _ = msgf_mod_seq(segs[9])
        seqs.add(seqence)
        codes.add(GetGodel(seqence))
    print(f"MSGF+: {len(seqs)}")
    return seqs, "MSGF+", codes

def MaxquantPep(res_path):
    seqs = set()
    codes = set()
    with open(res_path, 'r') as f:
        ff = f.readlines()
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        seqs.add(segs[3])
        codes.add(GetGodel(segs[3]))
    print(f"Maxquant: {len(seqs)}")
    return seqs, "Maxquant", codes

def pFindPep(res_path):
    seqs = set()
    codes = set()
    with open(res_path, 'r') as f:
        ff = f.readlines()
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        seqs.add(segs[5])
        codes.add(GetGodel(segs[5]))
    print(f"pFind: {len(seqs)}")
    return seqs, "pFind", codes

def MSFraggerPep(res_path):
    seqs = set()
    codes = set()
    with open(res_path, 'r') as f:
        ff = f.readlines()
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        seqs.add(segs[0])
        codes.add(GetGodel(segs[0]))
    print(f"MSFragger: {len(seqs)}")
    return seqs, "MSFragger", codes

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
    res = {} # spec 2 pep
    pep2code = {} # pep to godel code
    with open(pfind_res_path, 'r') as f:
        lines = f.readlines()
    del lines[0]
    for line in tqdm(lines, ascii=True):
        segs = line.strip().split('\t')
        if float(segs[4]) >= 0.01 or len(segs) < 6:
            break
        res[segs[0].strip()] = segs[5].strip()
        pep2code[segs[5].strip()] = GetGodel(segs[5].strip())
    return res, pep2code

def BuildCreSpec(res_path_list, cre_pepcodes):
    """构造可信谱图集"""
    cre_specs = set()
    for path in res_path_list: # 读取所有引擎的结果
        psms, pep2code = ReadPSM(path)
        for spec, pep in psms.items(): # 存储可信肽段对应到的谱图
            if pep2code[pep] in cre_pepcodes:
                cre_specs.add(spec)
    return cre_specs

def WriteCreSpec(out_path, cre_specs):
    """可信谱图名写出文件"""
    fout = open(out_path, 'w')
    for spec in cre_specs:
        fout.write(spec + "\n")
        
def Get1idPep(pep1, pep2, pep3, pep4, pep5):
    """得到至少1个引擎鉴定的肽段"""
    union = pep1 | pep2 | pep3 | pep4 | pep5
    print(f"1id pep: {len(union)}")
    return union

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

def Get3idPep(pep1, pep2, pep3, pep4, pep5):
    """得到至少3个引擎鉴定的肽段"""
    inter123 = pep1 & pep2 & pep3
    inter124 = pep1 & pep2 & pep4
    inter125 = pep1 & pep2 & pep5
    inter134 = pep1 & pep3 & pep4
    inter135 = pep1 & pep3 & pep5
    inter145 = pep1 & pep4 & pep5
    inter234 = pep2 & pep3 & pep4
    inter235 = pep2 & pep3 & pep5
    inter245 = pep2 & pep4 & pep5
    inter345 = pep3 & pep4 & pep5
    #求并集
    union = inter123 | inter124 | inter125 | inter134 | inter135 | inter145 | inter234 | inter235 | inter245 | inter345
    print(f"3id pep: {len(union)}")
    return union

def Get4idPep(pep1, pep2, pep3, pep4, pep5):
    """得到至少4个引擎鉴定的肽段"""
    inter1234 = pep1 & pep2 & pep3 & pep4
    inter1235 = pep1 & pep2 & pep3 & pep5
    inter1245 = pep1 & pep2 & pep4 & pep5
    inter1345 = pep1 & pep3 & pep4 & pep5
    inter2345 = pep2 & pep3 & pep4 & pep5
    #求并集
    union = inter1234 | inter1235 | inter1245 | inter1345 | inter2345
    print(f"4id pep: {len(union)}")
    return union
    
def Get5idPep(pep1, pep2, pep3, pep4, pep5):
    """得到5个引擎鉴定的肽段"""
    inter12345 = pep1 & pep2 & pep3 & pep4 & pep5
    print(f"5id pep: {len(inter12345)}")
    return inter12345
    
# def Draw(set1, set2, set3, set4, set5):
#     # 创建子图对象和绘图的布局
#     fig, ax = plt.subplots(nrows=2, ncols=5)
#     fig.set_size_inches(18, 8)
#     # 定义集合和标签的对应关系
#     set_labels = {
#         (0, 0): (set1[1], set2[1]),
#         (0, 1): (set1[1], set3[1]),
#         (0, 2): (set1[1], set4[1]),
#         (0, 3): (set1[1], set5[1]),
#         (0, 4): (set2[1], set3[1]),
#         (1, 0): (set2[1], set4[1]),
#         (1, 1): (set2[1], set5[1]),
#         (1, 2): (set3[1], set4[1]),
#         (1, 3): (set3[1], set5[1]),
#         (1, 4): (set4[1], set5[1])
#     }
#     tmp_sets = [set1, set2, set3, set4, set5]
#     row = 0
#     col = 0
#     for i in range(5):
#         for j in range(i+1, 5):
#             current_labels = set_labels[(row, col)]
#             venn_diagram = venn2([tmp_sets[i][0], tmp_sets[j][0]], set_labels=current_labels, ax=ax[row][col])
#             ax[row][col].set_title('{} vs. {}'.format(*current_labels), fontsize=16)
#             for label in venn_diagram.set_labels:
#                 label.set_fontsize(14)  # 设置集合名字体大小
#             for label in ax[row][col].texts:  # 遍历每个子图中的文本
#                 label.set_fontsize(14)  # 设置文本字体大小
#             col += 1
#             if col == 5:
#                 row += 1
#                 col = 0
#     # 调整子图之间的间距
#     plt.subplots_adjust(hspace=0.4, wspace=0.25)
#     # 保存高分辨率图片
#     plt.savefig('venn_diagrams.png', format='png', dpi=300, bbox_inches='tight')  # 保存为高分辨率的图片
#     # 显示图形
#     plt.show()

def DrawCol():
    """输入是一个4 * 5的矩阵"""
    # Ecoli, yeast, human, 拟南芥
    data = [[6892, 1428, 1468, 2818, 11680], [19982, 4153, 5300, 11306, 37173], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
    name = ['Dong-E.coli', 'Xu-Yeast', 'Gygi-Human', 'Arabidopsis thaliana']
    # 设置子图和柱子的数量
    num_subplots = 4
    num_bars = 5
    # 创建子图和 ax 对象
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
    # 设置使用 seaborn 的色卡
    colors = sns.color_palette('Set2')
    # 遍历子图，画柱状图
    for i in range(num_subplots):
        row = i // 2
        col = i % 2
        ax = axes[row][col]
        # 设置 x 轴标签和标题
        # ax.set_xlabel('X Label')
        ax.set_ylabel('#Peptides')
        ax.set_title(name[i])
        # 画柱状图
        sns.barplot(x=list(range(1, num_bars + 1)), y=data[i], ax=ax, palette=colors)
    # 调整子图之间的间距
    plt.tight_layout()
    # 保存高分辨率图片
    plt.savefig('venn_diagrams.png', format='png', dpi=300, bbox_inches='tight')  # 保存为高分辨率的图片
    # 显示图形
    plt.show()    
    
if __name__ == "__main__":
    GodelInitialize()
    comet_path = r"E:\wkf\Ecoli\credible_pep\Comet\res\result.txtfilter"
    msgf_path = r"E:\wkf\Ecoli\credible_pep\MSGF+\tsv\all_result.tsvfilter"
    maxquant_path = r"E:\wkf\Ecoli\credible_pep\MaxQuant\combined\txt\msms.txt"
    pFind_path = r"E:\wkf\Ecoli\credible_pep\pFind\pFind-Filtered.spectra"
    msfragger_path = r"E:\wkf\Ecoli\credible_pep\MSFragger\peptide.tsv"
    outpath = r"E:\wkf\Ecoli\credible_pep\credible_pep.txt"
    comet_pep, comet_name, comet_code = CometPep(comet_path)
    msgf_pep, msgf_name, msgf_code = MSGFPep(msgf_path)
    maxquant_pep, maxquant_name, maxquant_code = MaxquantPep(maxquant_path)
    pfind_pep, pfind_name, pfind_code = pFindPep(pFind_path)
    msfragger_pep, msfragger_name, msfragger_code = MSFraggerPep(msfragger_path)
    cre_code = GetCredCodes(comet_code, msgf_code, maxquant_code, pfind_code, msfragger_code) # 可信肽段的编码
    
    # 画韦恩图
    # Draw([comet_pep, comet_name], [msgf_pep, msgf_name], [maxquant_pep, maxquant_name], [pfind_pep, pfind_name], [msfragger_pep, msfragger_name])
    # Get1idPep(comet_pep, msgf_pep, maxquant_pep, pfind_pep, msfragger_pep)
    # 得到可信肽段
    GetCredPep(comet_pep, msgf_pep, maxquant_pep, pfind_pep, msfragger_pep, outpath, cre_code)
    # Get3idPep(comet_pep, msgf_pep, maxquant_pep, pfind_pep, msfragger_pep)
    # Get4idPep(comet_pep, msgf_pep, maxquant_pep, pfind_pep, msfragger_pep)
    # Get5idPep(comet_pep, msgf_pep, maxquant_pep, pfind_pep, msfragger_pep)
    # DrawCol()