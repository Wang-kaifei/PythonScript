'''
Descripttion: 
version: 
Author: Kaifei
Date: 2023-04-11 22:22:38
LastEditors: Kaifei
LastEditTime: 2023-04-26 17:06:34
'''
# -*- coding: utf-8 -*-

"""
划分RNA回贴蛋白的训练样本
提取已注释蛋白 和 第一轮结果回贴RNA的蛋白，计算两个集合间蛋白质的两两归一化得分
以已注释蛋白为准，先根据与某已注释蛋白与其他RNA蛋白的得分分布，判断其是否被召回
1. 未被召回：过滤掉
2. 已被召回：将与该已注释蛋白得分最高的，以及得分超过 90% * max_score的RNA蛋白，判定为正样本
"""

from Bio import Align
from Bio.Align import substitution_matrices
from Bio import SeqIO
import numpy as np
import multiprocessing
from functools import partial
from queue import Queue
from threading import Thread
import concurrent.futures
import editdistance
import heapq
from tqdm import tqdm
from simhash import Simhash
# from tqdm.contrib import multiprocessing

hashbits = 128.0
AATABLE = "ARNDCQEGHILKMFPSTWYV" # 氨基酸字符表
aa_weights = {'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
                  'Q': 146.15, 'E': 147.13, 'G': 75.07, 'H': 155.16, 'I': 131.17,
                  'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
                  'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15}                  
class MaxHeap:
    def __init__(self, capacity):
        self.capacity = capacity
        self.data = []
    
    def push(self, num, index):
        if len(self.data) < self.capacity:
            heapq.heappush(self.data, (num, index))
        elif num > self.data[0][0]:
            heapq.heapreplace(self.data, (num, index))
    
    def get_top(self):
        return sorted(self.data, key=lambda x: x[0], reverse=True)
    
def DeleILL(seqs):
    """删除非氨基酸非法字符"""
    res = []
    for seq in seqs:
        res.append("".join(char for char in str(seq) if char in AATABLE))
    return res

def ReadFastaFile(file_path):
    """
    读取fasta格式的文件，并返回序列列表

    Args:
        file_path: fasta文件路径，字符串类型

    Returns:
        序列列表，包含多个序列，每个序列由protein_name、id和seq组成的元组
    """
    seq_list = []
    with open(file_path, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq_list.append(str(record.seq))      
    print(seq_list[52253])
    return seq_list

def ReadORFs(file_path):
    """
    读取.panno格式的文件，并返回序列列表

    Args:
        file_path: .panno文件路径，字符串类型

    Returns:
        序列列表，包含多个orf的文本信息
    """
    orf_list = []
    seqs = []
    with open(file_path, "r") as f:
        orf_list = f.readlines()
    del orf_list[0]
    for orf in orf_list:
        seqs.append(orf.split('\t')[17])
    print(orf_list[760])
    return orf_list, seqs

def CheckLenDiff(pros, orf_pros, i):
    """
    为第i条已注释蛋白计算长度flag

    Args:
        pros: 已注释蛋白质序列列表
        orf_pros: orf对应的蛋白质列表
        i: 本次考察的已注释的下标

    Returns:
        bool矩阵，如果j处为true，则代表第i个已注释蛋白和第j个orf蛋白长度差异不超过40%
    """
    res = np.zeros((1, len(orf_pros)), dtype=bool)
    len_a = len(pros[i])
    for j in range(len(orf_pros)):
        len_diff = abs(len_a - len(orf_pros[j]))
        diff = len_diff / len_a
        res[0][j] = diff < 0.3
    return res, i

def MolecularWeight(seqs):
    """
    计算list中氨基酸序列的质量
    """
    mass = []
    for seq in seqs:
        seq_weight = 0
        for aa in seq: # 遍历氨基酸序列，计算质量
            if aa not in aa_weights:
                seq_weight += 128
            else:
                seq_weight += aa_weights.get(aa, 0)
        mass.append(seq_weight)
    return mass            

def GetSimHash(seqs):
    # 类似于计算mass，为所有seqs计算simhash存储
    hash_list = []
    for s in tqdm(seqs, ascii=True):
        hash_code = Simhash(s, f=128)
        hash_list.append(hash_code)
    return hash_list

def MoleWeightRatio(pros_mass, orf_pros_mass, pros, orf_pros, n, i):
    """
    为第i条已注释蛋白计算该行的质量flag

    Args:
        pros_mass: 已注释蛋白质序列列表
        orf_pros_mass: orf对应的蛋白质列表
        n: 保留前n名
        i: 本次考察的已注释的下标

    Returns:
        质量相似度矩阵 1 * n
    """    
    res = np.zeros((1, len(orf_pros_mass))) # 质量相似度
    res_len = np.zeros((1, len(orf_pros_mass))) # 长度相似度
    flag = np.zeros((1, len(orf_pros_mass)), dtype=bool)
    rec_cols = [] # 存储被召回的列号
    weight1 = pros_mass[i]
    len1 = len(pros[i])
    for j in range(len(orf_pros_mass)):
        if weight1 == 0:
            res[0][j] = 0
        else:
            res[0][j] = (-abs(weight1 - orf_pros_mass[j])) / weight1 # 得分越大，相似度越高
    for j in range(len(orf_pros)):
        res_len[0][j] = (-abs(len(orf_pros[j]) - len(orf_pros[j]))) / len1 # 得分越大，相似度越高
    # 只保留前n个
    heap = MaxHeap(n) # 容量为n的大顶堆
    for j in range(len(orf_pros_mass)): # 遍历orf蛋白构建候选集
        if (res[0][j] > -0.4 and res_len[0][j] > -0.4): # 只有质量、长度差均不超过40%才纳入
            heap.push(res[0][j] + res_len[0][j], j)
    top_n = heap.get_top()
    for num, index in top_n:
        flag[0][index] = True
        rec_cols.append(index)
    return rec_cols, i

# def BaseScore(pros_mass, orf_pros_mass, pros, orf_pros, n):
#     """
#     最基础的筛选，返回flag矩阵，1为长度差异小于30%

#     Args:
#         pros_mass: 已注释蛋白质序列列表
#         orf_pros_mass: orf对应的蛋白质列表
#         n: 基础打分保留的候选数
        
#     Returns:
#         二维打分矩阵，行表示已注释蛋白，列表示orf对应的蛋白序列
#         (i, j)为pros中第i个蛋白与orfs中第j个蛋白的flag
#         flag为0或1，1为接下来进行粗打分的标记
#     """
#     matrix = np.zeros((len(pros), len(orf_pros)))
#     flag_dic = {} # key: row index; value: recalled cols for this row 
#     MoleWeightRatio_1 = partial(MoleWeightRatio, pros_mass, orf_pros_mass, pros, orf_pros, n)
#     pool = multiprocessing.Pool(processes=50)
#     results = list(tqdm(pool.imap(MoleWeightRatio_1, [i for i in range(len(pros))], chunksize=100), total=len(pros), ascii=True))
#     for pair in results:
#         flag_dic[pair[1]] = pair[0]
#     print(f"Base score: {177 in flag_dic[1686]}")
#     return flag_dic


def HashSimRatio(pros_hash, orf_pros_hash, n, i):
    """
    为第i条已注释蛋白计算该行的质量flag

    Args:
        pros_mass: 已注释蛋白质序列的hash值列表
        orf_pros_mass: orf对应的蛋白质的hash值列表
        n: 基础打分保留的候选数
        i: 本次考察的已注释的下标

    Returns:
        质量相似度矩阵 1 * n
    """    
    res = np.zeros((1, len(orf_pros_hash))) # hash相似度
    rec_cols = [] # 存储被召回的列号
    hash1 = pros_hash[i]
    for j in range(len(orf_pros_hash)):
        distance = hash1.distance(orf_pros_hash[j])
        res[0][j] = 1 - distance / hashbits # 得分越大，相似度越高[0, 
    # 只保留前n个
    heap = MaxHeap(n) # 容量为n的大顶堆
    for j in range(len(orf_pros_hash)): # 遍历orf蛋白构建候选集
        heap.push(res[0][j], j)
    top_n = heap.get_top()
    for num, index in top_n:
        #if num > 0.5:
        rec_cols.append(index)
    return rec_cols, i

def BaseScore(pros_hash, orf_pros_hash, pros, orf_pros, n):
    """
    SimHash打分，存储前n个

    Args:
        pros_mass: 已注释蛋白质序列的hash值列表
        orf_pros_mass: orf对应的蛋白质的hash值列表
        n: 基础打分保留的候选数
        
    Returns:
        二维打分矩阵，行表示已注释蛋白，列表示orf对应的蛋白序列
        (i, j)为pros中第i个蛋白与orfs中第j个蛋白的flag
        flag为0或1，1为接下来进行粗打分的标记
    """
    matrix = np.zeros((len(pros), len(orf_pros)))
    flag_dic = {} # key: row index; value: recalled cols for this row 
    HashSimRatio_1 = partial(HashSimRatio, pros_hash, orf_pros_hash, n)
    pool = multiprocessing.Pool(processes=50)
    results = list(tqdm(pool.imap(HashSimRatio_1, [i for i in range(len(pros))], chunksize=100), total=len(pros), ascii=True))
    for pair in results:
        flag_dic[pair[1]] = pair[0]
    print(f"Base score: {54 in flag_dic[52253]}")
    return flag_dic

def SeqEditDisSim(pros, orf_pros, n, enum_flag):
    """
    为第i条已注释蛋白计算基于编辑距离的相似度

    Args:
        pros: 已注释蛋白质序列列表
        orf_pros: orf对应的蛋白质列表
        n: 粗打分保留的候选数
        i: 本次考察的已注释的下标

    Returns:
        此行的粗打分结果、这次计算对应的已注释蛋白下标
    """
    rec_cols = [] # 存储粗打分被召回的列号
    seq1 = pros[enum_flag[0]]
    heap = MaxHeap(n) # 容量为n的大顶堆
    for col in enum_flag[1]: # 遍历候选列
        # 计算两个序列的编辑距离
        distance = editdistance.eval(seq1, orf_pros[col])
        # 计算相似度
        similarity = 1 - distance / max(len(seq1), len(orf_pros[col]))
        if similarity > 0: # 编辑距离相似度大于0才纳入 
            heap.push(similarity, col)
    top_n = heap.get_top()
    for num, index in top_n:
        rec_cols.append(index)
    return rec_cols, enum_flag[0]

def RoughScore(pros, orf_pros, n, base_flag):
    # TODO: python自身的bug，不能传超2GB的给子进程
    """
    已注释蛋白和ORF对应的蛋白两两打分，返回得分矩阵

    Args:
        pros: 已注释蛋白质序列列表
        orf_pros: orf对应的蛋白质列表
        n: 粗打分保留的候选数
        base_flag: 根据长度粗筛后的flag矩阵
        
    Returns:
        二维打分矩阵，行表示已注释蛋白，列表示orf对应的蛋白序列
        (i, j)为pros中第i个蛋白与orfs中第j个蛋白的flag
        flag为0或1，1为接下来进行细打分的标记
        对每一行，保留前n名候选
    """
    flag_dic = {} # key: row index; value: recalled cols for this row 
    enum_flags = []
    # 使用enumerate()函数遍历矩阵中的每一行
    for key, value in base_flag.items():
        # 将每一行和其对应的召回列存入
        enum_flags.append([key, value])
    SeqEditDisSim_1 = partial(SeqEditDisSim, pros, orf_pros, n)
    # 基于编辑距离的相似度矩阵
    pool = multiprocessing.Pool(processes=50)
    results = list(tqdm(pool.imap(SeqEditDisSim_1, [enum_flag for enum_flag in enum_flags], chunksize=100), total=len(pros), ascii=True))
    for pair in results:
        flag_dic[pair[1]] = pair[0]
    print(f"Rough score: {54 in flag_dic[52253]}")
    return flag_dic

def SequenceAlignment(pros, orf_pros, enum_flag, matrix_name="BLOSUM62", open_gap_score=-10, extend_gap_score=-0.5):
    """
    第i个已注释蛋白和对应的orf蛋白两两打分，返回得分矩阵

    Args:
        pros: 已注释蛋白质序列列表
        orf_pros: orf对应的蛋白质列表
        rough_flag: 粗打分flag
        i: 已注释蛋白下标
        matrix_name: 氨基酸替换矩阵名称，默认为"BLOSUM62"
        open_gap_score: 开启gap的惩罚得分，默认为-10
        extend_gap_score: 继续gap的惩罚得分，默认为-0.5

    Returns:
        第i个已注释蛋白的得分矩阵size = 1 * n
    """
    res = np.zeros((1, len(orf_pros)), dtype=np.float32)
    seq1 = pros[enum_flag[0]]
    # 选择一个氨基酸替换矩阵
    aamatrix = substitution_matrices.load(matrix_name)
    # 创建一个 Align 对象
    aligner = Align.PairwiseAligner()
    # 设置比对参数
    aligner.substitution_matrix = aamatrix
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    for col in enum_flag[1]: # 遍历候选列
        # 进行全局比对
        alignment = aligner.align(seq1, orf_pros[col])
        # 计算得分和归一化得分(因为较长的序列有更多的可能性，所以需要归一化)
        score = alignment.score
        norm_score = score / min(len(seq1), len(orf_pros[col]))
        res[0][col] = norm_score
    return res, enum_flag[0]

def FineScore(pros, orf_pros, rough_flag):
    """
    已注释蛋白和候选ORF对应的蛋白两两打分，返回得分矩阵

    Args:
        pros: 已注释蛋白质序列列表
        orf_pros: orf对应的蛋白质列表
        rough_flag: 粗打分结果，1为候选
        
    Returns:
        二维打分矩阵，行表示已注释蛋白，列表示orf对应的蛋白序列
        (i, j)为pros中第i个蛋白与orfs中第j个蛋白的归一化相似度得分
    """    
    enum_flags = []
    for key, value in rough_flag.items():
        # 将每一行和其对应的召回列存入
        enum_flags.append([key, value])
    matrix = np.zeros((len(pros), len(orf_pros)), dtype=np.float32)
    SequenceAlignment_1 = partial(SequenceAlignment, pros, orf_pros)
    pool = multiprocessing.Pool(processes=50)
    results = list(tqdm(pool.imap(SequenceAlignment_1, [enum_flag for enum_flag in enum_flags], chunksize=100), total=len(pros), ascii=True))
    for pair in results:
        matrix[pair[1]] = pair[0]
    return matrix
        
def PartTrainORF(pros, orfs, matrix, rough_flag):
    """
    根据相似度打分结果，划分训练数据集

    Args:
        pros: 已注释蛋白质序列列表
        orfs: 第一轮肽段回贴RNA得到的orf
        matrix: 二维打分矩阵，(i, j)为pros中第i个蛋白与orfs中第j个蛋白的归一化相似度得分
        pos_file: 正样本输出路径
        neg_file: 负样本输出路径
        
    Returns:
        正负样本划分结果，1为正样本，0为负样本
    """
    flag_res = np.zeros(len(orfs), dtype=bool) # ORF标记，1为正样本，0为负样本
    m, n = matrix.shape
    for i in range(m): # 遍历所有已注释蛋白
        tmp = [] # 对于第i个已注释蛋白，将其与orf的打分抽取，第0列是得分，第1列是orf的下标
        for col in rough_flag[i]:
            tmp.append([matrix[i][col], col])
        tmp = sorted(tmp, key=lambda x: x[0], reverse=True)
        if i == 52253:
            print(f'{tmp[0][1]} : {tmp[0][0]:.4f}\n {tmp[1][1]} : {tmp[1][0]:.4f}')
        if tmp[0][0] <= 2: # 暂时：将最大得分小于1的作为未召回蛋白
            continue
        else: # 如果该已注释蛋白被召回
            max_score = tmp[0][0]
            flag_res[int(tmp[0][1])] = True
            k = 1
            while k < len(tmp) and tmp[k][0] >= max_score * 0.9:
                flag_res[int(tmp[k][1])] = True
                k += 1
    return flag_res

def WriteRes(orfs, flags, pos_file, neg_file):
    title = "Is_rna\tIs_mutation\tIs_splice\tIs_subset\tChr\tID\tFrame\tpep_q_value\tNovel_pep\tOld_pep\tStart_codon\tStart\tEnd\tStrand\tMutation_info\tSplice_info\tAnno_pro\tPro_seq\tGene_seq\tGene_q_value\tTarget/Decoy\n"
    with open(pos_file, 'w') as f1, open(neg_file, 'w') as f2:
        f1.write(title)
        f2.write(title)
        for i in range(len(orfs)):
            if flags[i]:
                f1.write(orfs[i])
            else:
                f2.write(orfs[i])

if __name__ == "__main__":
    multiprocessing.set_start_method('spawn')
    fasta_file = r"G:\kfwang\human\database\target100%.fasta" # 模拟已注释库
    orf_file = r"G:\kfwang\human\PaceEva\3\output\rna_anno_fir-filtered.panno" # 模拟orf
    pos_file = r"G:\kfwang\human\PaceEva\3\output\pos.panno" # 正样本
    neg_file = r"G:\kfwang\human\PaceEva\3\output\neg.panno" # 负样本
    pros = DeleILL(ReadFastaFile(fasta_file))
    orfs, orf_pros = ReadORFs(orf_file)
    orf_pros = DeleILL(orf_pros)
    # pros_mass = MolecularWeight(pros) # 已注释蛋白的质量
    # orf_pros_mass = MolecularWeight(orf_pros) # orf_pros质量
    pros_hash = GetSimHash(pros) # 已注释蛋白hash编码
    orf_pros_hash = GetSimHash(orf_pros) # orf_pros hash编码
    base_flag = BaseScore(pros_hash, orf_pros_hash, pros, orf_pros, 200)
    rough_flag = RoughScore(pros, orf_pros, 20, base_flag)
    fine_flag = FineScore(pros, orf_pros, rough_flag)
    flags = PartTrainORF(pros, orf_pros, fine_flag, rough_flag)
    WriteRes(orfs, flags, pos_file, neg_file)
    
