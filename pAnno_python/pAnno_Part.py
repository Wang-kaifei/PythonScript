# -*- coding: utf-8 -*-

"""脚本功能：考察pAnno鉴定结果与被删除库内容的一致性
已知：原始全库、被保留库
1. 构造被删除库
2. 考察pAnno的每一个基因对应的蛋白质序列，分为：
   a. pAnno蛋白质完全不存在于被删除库[需要再检查一下是否存在于原始全库，以验证正确性]
   b. 被删除库中能找到与pAnno蛋白质完全一致的序列
   c. 被删除库中能找到pAnno蛋白质的父串
   d. 被删除库中能找到pAnno蛋白质的子串
"""

import numpy as np
from tqdm import tqdm
import sys

"""读取蛋白质序列并存入set"""
def GetProtein (fasta_file):
    pro_set = set()
    with open(fasta_file) as file:
        file = file.readlines()  #以行形式读取fasta文件
        seq = ""
        for i in range(1, len(file)):  #遍历每一行
            if file[i][0] == ">":
                pro_set.add(seq)
                seq = ""
            else:
                seq += file[i].strip()
        #存储最后一个蛋白
        pro_set.add(seq)
    return pro_set

"""获取被删除的蛋白质序列"""
def GetDelProtein (all_data, retain_data):
    all_pro_set = GetProtein(all_data)
    re_pro_set = GetProtein(retain_data)
    de_pro_set = all_pro_set - re_pro_set
    print("ALL :  " + str(len(all_pro_set)) + " Re :  " + str(len(re_pro_set)) + " DEl :  " + str(len(de_pro_set)))
    return de_pro_set, re_pro_set

"""判断蛋白质序列是否为起始子错误提前的情况"""
def TestBig (seq, de_pro_set):
    for de in de_pro_set:
        if de in seq:
            return True
    return False

"""判断蛋白质序列是否为子串"""
def TestSmall (seq, de_pro_set):
    for de in de_pro_set:
        if seq in de:
            return True
    return False

"""判断蛋白质序列是否为子串且终止子能匹配上"""
def TestSmallSameEnd (seq, de_pro_set):
    for de in de_pro_set:
        if len(de) < len(seq):
            continue
        if de[len(de) - len(seq):] == seq:
            return True
    return False

def pAnnoResPart (filename, de_pro_set, not_right_file, right_file, max_qvalue):
    """将pAnno报告的新基因分组并输出被判定为不正确的行"""
    nr_file = open(not_right_file, "w") # 错误基因文件
    r_file = open(right_file, "w") # 正确基因文件
    with open(filename) as file: # 读取pAnno输出文件
        lines = file.readlines()
    nr_file.write(lines[0])
    r_file.write(lines[0])
    cnt_true = 0
    cnt_false = 0
    cnt_report = 0
    true_seq_set = set() # 正确蛋白质集合
    del lines[0]
    for line in tqdm(lines, ascii=True): # 遍历每一个基因
        segs = line.strip().split('\t')
        seq = segs[17].strip() #蛋白质序列
        if len(segs) == 22 and float(segs[19]) < 0.5:
            break
        # if segs[10] == "ZZZ":
        #     continue
        cnt_report += 1
        if segs[16] == ".": # 如果没有注释信息，则被直接判定为错误发现
            nr_file.write(line)
            cnt_false += 1
            continue
        if TestSmallSameEnd (seq, de_pro_set): # 该基因被判定为正确
            r_file.write(line)
            cnt_true += 1
            if not TestSmallSameEnd (seq, true_seq_set): # 如果蛋白质不是已存储蛋白质的同终止子串，才被判定为新召回蛋白。
                true_seq_set.add(seq)
        else:
            nr_file.write(line)
            cnt_false += 1
    nr_file.close()
    r_file.close()
    recall = len(true_seq_set) / len(de_pro_set)
    precision =  cnt_true / cnt_report
    print(f"true identified num {cnt_true}")
    print(f"Recall: {len(true_seq_set) / len(de_pro_set)}") #召回率
    print(f"recall denominator  {len(de_pro_set)}")
    print(f"Pre: {precision}") #准确率
    print(f"Pre denominator  {cnt_report}")
    return recall, precision, cnt_true


def TestSplice(end, anno_info:str):
    segs = anno_info.strip().split('\t')
    if len(segs) < 4:
        return False
    anno_start = int(segs[1].split("E")[0])
    anno_end = int(segs[2].split("S")[0])
    anno_strand = segs[3]
    if end == 87438:
        print(f"{anno_start}  {anno_end}  {anno_strand}")
    if anno_strand == '+' and anno_end < end:
        return True
    if anno_strand == '-' and anno_start > end:
        return True
    return False


def NotMatchPos(not_right_file, all_raw_pro_set):
    """考察不能匹配的新基因，是否能被原始已注释库覆盖"""
    match_set = set()
    cnt_match = 0
    small_set = set()
    cnt_small = 0
    big_set = set()
    cnt_big = 0
    cnt_splice = 0
    with open(not_right_file) as file:
        lines = file.readlines()
    for i in range(1, len(lines)):
        segs = lines[i].strip().split('\t')
        seq = segs[16].strip()

        if seq in all_raw_pro_set:
            match_set.add(seq)
            cnt_match = cnt_match + 1
            #print("match  ", lines[i])
        elif TestSmall (seq, all_raw_pro_set):
            small_set.add(seq)
            cnt_small = cnt_small + 1
            #print("small  ", lines[i])
        elif TestSplice(int(segs[11]), segs[15]):
            cnt_splice = cnt_splice + 1
        elif TestBig (seq, big_set):
            big_set.add(seq)
            cnt_big = cnt_big + 1
            #print("big  ", lines[i])
    print(f"Not Recall proteins set:  match: {len(match_set)}\nbig: {len(big_set)}\nsmall: {len(small_set)}\n")
    print(f"Not Recall proteins cnt:  match: {cnt_match}\nbig: {cnt_big}\nsmall: {cnt_small}\n")
    print(f"Not Recall CNT splice : {cnt_splice}")


def GetProLess(length, anno_file):
    """考察.anno文件中，长度 < length的行所占比例"""
    cnt_less = 0
    with open(anno_file) as file:
        lines = file.readlines()
    for i in range(1, len(lines)):
        segs = lines[i].split('\t')
        if len(segs[16].strip()) < length:
            cnt_less = cnt_less + 1
    print(f"Less {length} cnt : {cnt_less}  Pro : {cnt_less / (len(lines) - 1)}")


def WriteDelFasta(de_pro_set, out_del_fasta):
    i = 0
    with open(out_del_fasta, "w") as f:
        for pro in de_pro_set:
            f.write(">" + str(i) + "\n")
            i += 1
            f.write(pro + "\n")
    f.close()


def GetAnnoPro(pAnno_res):
    res = set()
    with open(pAnno_res) as file:
        lines = file.readlines()
    for line in lines:
        res.add(line.split('\t')[16].strip())
    return res


def GetnRecallPro(de_pro_set, pAnno_res, out_nrecall_fasta):
    anno_pro = GetAnnoPro(pAnno_res)
    i = 0
    with open(out_nrecall_fasta, "w") as f:
        for pro in de_pro_set:
            if pro in anno_pro or TestBig (pro, anno_pro): # 如果与anno结果 match 或者能以anno父串的形式存在，则被recall
                continue
            f.write(">" + str(i) + "\n")
            i += 1
            f.write(pro + "\n")
    f.close()

    
if __name__ == "__main__":
    all_data = r"D:\users\kfwang\pAnno\yeast\database\target100%.fasta"
    retain_data = r"D:\users\kfwang\pAnno\yeast\database\target_group_pace2.fasta"
    pAnno_res = r"D:\users\kfwang\pAnno\yeast\test\output\dna_anno-filtered.panno"
    not_right_file = r"D:\users\kfwang\pAnno\yeast\test\output\old_nr.panno"
    right_file = r"D:\users\kfwang\pAnno\yeast\test\output\old_r.panno"
    de_pro_set, re_pro_set = GetDelProtein (all_data, retain_data)
    #de_pro_set = GetProtein(retain_data)
    recalls = []
    precisions = []
    true_identified_nums = []
    q_value = 1
    recall, precision, true_identified_num = pAnnoResPart (pAnno_res, de_pro_set, not_right_file, right_file, q_value)
    recalls.append(recall)
    precisions.append(precision)
    true_identified_nums.append(true_identified_num)
    print(recalls)
    print(precisions)
    print(true_identified_nums)