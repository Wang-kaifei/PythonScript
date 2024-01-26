'''
Descripttion: 
version: 
Author: sueRimn
Date: 2022-09-16 15:09:19
LastEditors: Kaifei
LastEditTime: 2024-01-26 11:35:46
'''

"""需要改成：输入多个染色体序列【包括正+反】"""
from platform import node
from typing import Sequence
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np
import math

class Site:
    def __init__(self, pos, type, score):
        self.pos = pos # 位点在染色体链上的index，从0开始
        self.type = type # 位点类型
        self.score = score # 位点分数（三个分数）
    def __str__(self) -> str:
        scores =  ""
        for s in self.score:
            scores += '\t' + ("%.3f" % s) # 保留三位小数
        return str(self.pos) + '\t' + str(self.type) + scores + '\n'

class SingleNode: 
    """单链表的结点，由于内存限制，使用链表存储预测结果"""
    def __init__(self, item):
        self.item = item
        self.next = None

def ReadOldName(res_file):
    """读取已跑完结果的染色体名"""
    names = []
    with open(res_file) as f:
        lines = f.readlines()
    for line in lines:
        if line[0] == '>':
            names.append(line.strip()) # 也存储了正负
    return names

def CutSeq(sequence, max_length):
    """由于内存不够需要切分序列，不能存在 > 1e8长度的序列，以_cut1为标记"""
    cnt = math.ceil(len(sequence) / max_length) # 切分数量
    cut_pos = len(sequence) // cnt # 切分位置标记
    seqs = [sequence[0 : cut_pos + 5000]] # 最左端段
    for i in range(1, cnt - 1): # 遍历中间的切分
        l = cut_pos - 5000 # 左端延伸5k感受野
        r = min(len(sequence), cut_pos + cut_pos + 5000) 
        seqs.append(sequence[l : r])
        cut_pos += cut_pos
    seqs.append(sequence[cut_pos - 5000:])
    return seqs

def StoreSeq(seqs_dic, seq, name):
    """存储DNA序列"""
    if ">NC_000006.12" in name or ">NC_000011.10" in name:
        return
    if len(seq) > 1e8: # 需要做切分
        print(f"Cut the sequence: {name}")
        tmp_seqs = CutSeq(seq, 1e8)
        for i in range(len(tmp_seqs)):
            seqs_dic[name + f"_cut{i}"] = tmp_seqs[i]
    else:
        seqs_dic[name] = seq

def ReadChrSeq(dna_file):
    """读取DNA文件中的所有正链
    Args:
        dna_file (_str_): dna文件路径
    """
    seqs_dic = {} # key = name; value = sequence
    with open(dna_file) as f:
        lines = f.readlines()
    cnt = 0
    seq = ""
    name = ""
    for line in lines:
        if line[0] != '>':
            seq += line.strip()
        else:
            if seq != "":
                StoreSeq(seqs_dic, seq, name)
            name = line.strip() + f"_{cnt}"
            cnt += 1
            seq = ""
    if seq != "":
        StoreSeq(seqs_dic, seq, name)
    return seqs_dic

def FilterSite(pre_res):
    """将模型预测结果过滤，剔除普通位点

    Args:
        pre_res (_type_): 模型预测结果
    """
    prehead = SingleNode(None)
    p = prehead
    for i in range(0, len(pre_res)):
        if pre_res[i][0] > pre_res[i][1] and pre_res[i][0] > pre_res[i][2]: # 不存储普通位点
            continue
        node = SingleNode(str(Site(i, 0, pre_res[i]))) if pre_res[i][1] > pre_res[i][2] else SingleNode(str(Site(i, 1, pre_res[i])))
        p.next = node
        p = node
    return prehead.next

def SinglePre(name, seq, models, f_filter, context = 10000):
    """对seq单链的预测， 并写出结果"""
    f_filter.write(name + "\n")
    print(f"sequence length: {len(seq)}")
    x = one_hot_encode('N'*(context//2) + seq + 'N'*(context//2))[None, :] # 这个在长序列前后补充N，相当于一种判断的padding，保证测试位点的有效区域都能有10k的感受野
    print("Encoding end!")
    y = models[0].predict(x, verbose=1)
    print("Pre end!")
    # y = np.mean([models[m].predict(x, verbose=1) for m in range(1)], axis=0) # 会生成一个1 * n * 3维的矩阵  第三维依次为普通位点、受体、供体的得分
    nodes = FilterSite(y[0]) # 预测的所有剪接位点
    while nodes is not None:
        f_filter.write(nodes.item)
        nodes = nodes.next
    return
 

def PreSite(seqs_dic, splice_file, old_name):
    """对列表中的每一个序列预测其碱基位点是剪接位点的可能性

    先输出每个位点的得分情况（整体，n个三维数组）
    根据三维数组将位点归类
    输出归类后的位点信息（acceptor和donor） 包含：位置、碱基、对应项分数
    
    Args:
        seqs_dic (name:seq): 存储染色体序列的dict
    """
    f_filter = open(splice_file, 'w', 1)
    paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 2)) # 本来是5个模型做集成的，现在只跑一个
    models = [load_model(resource_filename('spliceai', x)) for x in paths]
    for name, sequence in seqs_dic.items():
        print(name)
        name_pos = name + "_pos"
        name_neg = name + "_neg"
        if name_pos not in old_name:
            print(name_pos)
            SinglePre(name_pos, sequence, models, f_filter) # 正链
        if name_neg not in old_name:
            print(name_neg)
            SinglePre(name_neg, sequence[::-1], models, f_filter) # 负链
        print("end111")
    print("end222")
    f_filter.close()


if __name__ == "__main__":
    dna_file = "/home/kfwang/GCF_000001405.40_GRCh38.p14_genomic_.fna"
    splice_file = "/home/kfwang/splice.res"
    old_res = "/home/kfwang/old.res"
    old_name = ReadOldName(old_res)
    seqs_dic = ReadChrSeq(dna_file)
    print(f"Read end: {len(seqs_dic)}!")
    PreSite(seqs_dic, splice_file, old_name)
    