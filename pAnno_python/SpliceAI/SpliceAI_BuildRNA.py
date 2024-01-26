'''
Descripttion: 
version: 
Author: sueRimn
Date: 2022-09-16 15:09:19
LastEditors: sueRimn
LastEditTime: 2022-10-09 19:11:12
'''

"""根据预测结果生成成熟的RNA转录本"""
from platform import node
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np


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
        
def ReadOldName(res_file):
    """读取已跑完结果的染色体名"""
    names = []
    with open(res_file) as f:
        lines = f.readlines()
    for line in lines:
        if line[0] == '>':
            names.append(line.strip().rsplit('_', 1)[0])
    return names

def ReadChrSeq(dna_file):
    """读取DNA文件中的所有正链
    Args:
        dna_file (_str_): dna文件路径
    """
    seqs_dic = {} # key = name; value = sequence
    with open(dna_file) as f:
        lines = f.readlines()
    seq = ""
    name = ""
    cnt = 0
    for line in lines:
        if line[0] != '>':
            seq += line.strip()
        else:
            if seq != "":
                seqs_dic[name] = seq
                seq = ""
            name = line.strip() + f"_{cnt}"
            cnt += 1
    if seq != "":
        seqs_dic[name] = seq
    return seqs_dic

def FilterSite(pre_res):
    """将模型预测结果过滤，剔除普通位点

    Args:
        pre_res (_type_): 模型预测结果
    """
    nodes_list = []
    for i in range(0, len(pre_res)):
        if pre_res[i][0] > pre_res[i][1] and pre_res[i][0] > pre_res[i][2]: # 不存储普通位点
            continue
        if pre_res[i][1] > pre_res[i][2]: # 供体位点
            nodes_list.append(Site(i, 0, pre_res[i]))
        else: # 受体位点
            nodes_list.append(Site(i, 1, pre_res[i]))
    return nodes_list
    
def SinglePre(seq, models, context = 10000):
    """对seq单链的预测"""
    x = one_hot_encode('N'*(context//2) + seq + 'N'*(context//2))[None, :] # 这个在长序列前后补充N，相当于一种判断的padding，保证测试位点的有效区域都能有10k的感受野
    y = np.mean([models[m].predict(x, verbose=1) for m in range(1)], axis=0) # 会生成一个1 * n * 3维的矩阵  第三维依次为普通位点、受体、供体的得分
    nodes_list = FilterSite(y[0]) # 返回预测的所有剪接位点
    return nodes_list

def PreSite(seqs_dic, splice_file, old_name):
    """对列表中的每一个序列预测其碱基位点是剪接位点的可能性

    先输出每个位点的得分情况（整体，n个三维数组）
    根据三维数组将位点归类
    输出归类后的位点信息（acceptor和donor） 包含：位置、碱基、对应项分数
    
    Args:
        seqs_dic (name:seq): 存储染色体序列的dict
    """
    f_filter = open(splice_file, 'w')
    context = 10000
    paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 2)) # 本来是5个模型做集成的，现在只跑一个
    models = [load_model(resource_filename('spliceai', x)) for x in paths]
    for name, sequence in seqs_dic.items():
        print(name)
        if name in old_name:
            continue
        f_filter.write(name + "_pos\n")
        nodes_list_pos = SinglePre(sequence, models) # 正链
        for node in nodes_list_pos:
            f_filter.write(str(node))
        f_filter.write(name + "_neg\n")
        nodes_list_neg = SinglePre(sequence[::-1], models) # 负链
        for node in nodes_list_neg:
            f_filter.write(str(node))
    f_filter.close()


if __name__ == "__main__":
    dna_file = "/home/kfwang/GCF_000001405.40_GRCh38.p14_genomic_.fna"
    splice_file = "/home/kfwang/splice.res"
    old_res = "/home/kfwang/old.res"
    seqs_dic = ReadChrSeq(dna_file)
    print(f"Read end: {len(seqs_dic)}!")
    old_name = ReadOldName(old_res)
    PreSite(seqs_dic, splice_file, old_name)
    