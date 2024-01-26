'''
Descripttion: Balrog预测结果转化成蛋白库
version: 
Author: Kaifei
Date: 2024-01-23 21:43:35
LastEditors: Kaifei
LastEditTime: 2024-01-26 21:28:52
'''
# -*- coding: utf-8 -*-

"""
Balrog输出结果是gff文件，转化有两步: 1. 通过gff的position信息找到对应的碱基序列
2. 碱基序列翻译成氨基酸序列
3. 输出到fasta文件
"""
from Bio import Seq

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
                seqs_dic[name] = seq
            name = line.split(' ')[0][1:]
            if name in seqs_dic.keys():
                print(f"Warning: {name} has existed.")
            cnt += 1
            seq = ""
    if seq != "":
        seqs_dic[name] = seq
    print(f"Read {cnt} sequences from {dna_file}")
    return seqs_dic

def ReadGFF(gff_file):
    """读取gff文件
    Args:
        gff_file (_str_): gff文件路径
    """
    with open(gff_file) as f:
        lines = f.readlines()
    cnt = 0
    gff_dic = {} # key = name; value = [start, end, strand]
    for line in lines:
        if line[0] == '#':
            continue
        cnt += 1
        line = line.strip().split('\t')
        name = line[0]
        start = int(line[3])
        end = int(line[4])
        strand = line[6]
        if name not in gff_dic.keys():
            gff_dic[name] = []
        gff_dic[name].append([start, end, strand])
    print(f"Read {cnt} lines from {gff_file}")
    return gff_dic

def GFF2bases(gff_dic, seqs_dic):
    """将gff行信息转化成碱基序列
    """
    for name, cdses in gff_dic.items():
        for cds in cdses: # cds = [start, end, strand]
            start = cds[0]
            end = cds[1]
            cds.append(seqs_dic[name][start-1:end])
    return gff_dic

def TranslateOne(bases, strand):
    """将碱基序列翻译成氨基酸序列返回"""
    dna_seq = Seq.Seq(bases)
    if strand == '-':
        dna_seq = dna_seq.reverse_complement()
    aa_seq = dna_seq.translate()
    if aa_seq[-1] == '*':
        aa_seq = aa_seq[:-1]
    return aa_seq
    
def Translate(gff_dic, out_path):
    """将CDS序列翻译成氨基酸序列，并写出fasta文件
    """
    fopen = open(out_path, 'w')
    for name, cdses in gff_dic.items():
        for cds in cdses:
            start = cds[0]
            end = cds[1]
            strand = cds[2]
            bases = cds[3]
            aa = TranslateOne(bases, strand)
            fopen.write(f">{name}_{start}_{end}_{strand}\n")
            fopen.write(f"{aa}\n")
    fopen.close()

if __name__ == "__main__":
    dna_path = r"E:\FTPtrans\pAnno\DBBuild\Balrog\ecoli\ncbi_dataset\ncbi_dataset\GCA_000005845.2_ASM584v2_genomic.fna"
    gff_path = r"E:\FTPtrans\pAnno\DBBuild\Balrog\ecoli\GCA_000005845.2_ASM584v2_genomic.fna__.gff"
    out_path = r"E:\FTPtrans\pAnno\DBBuild\Balrog\ecoli\Balrog_ecoli.fasta"
    dna_seqs_dic = ReadChrSeq(dna_path) # 碱基序列
    gff_dic = ReadGFF(gff_path) # 读取每个CDS行
    gff_dic = GFF2bases(gff_dic, dna_seqs_dic) # 将碱基序列加入到CDS行中
    Translate(gff_dic, out_path)