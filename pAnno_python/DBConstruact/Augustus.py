'''
Descripttion: Augustus的预测结果转化成fasta文件 
version: 
Author: Kaifei
Date: 2024-01-23 21:43:35
LastEditors: Kaifei
LastEditTime: 2024-01-24 22:39:21
'''
# -*- coding: utf-8 -*-

def ExtractPro(file_path, out_path):
    # 打开文件
    with open(file_path, 'r') as file:
        # 逐行读取文件内容
        lines = file.readlines()
    fopen = open(out_path, 'w')
    cnt = 0
    i = 0
    seqs = ""
    while i < len(lines):
        if lines[i].startswith("# protein sequence = ["):
            if seqs != "":
                fopen.write(f">AugustusPro_{cnt}\n")
                fopen.write(seqs[:-1] + "\n")
            cnt += 1
            seqs = lines[i].strip().split("[")[1]
            if ']' in lines[i]: # 如果这个蛋白质只有一行
                i += 1
                continue
            i += 1
            while i < len(lines) and ']' not in lines[i]:
                seqs += lines[i].strip().split(' ')[1]
                i += 1
            seqs += lines[i].strip().split(' ')[1] # 添加最后一行
        i += 1
    #写出最后一个
    fopen.write(f">AugustusPro_{cnt}\n")
    fopen.write(seqs[:-1] + "\n")
    fopen.close()

if __name__ == "__main__":
    inpath = r"E:\FTPtrans\pAnno\DBBuild\Augustus\Ecoli\pre_res"
    outpath = r"E:\FTPtrans\pAnno\DBBuild\Augustus\Ecoli\pre.fasta"
    ExtractPro(inpath, outpath)
