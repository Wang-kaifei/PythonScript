'''
Descripttion: 
version: 
Author: sueRimn
Date: 2022-10-08 17:27:03
LastEditors: sueRimn
LastEditTime: 2022-10-08 17:37:44
'''
# -*- coding: utf-8 -*-
# @author Wang kaifei

"""思路：
1. 得到path下的所有文件名，遍历每个文件
2. 单个文件的处理：S开头的行是这个谱图的起始tag，只需要作为flag即可，这一行和下一行都不需要具体的处理；
   之后的行是肽段行，segs[13]是蛋白质描述块，需要修改。分割符是;如果分隔后list的长度==4，则需要修改。
   每个谱图只保留第一个肽段
   
   问题：蛋白质描述不一定是segs[13]，因为可能有很多蛋白质，所有有;的seg都应该被考虑"""
import sys
import os
from tqdm import tqdm

def CleanRes(input_file, output_file):
    print(f"Processing file {input_file}\n")
    with open(input_file, 'r') as f:
        lines = f.readlines()
    f = open(output_file, 'w')
    if_first = True # 是否为第一个排名为1的
    for line in lines:
        segs = line.strip().split('\t')
        if segs[0].isdigit(): # 如果是整数，说明是肽段行
            if int(segs[0]) == 1 and if_first: # top1
                if_first = False
                if float(segs[5]) > 0.05: # 不符合q-value要求的行被过滤了
                    continue
                for i in range(len(segs)):
                    if ';' in segs[i]: # 如果segs[i]是蛋白质描述块
                        mini_segs = segs[i].strip().split(';')
                        if len(mini_segs) == 4: # 旧版pFind结果，需要修改
                            is_decoy = "1" if int(mini_segs[0]) % 2 == 0 else "0"
                            mini_segs.insert(1, is_decoy)
                            segs[i] = ';'.join(mini_segs)
                out_line = '\t'.join(segs) + '\n' # 输出结果
                f.write(out_line)
            else:
                if_first = True
        else: # 如果不是肽段行，则直接输出
            f.write(line)
            if_first = True


if __name__ == "__main__":
    if (len(sys.argv) != 3):
        print("Error: Wrong param number!")
        exit(0)
            
    input_path = sys.argv[1] # 输入文件夹
    output_path = sys.argv[2] # 输出文件夹
    print(input_path)
    print(output_path)
    files= os.listdir(input_path) #得到文件夹下的所有文件名称
    s = []
    for file in files: #遍历文件夹
        if not os.path.isdir(file): #判断是否是文件夹，不是文件夹才打开
            if (file.strip().split('.')[-1] == "res"): # 处理.res文件
                CleanRes(input_path + "\\" + file, output_path + "\\" + file)


