'''
Descripttion: 
version: 
Author: Kaifei
Date: 2024-01-25 17:15:11
LastEditors: Kaifei
LastEditTime: 2024-03-04 20:26:26
'''
# -*- coding: utf-8 -*-
# @author Wang kaifei

"""简化.res文件，只保留排名第一的，以减轻内存的负担
"""
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
                if float(segs[5]) > 0.005: # 不符合q-value要求的行被过滤了
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


