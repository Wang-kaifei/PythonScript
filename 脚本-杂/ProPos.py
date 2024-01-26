'''
Descripttion: 找到以肽段修饰在蛋白质中的位置，并以protein name作为主键，输出所有的修饰位点
version: 
Author: kfwang
LastEditTime: 2022-09-08 10:36:27
'''
# -*- coding: utf-8 -*-
# @author Wang kaifei

import re
import sys

def ProcessLine(pros, poses, mods):
    """处理每一行肽段

    Args:
        pros (_type_): 肽段对应到的蛋白质
        poses (_type_): 蛋白质起始位点
        mods (_type_): 修饰
    """
    pro_dict = {} # 返回该行的修饰位点信息
    for i in range(len(pros)):
        if pros[i] == "":
            continue
        for mod in mods:
            mini_segs = mod.split(',')
            if len(mini_segs) == 2: # 合法的修饰
                if pros[i] in pro_dict.keys():
                    pro_dict[pros[i]].add(str(int(poses[i].split(',')[0]) + int(mini_segs[0])) + ',' + mini_segs[1])
                else:
                    pro_dict[pros[i]] = set()
                    pro_dict[pros[i]].add(str(int(poses[i].split(',')[0]) + int(mini_segs[0])) + ',' + mini_segs[1])
            
    return pro_dict

def ProPos(input_file, output_file):
    """在.protien文件的基础上，为每一肽段行最后添加一列，标注修饰的位点
    考虑一下对应到多个蛋白质的情况。

    Args:
        input_file (_type_): 输入文件路径(.protein文件)
        output_file (_type_):
    """
    pro_dict = {} #key: protein name; value: modifications [set]
    with open(input_file) as f:
        lines = f.readlines()
    del lines[0]
    cnt = 0
    for line in lines:
        cnt += 1
        segs = re.split("\t", line)
        if len(segs) > 13:
            pros = segs[12].split('/')
            poses = segs[13].split('/')
            mods = segs[10].split(';')
            mini_dict = ProcessLine(pros, poses, mods)
            for key, value in mini_dict.items():
                if cnt < 10:
                    print(key)
                    print(value)
                if key in pro_dict.keys():
                    pro_dict[key] = pro_dict[key] | value
                else:
                    pro_dict[key] = set()
    ofile = open(output_file, 'w')
    for key, value in pro_dict.items():
        ofile.write(key + '\t' + ';'.join(value) + '\n')
    ofile.close()
            
if __name__ == "__main__":
    if (len(sys.argv) != 3):
        print("Error: Wrong param number!")
        exit(0)      
    input_file = sys.argv[1] # 输入文件
    output_file = sys.argv[2] # 输出文件
    ProPos(input_file, output_file)