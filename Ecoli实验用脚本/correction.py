#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 12:15:00 2020

@author: kaifeiwang
"""
"""脚本功能：改正错误参数配置"""
import os

def file_name(file_dir):
    """函数功能：提取file_dir文件夹下的所有pQuant.cfg文件名返回"""
    L = []
    for dirpath, dirnames, filenames in os.walk(file_dir):
        for file in filenames:
            if file == 'pQuant.cfg':
                L.append(os.path.join(dirpath, file))
    return L

def alter(filename, old_str, new_str):
    """函数功能：在filename文件中将old_str替换为new_str"""
    file_data = ""
    with open(filename, "r") as f:
        for line in f:
            if old_str in line:
                line = line.replace(old_str, new_str)
            file_data += line
    with open(filename, "w") as f:
        f.write(file_data)


if __name__ == "__main__":
    path = file_name("C:\\test_wkf\\two_step_test\\Ecoli\\after_two_step\\")
    print(path)
    for filename in path:
        alter(filename, "3|none|15N_Labeling|13C_Labeling", "3|none|R:*{N,15N}M:Deamidated[N]{N,15N}M:Gln->pyro-Glu[AnyN-termQ]{N,15N}|R:*{C,13C}M:Deamidated[N]{C,13C}M:Gln->pyro-Glu[AnyN-termQ]{C,13C}")
