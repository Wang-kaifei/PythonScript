# -*- coding: utf-8 -*-
import os

"""脚本功能：使用处理config文件并重写"""
def build_list(configs):
    #处理每一行元素，最终返回两个list，分别是变量和值
    value = []
    parameter = []
    for config in configs:
        segs = config.split('(')[1].split(')')[0]  #取出括号中的元素
        segs1 = segs.split(',') #将括号中的元素分割
        parameter.append(segs1[1].strip())
        value.append(segs1[2].strip())
    return parameter, value

def write_config(parameter, value, out_file):
    #将提取出的参数写出
    line = ""
    for i in range(0, len(parameter)):
        line += parameter[i] + " = " + value[i] + "\n"
    with open(out_file, 'w') as f:
        f.write(line)

if __name__ == "__main__":
    path = "C:\\Users\\kfwang\\test_wkf\\bumbershoot\\default.cfg"
    out_file = "C:\\Users\\kfwang\\test_wkf\\bumbershoot\\default_new.cfg"
    with open(path) as f:
        configs = f.readlines()
    parameter, value = build_list(configs)
    write_config(parameter, value, out_file)
    
