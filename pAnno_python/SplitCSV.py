'''
Descripttion: 切割pAnno的结果文件
version: 
Author: kfwang
Date: 2022-03-27 15:50:22
LastEditors: Kaifei
LastEditTime: 2023-01-08 13:05:23
'''
import math

def GetPep (deal_file, out_path, part_num):
    """将大文件读入并划分为小文件"""
    with open(deal_file, 'r') as f:
        lines = f.readlines()
    title = lines[0]
    del lines[0]
    cnt_line = math.ceil(len(lines) / part_num) # 每个小文件包含行数
    now_cnt = 0
    cnt_file = 1
    fo = open(out_path + f"\\test_panno_{cnt_file}.panno", "w")
    fo.write(title)
    for line in lines:
        if now_cnt > cnt_line: # 写满一个小文件，打开下一个文件
            fo.close()
            cnt_file += 1
            fo = open(out_path + f"\\test_panno_{cnt_file}.panno", "w")
            fo.write(title)
            fo.write(line)
            now_cnt = 1
        else:
            fo.write(line)
            now_cnt += 1

def GetDNAMap (deal_file, out_path):
    """提取DNA回贴结果写出"""
    with open(deal_file, 'r') as f:
        lines = f.readlines()
    title = lines[0]
    del lines[0]
    fo = open(out_path, "w")
    fo.write(title)
    for line in lines:
        segs = line.split('\t', 1)
        if segs[0] == "0":
            fo.write(line)

if __name__ == "__main__":
    infile = r"D:\pAnno_out\test_anno_old.panno"
    out_path = r"D:\pAnno_out\DNA.panno"
    GetDNAMap (infile, out_path)