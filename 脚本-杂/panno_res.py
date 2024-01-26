'''
Descripttion: 
version: 
Author: Kaifei
Date: 2023-03-26 18:08:28
LastEditors: Kaifei
LastEditTime: 2023-03-26 18:11:44
'''
import csv

# 打开输入文件和输出文件
with open(r'C:\Users\pFind\Desktop\pAnno results\Human\group-filtered-merged.panno', 'r', newline='', encoding='utf-8') as infile, \
     open(r'C:\Users\pFind\Desktop\pAnno results\Human\group-filtered-merged-mini.panno', 'w', newline='', encoding='utf-8') as outfile:


    # 逐行处理输入文件
    for line in infile:
        row = line.rstrip('\n').split('\t')  # 去除行末的换行符并按Tab键分割单元格

        # 处理每个单元格
        for i, cell in enumerate(row):
            if len(cell) > 200:
                row[i] = cell[:200] + '$'  # 在截断后的单元格后面添加$标记

        # 将处理后的行写入输出文件
        new_line = '\t'.join(row) + '\n'
        outfile.write(new_line)