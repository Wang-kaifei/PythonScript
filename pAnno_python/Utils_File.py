'''
Descripttion: 切割pAnno的结果文件
version: 
Author: kfwang
Date: 2022-03-27 15:50:22
LastEditors: Kaifei
LastEditTime: 2023-01-08 13:05:23
'''
import math

def tail(filename, n=10):
    with open(filename, 'r') as file:
        lines = file.readlines()
        last_lines = lines[-n:]
        for line in last_lines:
            print(line.strip())


if __name__ == "__main__":
    infile = r"G:\kfwang\human\PaceEva\3\trans_res\database.fasta"
    tail(infile, n=10)