# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
"""脚本功能:考察.protein文件中的所有蛋白质group, 要存储肽段，若肽段都已经在之前存在，
则该group被舍弃。"""
import os

def test(group, peptides):
    """处理蛋白质group，判断该蛋白质group的肽段条件是否符合要求"""
    flag = False
    for line in group:
        segs = line.split('\t')
        if segs[0] == '' and segs[1] == '': #肽段信息
            pep = segs[3] + segs[8]  
            if pep not in peptides:  #若该肽段信息在之前没有出现过
                flag = True
                peptides.add(pep)
    return flag

def get_string(group, num):
    res = ""
    segs = group[0].split('\t')
    segs[0] = str(num)
    for i in range(0, len(segs) - 1):
        res += segs[i] + '\t'
    res += '\n'
    for i in range(1, len(group)):
        res += group[i]
    return res

def get_protein_filter(Filename):
    """处理后缀为.protein的文件，将每个protein group存储依次判断是否被过滤
    返回最终输出字符串"""
    res = ""
    peptides = set()  #存储目前已经出现过的肽段
    with open(Filename) as f:
        f = f.readlines()
    group_now = []
    num = 1
    for i in range(2, len(f)):
        segs = f[i].split('\t')
        if segs[0].isdigit():  #group头
            if len(group_now) != 0 and test(group_now, peptides):
                print num
                res += get_string(group_now, num)  #合法，存入结果字符串
                num += 1
            group_now = []  #开始记录新的group
            group_now.append(f[i])
        elif "--" in segs[0]:
            break
        else:
            group_now.append(f[i])
    return res


if __name__ == "__main__":
    protein = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\R1\\pFind.protein"
    out = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\R1\\pFind_filtered.protein"
    res = ""
    res = get_protein_filter(protein)
    with open(out, 'w') as f:
        f.write(res)
        f.close()

