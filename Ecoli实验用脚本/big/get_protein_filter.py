# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
"""脚本功能:考察.protein文件中的所有蛋白质group，看一下每个group中的支持肽段，如果这些肽段
都是有非常规修饰的，则把该group删除掉"""
import os

def test(group):
    """处理蛋白质group，判断该蛋白质group是否符合要求"""
    mod = ['Carbamidomethyl[C]', 'Oxidation[M]', 'Acetyl[ProteinN-term]', 'Gln->pyro-Glu[AnyN-term]']
    flag = False
    for i in range(0, len(group)):
        segs = group[i].split('\t')
        if segs[0] == '' and segs[1] == '':
            try:
                mod_temp = segs[8].split(';')
                if len(mod_temp) == 1:
                    return True
                flag = False
                for i in range(0, len(mod_temp) - 1):
                    if mod_temp[i].split(',')[-1] in mod:
                        flag = True
            except IndexError:
                print group[i - 1]
    return flag

def get_string(group):
    res = ""
    for line in group:
        res += line
    return res

    
def get_protein_filter(Filename):
    """处理后缀为.protein的文件，将每个protein group存储依次判断是否被过滤
    返回最终输出字符串"""
    res = ""
    with open(Filename) as f1:
        f = f1.readlines()
    group_now = []
    num = 0
    for i in range(2, len(f)):
        segs = f[i].split('\t')
        if segs[0] != '':  #group头
            if len(group_now) != 0 and test(group_now):
                res += get_string(group_now)  #合法，存入结果字符串
            print num
            num += 1
            if num > 5000:
                break
            group_now = []  #开始记录新的group
            group_now.append(f[i])
        else:
            group_now.append(f[i])
    return res

if __name__ == "__main__":
    protein = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\C\\pFind.protein"
    out = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\big\\C\\filtered.protein"
    res = ""
    res = get_protein_filter(protein)
    with open(out, 'w') as f:
        f.write(res)
        f.close()
