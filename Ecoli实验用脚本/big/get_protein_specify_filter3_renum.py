# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
"""脚本功能:考察filtered.protein文件中的所有蛋白质group, 要存储肽段，若之前没有出现过的肽段
都是 带非常规修饰或非特异性酶切，则该group被舍弃。"""
import os

def test(group, peptides, num):
    """处理蛋白质group，判断该蛋白质group的肽段条件是否符合要求"""
    MOD_FILTER = ['Carbamidomethyl[C]', 'Oxidation[M]', 'Acetyl[ProteinN-term]', 'Gln->pyro-Glu[AnyN-term]']
    flag = False
    for line in group:
        segs = line.split('\t')
        if segs[0] == '' and segs[1] == '': #肽段信息
            pep = segs[3] + segs[8]  
            if pep not in peptides:  #若该肽段信息在之前没有出现过
                peptides.add(pep)  #存储肽段
                mods = segs[8].split(';') #得到修饰名称
                if len(mods) == 1:  #若没有携带修饰
                    flag = True
                for i in range(0, len(mods) - 1):
                    if mods[i].split(',')[1] in MOD_FILTER and segs[9] == 3:
                       # print mods[i].split(',')[1]
                        flag = True
                        break
    print num, flag
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
    num2 = 1
    for i in range(0, len(f)):
        segs = f[i].split('\t')
        if segs[0].isdigit():  #group头
            if len(group_now) != 0 and test(group_now, peptides, num2):
                #print num
                res += get_string(group_now, num)  #合法，存入结果字符串
                num += 1
            group_now = []  #开始记录新的group
            group_now.append(f[i])
            num2 += 1
        elif "--" in segs[0]:
            break
        else:
            group_now.append(f[i])
    return res


if __name__ == "__main__":
    protein = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\R1\\pFind_filtered.protein"
    out = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\R1\\pFind_filtered3.protein"
    res = ""
    res = get_protein_filter(protein)
    with open(out, 'w') as f:
        f.write(res)
        f.close()

