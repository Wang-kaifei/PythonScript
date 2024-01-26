# -*- coding: utf-8 -*-
"""主要功能为判断差集sequence是否能在对方sequence集中找到类似项。类似项判断规则：
1. L、I同构
2. 作为子串存在于对方集合
判断：差集中的每个seq都进行同构体扩展，如果和对方集合有交集则计数++"""

import math

def ChangeStr(input_str, index, char):
    if(index >= len(input_str)):
        print("Error, index out of the length!")
        return input_str
    s1 = list(input_str)
    s1[index] = char
    return ''.join(s1)
    
def GetLI(sequence):
    """由输入序列产生所有的同构体"""
    res = set()
    index = []
    #记录同构体位置并计数
    str_start = sequence  #同构体初始形态，全为L
    for i in range(0, len(sequence)):
        if sequence[i] == 'L' or sequence[i] == 'I':
            str_start = ChangeStr(str_start, i, 'L')
            index.append(i)
    cnt = int(math.pow(2, len(index)))
    #遍历同构体所有的情况  L为0，I为1 
    for i in range(0, cnt):
        str_tmp = str_start
        j = 0
        while i > 0:
            if i % 2 == 1:
                str_tmp = ChangeStr(str_tmp, index[j], 'I')
            j += 1
            i = i // 2
        res.add(str_tmp)
    return res

def GetSeq(input_path):
    """读取sequence"""
    res = set()
    with open(input_path) as f1:
        f11 = f1.readlines()
    for i in range(0, len(f11)):
        res.add(f11[i])
    return res

def Test_LIcnt(all_set, sub_path):
    """all_set为对方集合（被减）
    sub_path为差集文件
    函数功能为遍历整个差集文件，对每条sequence都进行同构体扩展，如果扩展后的结果与对方集合存在交集，则计数++"""
    cnt = 0
    with open(sub_path) as f1:
        f11 = f1.readlines()
    for i in range(0, len(f11)):
        LI_expand = GetLI(f11[i])
        if len(LI_expand & all_set) > 0:
            cnt += 1
    print(cnt)

def Test_Subcnt(all_set, sub_path):
    """all_set为对方集合（被减）
    sub_path为差集文件
    函数功能为遍历整个差集文件，对每条sequence都遍历all_set，若找到某all_set中的序列包含该差集序列，则计数++"""
    cnt1 = 0
    cnt2 = 0
    with open(sub_path) as f1:
        f11 = f1.readlines()
    for i in range(0, len(f11)):
        for seq in all_set:
            if f11[i] in seq:
                cnt1 += 1
            if seq in f11[i]:
                cnt2 += 1
            #因为提前取了差集，所以不会有完全相同的情况
    print(cnt1)
    print(cnt2)
    print(cnt1 + cnt2)

if __name__ == "__main__":
    gape_path = "Z:\\kfwang\\pAnno\\pt1\\All_GSSP_only.txt" #被减数集合
    res_gape_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\open_gape.txt" #差集
    Test_LIcnt(GetSeq(gape_path), res_gape_path)
    #Test_Subcnt(GetSeq(gape_path), res_gape_path)


