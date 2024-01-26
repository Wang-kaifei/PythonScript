# -*- coding: utf-8 -*-

"""脚本功能：将fasta文件随机按比例删除
删除规则为：从.protein文件中随机选取蛋白质，一旦被选中，则整个group中的蛋白质名均被存储到set中
与原始fasta文件取差集输出最终结果

正确性验证：从.protein文件中读取的蛋白名都能以某种规则匹配到fasta

问题：.protein文件是根据pFind的鉴定结果写出的，其实会包含很多target100%文件所未能覆盖的蛋白质，这样就会导致比例的不准确

"""
from math import remainder
from modulefinder import packagePathMap
import random
from tqdm import tqdm

def IsNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

def CMP (a):
    return len(a)

def BuildDic(fasta_path):
    """读取fasta文件，构造name2seq dict 同时存储name2ttname dict"""
    name_seq = {}
    name_ttname = {}
    with open(fasta_path) as file:
        lines = file.readlines()
    name = ""
    seq = ""
    total_name = ""
    for line in lines:
        if line[0] == ">":
            if seq != "":
                name_seq[name] = seq
                name_ttname[name] = total_name
                total_name = line.strip()
            line = line.strip()
            pos1 = line.find(" ")
            pos2 = line.find('|')
            pos3 = line.find('\t')
            if pos3 > 0 and (pos3 < pos1 or pos1 < 0):
                pos1 = pos3
            if pos2 < pos1 and pos2 > 15:
                pos1 = pos2
            if pos1 < 0:
                pos1 = len(line)
            name = line[1:pos1]
            name = name.split('/')[0]
            seq = ""
        else:
            seq += line.strip()
    # 存储最后一个
    if seq != "" and name != "":
        name_seq[name] = seq
        name_ttname[name] = total_name
    print("Build end, all:", len(name_seq), " proteins")
    return name_seq, name_ttname

def ReadGroup(filepath):
    """从自己生成的.protein文件中读取protein group list并按照包含蛋白个数排序"""
    res = []
    with open(filepath) as f:
        lines = f.readlines()
    del lines[0]
    group = [] # 存储当前group中包含的蛋白名
    for line in lines:
        segs = line.split('\t')
        if IsNumber(segs[0]):
            if len(group) != 0:
                res.append(group)
            group = [] # group清空
            if 'REV' != segs[1][0:3]: #填入lead name
                group.append(segs[1].strip())
        elif IsNumber(segs[1]):
            if 'REV' != segs[3].strip()[0:3]:
                group.append(segs[3].strip())
        else:
            print(segs)
    if len(group) != 0: # 存储最后一个group
        res.append(group)
    res = sorted(res, key = len, reverse=True)
    print(f"group length: {len(res)}\n")
    return res


def GetDelPace(group, pace:int, pro_cnt:int):
    res = set() # 记录被删除的蛋白质名
    del_cnt = 0
    remain = len(group) % pace
    for i in range(0, len(group) - remain, pace):
        index = random.randint(i, i + pace - 1)
        res = res.union(set(group[index]))
        if len(res) > (float(pro_cnt) / pace): # 使删除比例更准确
            print(len(res))
            break
    if remain > (pace / 2): #如果剩余量大于可提取的一半
        index = random.randint(len(group) - remain, len(group) - 1)
        res = res.union(set(group[index]))
    print(f"Delete num: {len(res)}")
    return res


def ReadNameDel (filepath, threshold):
    """从pFind.protein文件随机提取蛋白质名，按group提取"""
    proteins = set()
    with open(filepath) as f:
        f = f.readlines()
    print(len(f))
    i = 0
    cnt = 0
    while(i < len(f)):
        if '---' in f[i]:
            break
        segs = f[i].split('\t')
        if IsNumber(segs[0]): # group
            cnt = cnt + 1
            if random.random() <= threshold: #被选中
                if 'REV' != segs[1].strip()[0:3]: #填入lead name
                    proteins.add(segs[1].strip())
                i = i + 1 # 填入sameset、subset
                segs = f[i].split('\t')
                while (IsNumber(segs[3])):
                    if 'REV' != segs[2].strip()[0:3]:
                        proteins.add(segs[2].strip())
                    i = i + 1
                    segs = f[i].split('\t')
            else:
                i = i + 1
        else:
            i = i + 1
    print(".protein cnt: ", cnt)
    print(len(proteins))
    return proteins


def IfDelete(name, del_set):
    for de_name in del_set:
        if de_name in name:
            return True
    return False


def ProCnt(fasta_file):
    cnt = 0
    with open(fasta_file) as file:
        lines = file.readlines()  #以行形式读取fasta文件
    for line in lines:
        if line[0] == ">":
            cnt += 1
    print(f"Pro cnt: {cnt}")
    return cnt


def StorePro (fasta_file, del_set, out_file):
    """"只存储不在del_set中出现的蛋白名并写出"""
    with open(fasta_file) as file:
        lines = file.readlines()  #以行形式读取fasta文件
    cnt_re = 0
    with open(out_file, "w") as f:
        i = 0
        while(i < len(lines)):  #遍历每一行
            if lines[i][0] == ">": #蛋白名
                if not IfDelete(lines[i], del_set): #如果是不删除的蛋白
                    cnt_re = cnt_re + 1
                    f.write(lines[i])
                    i = i + 1
                    while(i < len(lines) and lines[i][0] != '>'):
                        f.write(lines[i])
                        i = i + 1
                else:
                    i = i + 1
                    while(i < len(lines) and lines[i][0] != '>'):
                        i = i + 1
    print("cnt retain: ", cnt_re)


def DelRandom(fasta_path, res_path, out_path, scale):
    """完全随机删除，没有考虑均匀"""
    del_set = ReadNameDel(res_path, scale)
    StorePro(fasta_path, del_set, out_path)


def DelPace(fasta_path, group_path , out_path, scale):
    group = ReadGroup(group_path)
    del_set = GetDelPace(group, scale, ProCnt(fasta_path))
    StorePro(fasta_path, del_set, out_path)

if __name__ == "__main__":
    fasta_path = r"D:\users\kfwang\pAnno_2\human\database\target100%.fasta"
    group_path = r"D:\users\kfwang\pAnno_2\human\database\pFind_Full\mygroup.protein"
    out_path = r"D:\users\kfwang\pAnno_2\human\database\target_group_pace2.fasta"
    pace = 2
    DelPace(fasta_path, group_path, out_path, pace)