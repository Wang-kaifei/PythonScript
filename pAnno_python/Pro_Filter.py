# -*- coding: utf-8 -*-

"""脚本功能：从pFind-Filtered.spectra中提取蛋白质名，并将不包含在已鉴定结果中的蛋白质从fasta文件删除

构建评测数据库流程:
1. 搜索完整已注释库(pFind_Full结果)
2. Pro_Filter.py pFind_Full结果生成已鉴定数据库(target100%.fasta)
3. ProteinGroup.py 对pFind_Full结果重新做蛋白质group分析, 只有肽段全部被覆盖才作为subset, -> mygroup.protein
4. Pro_Del_group.py以protein group层次生成不完整已注释库
"""

from tqdm import tqdm
import ahocorasick

def BuildACAutomaton(keywords):
    A = ahocorasick.Automaton()
    for idx, keyword in enumerate(keywords):
        A.add_word(keyword, (idx, keyword))
    A.make_automaton()
    return A

def CheckPro(pro_full_name, automaton):
    for end_index, (idx, keyword) in automaton.iter(pro_full_name):
        return True
    return False

def GetProteinname (filename):
    """读取已鉴定蛋白质"""
    pro_set = set()
    with open(filename) as file:
        file = file.readlines()  #以行形式读取文件
    for i in range(1, len(file)):  #遍历每一行
        segs = file[i].split('\t')
        proteins = segs[12].split('/')
        for p in proteins:
            if len(p) < 3:
                continue
            pro_set.add(p)
    print(len(pro_set))
    return pro_set

def StoreIdentified (filename, out_path, trie):
    """"读取fasta文件并只存储已鉴定蛋白质"""
    cnt = 0
    f = open(out_path, "w")
    with open(filename) as file:
        file = file.readlines()  #以行形式读取fasta文件
    seq = ""
    name = file[0].strip()
    for i in tqdm(range(1, len(file)), ascii=True):  #遍历每一行
        if file[i][0] == ">":
            if CheckPro(name, trie) and seq != "":
                f.write(name + "\n")
                f.write(seq + "\n")
                cnt += 1
            seq = ""
            name = file[i].strip()
        else:
            seq += file[i].strip()
    #存储最后一个蛋白
    if CheckPro(name, trie):
        f.write(name + "\n")
        f.write(seq + "\n")
        cnt += 1
    print(f"Stored: {cnt}")
    f.close()

def GetLeadName(in_file):
    """从.protein文件中读取lead + sameset蛋白质名"""
    pro_name = set()
    cnt = 0
    with open(in_file) as file:
        lines = file.readlines()  #以行形式读取fasta文件
    for line in lines:
        if "-----" == line[0 : 5]:
            break
        segs = line.split('\t')
        if segs[0].strip().isdigit() and segs[1][0:3] != "REV":
            pro_name.add(segs[1].strip())
            cnt += 1
        elif segs[1].strip() == "SameSet" and segs[2][0:3] != "REV":
            pro_name.add(segs[2].strip())
    print(f"#Protein Group: {cnt}\n#Protein: {len(pro_name)}")
    return pro_name

if __name__ == "__main__":
    fasta_path = r"D:\users\kfwang\pAnno_2\human\NCBI-download\protein.faa"
    pro_path = r"D:\users\kfwang\pAnno_2\human\database\pFind_Full\pFind.protein"
    pFind_res = r"D:\users\kfwang\pAnno_2\human\database\pFind_Full\pFind-Filtered.spectra"
    out_path = r"D:\users\kfwang\pAnno_2\human\database\target100%.fasta"
    pro_set = GetLeadName(pro_path)
    # pro_set = GetProteinname (pFind_res)
    trie = BuildACAutomaton(pro_set)
    StoreIdentified(fasta_path, out_path, trie)