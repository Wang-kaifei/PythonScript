# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""脚本功能：从pFind.spectra文件中提取蛋白，然后再到原fasta数据库中查找，并放入新的数据库
   输入：卡的FDR值，pFind.spectra、原fasta文件路径、输出新fasta的路径
"""
import os

def get_protein (threshold, spec_path):
    """函数功能：从路径为spec_path的.spectra文件中提取蛋白质名称放入set返回。
    .spectra文件特点：第5列为q-value,第13列为该肽段对应的蛋白质信息；
    q-value值小于threshold才予以考虑；
    由于一个肽段可能对应很多蛋白质，第13列内容以'/'为分隔符代表多个蛋白质；
    返回蛋白质名的集合。
    """
    list_mod = ['Carboxymethyl[C]', 'Gln->pyro-Glu[AnyN-termQ]', 'Oxidation[M]', 'Acetyl[ProteinN-term]']   #存储常规修饰
    with open(spec_path) as file:
        file = file.readlines()  #以行形式读取spectra文件
    proteins = set()
    for i in range(1, len(file)):  #遍历每一行
        flag = 0  #是否出现非常规修饰的标记
        line_part = file[i].split('\t', 13)
        mods = line_part[10].split(';')  #提取肽段修饰
        for k in range(0, len(mods) - 1):   #判断修饰是否存在非常规
            if mods[k].split(',')[1] not in list_mod:  #若该肽段存在非常规修饰，则把它过滤
                flag = 1
                break    
        if flag:    #若出现非常规修饰过滤该行
            continue
        protein = line_part[12].split('/')  #提取蛋白质名称
        for j in range(0, len(protein) - 1):
            proteins.add(protein[j])
        if float(line_part[4]) > threshold :  #超过阈值，之后的spectra都不予考虑
            break
    return proteins

def filt_protein(filename, proteins):
    """函数功能：过滤支持度小的subset"""
    res = []
    with open(filename) as file:
        file = file.readlines()
    for i in range(2, len(file)):
        list1 = file[i].split('\t', 5)
        if len(list1) > 4 and list1[2] in proteins:
            if list1[4] == '1\n':
                if list1[2] in res:
                    res.remove(list1[2])
                else:
                    res.append(list1[2])
    print len(res)
    res = set(res)
    return proteins.difference(res)

def get_sub (protein_names, fasta_path):
    """函数功能：从路径为fasta_path的.fasta文件中提取与protein_names集合中相同的
    蛋白质信息，写入新的fasta文件。
    .fasta文件的特点：
    每一行只有一列;
    代表蛋白质名的行开头为'3E'。第一个空格之前的部分（去除3E)，可与get_protein函数返回值匹配;
    实现：
    遍历fasta_path指向的文件，对于每一个蛋白质名，都搜索proteins。
    """
    res = ""
    with open(fasta_path) as file:
        file = file.readlines()
    for i in range(0, len(file)):  #遍历每一行
        if file[i][0 : 1] == ">":  #如果该行是蛋白质名
            protein_name = file[i].split(' ', 2)[0][1:] #提取与get_protein函数匹配的形式
            if protein_name in protein_names:  #如果该蛋白质名在候选集合中
                res += file[i]   #存储蛋白质信息
                i += 1
                while i < len(file) and file[i][0 : 1] != ">":
                    res += file[i]
                    i += 1
                i -= 1
    return res

def write_info(infos, filename):
    """将infos写入filename文件"""
    with open(filename, "w") as f:
        for info in infos:
            f.write(info)

"""protein_names = get_protein(0.05, "C:\\test_wkf\\two_step_test\\limit1\\pFind.spectra")
infos = get_sub(protein_names, "Y:\\chihao_dong111\\dong_111_raw\\ecoli_MG1655_20151014_NC_000913_con.fasta")  #0, 2 -> 0; 1, 3 -> 1 
write_info(infos, "C:\\test_wkf\\two_step_test\\limit1\\new_0.05.fasta")"""


if __name__ == "__main__":
    spectra_paths = ["C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2mod4\\pFind.spectra"] #"C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind.spectra"]
    protein_paths = ["C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2mod4\\pFind.protein"]
    thresholds = [0.01, 0.05, 1.00]
    output_paths = [["C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2mod4\\new_0.01.fasta", "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2mod4\\new_0.05.fasta", "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2mod4\\new_1.fasta"]
                     ]
    database_path = ["Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"]#["Y:\\chihao_dong111\\dong_111_raw\\ecoli_MG1655_20151014_NC_000913_con.fasta", "Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"]
    for i in range(0, len(spectra_paths)):
        for j in range(0, len(thresholds)):
            print("Build the ", i, "new database of ", thresholds[j], " FDR rate!")
            protein_names = get_protein(thresholds[j], spectra_paths[i])
            protein_names = filt_protein(protein_paths[i], protein_names)
            infos = get_sub(protein_names, database_path[i])  #0 -> 0; 1 -> 1 
            write_info(infos, output_paths[i][j])




            
