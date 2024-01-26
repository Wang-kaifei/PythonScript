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
        '''mods = line_part[10].split(';')  #提取肽段修饰
        for k in range(0, len(mods) - 1):   #判断修饰是否存在非常规
            if mods[k].split(',')[1] not in list_mod:  #若该肽段存在非常规修饰，则把它过滤
                flag = 1
                break    
        if flag:    #若出现非常规修饰过滤该行
            continue'''
        protein = line_part[12].strip().split('/')  #提取蛋白质名称
        for j in range(0, len(protein) - 1):
            if protein[j][0:3] != 'REV':
                proteins.add(protein[j])
        if float(line_part[4]) > threshold :  #超过阈值，之后的spectra都不予考虑
            break
    return proteins

def test(protein_names, protein):
    for pro in protein_names:
        if pro in protein:
            return True, pro
    return False, pro

def get_sub (protein_names, fasta_path, outfile):
    """函数功能：从路径为fasta_path的.fasta文件中提取与protein_names集合中相同的
    蛋白质信息，写入新的fasta文件。
    """
    used = set()
    res = ""
    count = 0
    print(len(protein_names))
    with open(fasta_path) as file:
        file = file.readlines()
    for i in range(0, len(file)):  #遍历每一行
        if file[i][0 : 1] == ">":  #如果该行是蛋白质名
            #protein_name = file[i].split(' ')[0][1:] #提取与get_protein函数匹配的形式
            protein_name = file[i]
            #if protein_name in protein_names or protein_name.split('|')[0] in protein_names:  #如果该蛋白质名在候选集合中
            a, b = test(protein_names, protein_name)
            if a:
                #print("right")
                if b in used:
                    print(b)
                used.add(b)
                res += file[i]   #存储蛋白质信息
                count += 1
                i += 1
                while i < len(file) and file[i][0] != ">":
                    res += file[i]
                    i += 1
                i -= 1
    print(count)
    print(len(used))
    with open(outfile, "w") as f:
        f.write(res)


if __name__ == "__main__":
    spectra_paths = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\pFind\\pFind.spectra"
    thresholds = [0.0005, 0.001, 0.01, 0.05]
    output_paths =["C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\database\\pFind_1\\0.0005.fasta","C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\database\\pFind_1\\0.001.fasta", "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\database\\pFind_1\\0.01.fasta", "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\one_step\\database\\pFind_1\\0.05.fasta"]
    database_path = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\ecoli_plus_human_con.fasta"  #["Y:\\chihao_dong111\\dong_111_raw\\ecoli_MG1655_20151014_NC_000913_con.fasta", "Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"]
    for i in range(0, len(thresholds)):
        print("Build the ", i, "new database of ", thresholds[i], " FDR rate!")
        protein_names = get_protein(thresholds[i], spectra_paths)
        get_sub(protein_names, database_path, output_paths[i])  #0 -> 0; 1 -> 1





            
