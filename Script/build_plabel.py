#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 11:10:07 2020

@author: kaifeiwang
"""
"""脚本功能：根据.spectra信息生成.plabel文件"""
import os

def depart(spectra_path):
    """将.spectra文件中的PSM按照mgf文件分类"""
    res = {}
    with open(spectra_path) as f:
        f = f.readlines()[1:]
    for line in f:
        name = line.split('\t', 1)[0].split('.')[0]
        if name in res:
            res[name].append(line)
        else:
            res[name] = [line]
    return res 


def build_plabel(infos, output_path, mgf_path):
    """函数功能：根据info中的信息（谱图肽段信息，字符串格式）生成.plabel文件"""
    dict_mod = {'Carbamidomethyl[C]' : 1, 'Gln->pyro-Glu[AnyN-termQ]' : 2, 'Oxidation[M]' : 3, 'Acetyl[ProteinN-term]' : 4}
    #修饰信息
    w_mod = "[Modification]\n1=Carbamidomethyl[C]\n2=Gln->pyro-Glu[AnyN-termQ]\n3=Oxidation[M]\n4=Acetyl[ProteinN-term]\n"
    #.mgf文件路径
    w_filepath = "[FilePath]\nFile_Path=" + mgf_path + infos[0].split('\t', 1)[0].split('.')[0] + "_HCDFT.mgf\n"
    w_xlink = "[xlink]\nxlink=NULL\n"
    #谱图数量
   
    #谱图信息
    w_spectra = ""
    low = 0
    mid = 0
    high = 0
    cnt = 0
    #得到mgf文件路径
    for i in range(0, len(infos)):   #遍历输入信息的每一行
        line = infos[i].split('\t')  #将该行信息分开存储
        if float(line[4]) > 0.05 or high == 2:
            break
        if float(line[4]) < 0.025 and low == 2 or mid == 2 and float(line[4]) < 0.03:
            continue
        try:
            title = "[Spectrum" + str(i + 1) + ']' + '\n'
            name = "name=" + line[0] + '\n'
            #肽段类型、序列、打分、修饰信息
            pep1 = "pep1=0 " + line[5] + " " + line[9] + " "
            per_mods = line[10].split(';')  #分别得到该行的所有修饰
            for i in range(0, len(per_mods) - 1):
                pos = per_mods[i].split(',')
                pep1 += (str(pos[0]) + ',' + str(dict_mod[pos[1]]) + ' ')  #0为修饰位点，1为修饰名
            pep1 += '\n'  #该psm记录完毕
            w_spectra += title + name + pep1
            if float(line[4]) < 0.025:
                low += 1
            elif float(line[4]) < 0.03:
                mid += 1
            else:
                high += 1
            cnt += 1
        except:
            continue
    #生成输出文件
    w_total = "[Total]\ntotal=" + str(cnt) + '\n'
    output_path += infos[0].split('\t', 1)[0].split('.')[0] + ".plabel"
    with open(output_path, 'w') as f:
        f.write(w_filepath)
        f.write(w_mod)
        f.write(w_xlink)
        f.write(w_total)
        f.write(w_spectra)
        f.close()
        

if __name__ == "__main__":
    if (len(sys.argv) != 4):
        print("Error: Wrong param number!")
        exit(0)
            
    spectra = sys.argv[1] #.spectra文件
    output_path = sys.argv[2] #.plabel文件输出路径
    mgf_path = sys.argv[3] #存储谱图文件的文件夹
    
    # spectra = "Y:\\pAnno\\PT1\\pFind_res\\Open\\result\\pFind-Filtered.spectra"
    # output_path = "D:\\pAnno\\pt1\\pFind-Filtered.plabel"
    # mgf_path = "Y:\\pAnno\\PT1\\raw\\"
    dicts = depart(spectra)
    for key in dicts:
        build_plabel(dicts[key], output_path, mgf_path)
    
    














            
            
        