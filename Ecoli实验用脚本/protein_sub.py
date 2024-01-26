# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""脚本功能：从pFind.spectra文件中提取蛋白，求两实验中蛋白鉴定的差集   差集中干扰蛋白的比例
"""
import os

def get_protein (threshold, spec_path):
    """函数功能：从路径为spec_path的.spectra文件中提取蛋白质名称放入set返回。
    .spectra文件特点：第5列为q-value,第13列为该肽段对应的蛋白质信息；
    q-value值小于threshold才予以考虑；
    由于一个肽段可能对应很多蛋白质，第13列内容以'/'为分隔符代表多个蛋白质；
    返回蛋白质名的集合。
    """
    with open(spec_path) as file:
        file = file.readlines()  #以行形式读取spectra文件
    proteins = set()
    for i in range(1, len(file)):  #遍历每一行
        line_part = file[i].split('\t', 13)
        protein = line_part[12].split('/')  #提取蛋白质名称
        for j in range(0, len(protein) - 1):
            if protein[j].find("REV") == -1:
                proteins.add(protein[j])
        if float(line_part[4]) > threshold :  #超过阈值，之后的spectra都不予考虑
            break
    return proteins

def write_info(infos, filename):
    """将infos写入filename文件"""
    with open(filename, "w") as f:
        for info in infos:
            if "gi" in info:
                f.write(info + '\n')

def get(filename):
    with open(filename) as file:
        file = file.readlines()  #以行形式读取spectra文件
    proteins = set()
    for i in range(0, len(file)):
        proteins.add(file[i])
    return proteins

def test_pro(proteins):
    res = 0
    for protein in proteins:
        if "gi" not in protein:
            res += 1
            #print protein
    print float(res) / len(proteins)

if __name__ == "__main__":
    patha = "C:\\test_wkf\\two_step_test\\Ecoli\\OR_new\\no_filter\\OR2_0.01\\result.spectra"
    pathb = "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\complete2\\pFind.spectra"
    proa = get_protein(0.01, patha)
    prob = get_protein(0.01, pathb)
    print len(proa.difference(prob)), len(prob.difference(proa))
    test_pro(proa)
    test_pro(prob)




    

'''if __name__ == "__main__":
    outa = "C:\\test_wkf\\temp.txt"
    outb = "C:\\test_wkf\\temp2.txt"
    a = get(outa)
    b = get(outb)
    print len(a.difference(b))
    print len(b.difference(a))'''




    

            
