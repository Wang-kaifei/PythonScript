# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""定量分析脚本，分析内容为Nan的比例
    对象为，每一个实验结果肽段之间的交集，Nan/交集大小"""

import os

def Unite(setA, setB):
    #求集合A和B的交集
    return setA.intersection(setB)

def Union(setA, setB):
    #求集合A和B的并集
    return setA.union(setB)

def Sub(setA, setB):
    #求差集，在A中但是不在B中的元素
    return setA.difference(setB)

def get_ratio(Filename, peps):
    """函数功能: 提取肽段对应的两个ratio值为Nan的个数
    Filename:.spectra.list文件的上上层路径
    peps:存储肽段的集合
    """
    quantfile = ""
    for dirpath, dirnames, filenames in os.walk(Filename):
        for file in filenames:
            if file == "pQuant_spectra.list_nan":
                quantfile = os.path.join(dirpath, file)
                break
    with open(quantfile) as f1:
        f11 = f1.readlines()
    res = 0  #记录Nan值个数
    for i in range(1, len(f11)):
        segs = f11[i].rstrip().split('\t')
        mods = segs[2].split('|')
        pep_test = segs[1]
        for i in range(1, len(mods)):
            pep_test += mods[i]
        for pep in peps:
            if pep_test == pep:
                res += 1
                break
    return res

def get_pep(Filename):
    """函数功能：提取文件中的肽段信息
    """
    Filename += '\\pFind-Filtered.spectra'
    with open(Filename) as f1:
        f11 = f1.readlines()
    res = set()
    for i in range(1, len(f11)):
        segs = f11[i].rstrip().split('\t')
        mods = segs[10].split(';')
        pep = segs[5]
        for mod in mods:
            pep += mod
        res.add(pep)
    return res 

if __name__ == "__main__":
    #查看一次搜索和二次搜索差集的ratio值
    #一次搜索文件路径
    dir_list1 = ["C:\\test_wkf\\two_step_test\\Ecoli\\database1\\complete1",
                 "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\complete2"]
    #二次搜索0.01文件路径
    dir_list2 = [["C:\\test_wkf\\two_step_test\\Ecoli\\OR_new\\no_filter\\OR1_0.01", "C:\\test_wkf\\two_step_test\\Ecoli\\OR_new\\no_filter\\OR1_0.05", "C:\\test_wkf\\two_step_test\\Ecoli\\OR_new\\no_filter\\OR1_1"],
                 ["C:\\test_wkf\\two_step_test\\Ecoli\\OR_new\\no_filter\\OR2_0.01", "C:\\test_wkf\\two_step_test\\Ecoli\\OR_new\\no_filter\\OR2_0.05", "C:\\test_wkf\\two_step_test\\Ecoli\\OR_new\\no_filter\\OR2_1"]]
    output = "C:\\test_wkf\\two_step_test\\Ecoli\\two_one_sub.txt" #结果输出路径
    res = "Test_Name\tNan_Ratio\n"
    for i in range(0, len(dir_list1)):
        print("start ", i)
        pep_path1 = get_pep(dir_list1[i])
        for j in range(0, len(dir_list2[i])):
            res += dir_list1[i].rsplit('\\', 1)[1] + " VS " + dir_list2[i][j].rsplit('\\', 1)[1] + '\t'
            pep_path2 = get_pep(dir_list2[i][j])
            #求path1-path2的pep集合  
            peps = Sub(pep_path1, pep_path2)
            print(len(peps))
            #根据肽段在path1的.spectra.list文件中得到ratio信息
            res += str(float(get_ratio(dir_list1[i], peps)) / float(len(peps))) + '\n'
            res += dir_list2[i][j].rsplit('\\', 1)[1] + " VS " + dir_list1[i].rsplit('\\', 1)[1] + '\t'
            #求path2-path1的ms2集合  
            peps = Sub(pep_path2, pep_path1)
            print(len(peps))
            #根据肽段在path2的.spectra.list文件中得到ratio信息
            res += str(float(get_ratio(dir_list2[i][j], peps)) / float(len(peps))) + '\n'
    with open(output, 'w') as f:
        f.write(res)
        f.close()



        
'''if __name__ == "__main__":
    #查看一次搜索和二次搜索差集的ratio值
    #二次搜索0.01文件路径
    dir_list = ["C:\\test_wkf\\two_step_test\\Ecoli\\OR_new\\no_filter\\OR1_0.01", "C:\\test_wkf\\two_step_test\\Ecoli\\OR_new\\no_filter\\OR1_0.05",
                 "C:\\test_wkf\\two_step_test\\Ecoli\\OR_new\\no_filter\\OR1_1"]
    output = "C:\\test_wkf\\two_step_test\\Ecoli\\two_two_sub.txt" #结果输出路径
    res = "Test_Name\tNan_Ratio\n"
    for i in range(0, len(dir_list) - 1):   #路径1
        pep_path1 = get_pep(dir_list[i])
        for j in range(i + 1, len(dir_list)):   #路径2
            res += dir_list[i].rsplit('\\', 1)[1] + " VS " + dir_list[j].rsplit('\\', 1)[1] + '\t'
            pep_path2 = get_pep(dir_list[j])
            #求path1-path2的pep集合  
            peps = Sub(pep_path1, pep_path2)
            print(len(peps))
            #根据肽段在path1的.spectra.list文件中得到ratio信息
            res += str(float(get_ratio(dir_list[i], peps)) / float(len(peps))) + '\n'
            res += dir_list[j].rsplit('\\', 1)[1] + " VS " + dir_list[i].rsplit('\\', 1)[1] + '\t'
            #求path2-path1的ms2集合  
            peps = Sub(pep_path2, pep_path1)
            print(len(peps))
            #根据肽段在path2的.spectra.list文件中得到ratio信息
            res += str(float(get_ratio(dir_list[j], peps)) / float(len(peps))) + '\n'
    with open(output, 'w') as f:
        f.write(res)
        f.close()'''
            

            
