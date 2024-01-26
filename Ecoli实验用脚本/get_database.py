# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""脚本功能：从pFind.spectra文件中提取蛋白，然后再到原fasta数据库中查找，并放入新的数据库
   输入：卡的FDR值，pFind.spectra、原fasta文件路径、输出新fasta的路径
"""
import os

def get_protein (threshold, spec_path, summary):
    """函数功能：从路径为spec_path的.spectra文件中提取蛋白质名称放入set返回。
    .spectra文件特点：第5列为q-value,第13列为该肽段对应的蛋白质信息；
    q-value值小于threshold才予以考虑；
    由于一个肽段可能对应很多蛋白质，第13列内容以'/'为分隔符代表多个蛋白质；
    返回蛋白质名的集合。
    """
    #得到建库所考虑的修饰
    list_mod = {'Carboxymethyl[C]', 'Gln->pyro-Glu[AnyN-termQ]', 'Oxidation[M]', 'Acetyl[ProteinN-term]'}   #存储常规修饰
    """with open(summary) as file:
        file = file.readlines()  #以行形式读取summary文件
    line_num = 0
    for i in range(0, len(file)):
        if "Modifications:" in file[i]:
            line_num = i + 2    #存储第一个修饰的位置
            break
    for i in range(line_num, line_num + 5):     #提取排名前五的修饰
        list_mod.add(file[i].split('\t', 1)[0])  #存储修饰名称"""
    print list_mod

    
    #提取符合条件肽段对应的蛋白质名
    with open(spec_path) as file:
        file = file.readlines()  #以行形式读取spectra文件
    proteins = set()
    for i in range(1, len(file)):  #遍历每一行
        flag = 0  #是否出现范围外修饰的标记
        line_part = file[i].split('\t', 13)
        mods = line_part[10].split(';')  #提取肽段修饰
        #若存在范围外修饰，则把它过滤
        for k in range(0, len(mods) - 1):   
            if mods[k].split(',')[1] not in list_mod:  
                flag = 1
                break    
        if flag:    
            continue
        
        protein = line_part[12].split('/')  #提取蛋白质名称
        for j in range(0, len(protein) - 1):
            proteins.add(protein[j])
        if float(line_part[4]) > threshold :  #超过阈值，之后的psm都不予考虑
            break
    return proteins

def filt_protein(filename):
    """函数功能：提取支持度小的subset
        将只有一条肽段支持的蛋白质提取"""
    legal = set()
    res = []
    with open(filename) as file:
        file = file.readlines()
    for i in range(2, len(file)):
        list1 = file[i].split('\t', 5)
        if len(list1) > 4:
            if list1[2] not in legal and list1[4] == '1\n':  #不在legal中且在group中只有一条肽段对应
                if list1[2] in res: #若在res中有记录，将其从res删除并存入合法
                    res.remove(list1[2])
                    legal.add(list1[2])
                else:               #若其不再legal中，存入res
                    res.append(list1[2])
    print len(res)
    res = set(res)
    return res

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



if __name__ == "__main__":
    #.spectra文件
    spectra_paths = "C:\\test_wkf\\two_step_test\\yanshou\\R1\\pFind.spectra" #"C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind.spectra"]
    #.protein文件
    protein_paths = "C:\\test_wkf\\two_step_test\\yanshou\\R1\\pFind.protein"
    #.summary文件
    summary_path = "C:\\test_wkf\\two_step_test\\yanshou\\R1\\pFind.summary"
    #分别采用3个阈值
    thresholds = [0.01, 0.05, 1.00]
    #输出路径，对应3个阈值
    output_paths = ["C:\\test_wkf\\two_step_test\\yanshou\\R1\\new_0.01.fasta", "C:\\test_wkf\\two_step_test\\yanshou\\R1\\new_0.05.fasta", "C:\\test_wkf\\two_step_test\\yanshou\\R1\\new_1.fasta"]
    #待搜索数据库
    database_path = "Y:\\chihao_dong111\\dong_111_raw\\ecoli_MG1655_20151014_NC_000913_con.fasta"#["Y:\\chihao_dong111\\dong_111_raw\\ecoli_MG1655_20151014_NC_000913_con.fasta", "Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"]


    dele_pro = filt_protein(protein_paths)    #得到需要过滤的蛋白质
    for i in range(0, len(thresholds)):
        proteins = get_protein(thresholds[i], spectra_paths, summary_path) #根据肽段得到蛋白质名
        proteins = proteins.difference(dele_pro)    #过滤
        res =  get_sub(proteins, database_path)    #新建库
        write_info(res, output_paths[i])    #写出fasta文件




'''if __name__ == "__main__":
    #.spectra文件
    spectra_paths = "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind.spectra" #"C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind.spectra"]
    #.protein文件
    protein_paths = "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind.protein"
    #.summary文件
    summary_path = "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind.summary"
    #分别采用3个阈值
    thresholds = [0.01, 0.05, 1.00]
    #输出路径，对应3个阈值
    output_paths = ["C:\\test_wkf\\two_step_test\\Ecoli\\vs_R\\database2\\new_0.01.fasta", "C:\\test_wkf\\two_step_test\\Ecoli\\vs_R\\database2\\new_0.05.fasta", "C:\\test_wkf\\two_step_test\\Ecoli\\vs_R\\database2\\new_1.fasta"]
    #待搜索数据库
    database_path = "Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"   #["Y:\\chihao_dong111\\dong_111_raw\\ecoli_MG1655_20151014_NC_000913_con.fasta", "Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"]


    dele_pro = filt_protein(protein_paths)    #得到需要过滤的蛋白质
    for i in range(0, len(thresholds)):
        proteins = get_protein(thresholds[i], spectra_paths, summary_path) #根据肽段得到蛋白质名
        proteins = proteins.difference(dele_pro)    #过滤
        res =  get_sub(proteins, database_path)    #新建库
        write_info(res, output_paths[i])    #写出fasta文件
        '''
    



            
