#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""MS-GF+结果合并并过滤"""

class PSM:
    def __init__(self, Spec_file = "", Spec_ID = "", ScanNum = 0, Scan_time = 0.0, FragMethod = "", Precursor = 0.0, IsotopeError = 0, Precursor_Error = 0.0, charge = 0, Peptide = "", Protein = "", DeNovoScore = 0, MSGFScore = 0, SpecEvalue = 0.0, Evalue = 0.0, target = "", q_value = 0.0):
        self.Spec_file = Spec_file
        self.Spec_ID = Spec_ID
        self.ScanNum = ScanNum
        self.Scan_time = Scan_time
        self.FragMethod = FragMethod
        self.Precursor = Precursor
        self.IsotopeError = IsotopeError
        self.Precursor_Error = Precursor_Error
        self.charge = charge
        self.Peptide = Peptide
        self.Protein = Protein
        self.DeNovoScore = DeNovoScore
        self.MSGFScore = MSGFScore
        self.SpecEvalue = SpecEvalue
        self.Evalue = Evalue
        self.target = False if target.lower() == "false" else True
        self.q_value = q_value
        
    def to_string(self):
        return str(self.Spec_file) + '\t' + str(self.Spec_ID) + '\t' + str(self.ScanNum) + '\t' + str(self.Scan_time) + '\t' + str(self.FragMethod) + '\t' + str(self.Precursor) + '\t' + str(self.IsotopeError) + '\t' + str(self.Precursor_Error) + '\t' + str(self.charge) + '\t' + str(self.Peptide) + '\t' + str(self.Protein) + '\t' + str(self.DeNovoScore) + '\t' + str(self.MSGFScore) + '\t' + str(self.SpecEvalue) + '\t' + str(self.Evalue) + '\t' + str(self.target) + '\t' + str(self.q_value) + '\n'
        #return '\t'.join(('%s' % item for item in self.__dict__.values())) + '\n'

def add_flag(line):
    """判断该PSM是否为target"""
    proteins = line.strip().split('\t')[10].split(',')
    for protein in proteins:
        if len(protein.strip()) == 0:
            continue
        if protein.strip()[0 : 3].lower() != "rev":
            return line.strip() + '\tTrue'
    return line.strip() + '\tFalse'

def string_to_PSM(line):
    segs = line.strip().split('\t')
    return PSM(segs[0], segs[1], int(segs[2]), float(segs[3]), segs[4], float(segs[5]), int(segs[6]), float(segs[7]), int(segs[8]), segs[9], segs[10], int(segs[11]), int(segs[12]), float(segs[13]), float(segs[14]), segs[-1])

def read_MSGF_lines(MSGF_file):
    """处理单个MSGF结果，得到target标记
    返回由每行PSM格式化后组成的list""" 
    print(MSGF_file.split('\\')[-1])
    with open(MSGF_file, 'r') as f:
        ff = f.readlines()
    res = []
    cnt = 0
    for i in range(1, len(ff)):
        line = add_flag(ff[i]) #添加target标记
        if float(ff[i].split('\t')[15]) < 0.01:
            cnt += 1
        res.append(string_to_PSM(line))
    print(cnt)
    return res

def get_MSGF_result(folder_path, out_file, add1):
    """将文件夹中的结果文件都合并（按照e-value从小到大排序）
    重新打分，得到q_value"""
    title = "#Spec_file\tSpec_ID\tScanNum\tScan_time\tFragMethod\tPrecursor\tIsotopeError\tPrecursor_Error\tcharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecEvalue\tEvalue\ttarget\tq_value\n"
    import os
    res = []
    file_names = os.listdir(folder_path)
    for file_ in file_names:
        if file_.split('.')[-1] == 'tsv':
            res.extend(read_MSGF_lines(folder_path + file_))
    res.sort(key=lambda x: x.Evalue) # res存放按e-value排序后的psm
    fdr = [] # float
    t = 0
    d = 1.00001 if add1 else 0 #注意此处，有加一修正
    for gg in res:
        if gg.target:
            t += 1
        else:
            d += 1
        fdr.append(0 if t == 0 else float(d) / t)
    for j in range(len(res) - 2, -1, -1):
        if fdr[j] > fdr[j + 1]:
            fdr[j] = fdr[j + 1]
        res[j].q_value = fdr[j]
    res[len(res) - 1].q_value = fdr[len(res) - 1]  #写入q_value
    cnt = 0
    for j in range(0, len(res)):
        if fdr[j] < 0.01:
            cnt += 1
    print(cnt)
    with open(folder_path + out_file, 'w') as f:
        f.write(title)
        for line in res:
            f.write(line.to_string())

    with open(folder_path + out_file +  "filter", 'w') as f:
        f.write(title)
        for line in res:
            if line.q_value >= 0.01:
                break
            f.write(line.to_string())
    print(len(res))



folder_path = "C:\\Users\\kfwang\\results\\griffin\\two-step\\pFind\\0.1\\MS-GF+\\tsv\\"
out_file = "\\all_result.tsv"
get_MSGF_result(folder_path, out_file, False)
        





