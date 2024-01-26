#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""comet结果合并并过滤"""

class PSM:
    def __init__(self, scan = 0, num = 0, charge = 0, exp_neutral_mass = 0.0, calc_neutral_mass = 0.0, e_value = 0.0, xcorr = 0.0, delta_cn = 0.0, sp_score = 0.0, ions_matched = 0, ions_total = 0, plain_peptide = "", modified_peptide = "", prev_aa = "", next_aa = "", protein = "", protein_count = 0, modifications = "", retention_time_sec = 0.0, sp_rank = 0, target = "", file_name = "", q_value = 0):
        self.scan = scan
        self.num = num
        self.charge = charge
        self.exp_neutral_mass = exp_neutral_mass
        self.calc_neutral_mass = calc_neutral_mass
        self.e_value = e_value
        self.xcorr = xcorr
        self.delta_cn = delta_cn
        self.sp_score = sp_score
        self.ions_matched = ions_matched
        self.ions_total = ions_total
        self.plain_peptide = plain_peptide
        self.modified_peptide = modified_peptide
        self.prev_aa = prev_aa
        self.next_aa = next_aa
        self.protein = protein
        self.protein_count = protein_count
        self.modifications = modifications
        self.retention_time_sec = retention_time_sec
        self.sp_rank = sp_rank
        self.target = True if target == 'target' else False
        self.file_name = file_name
        self.q_value = q_value
    def to_string(self):
        res = [str(self.scan), str(self.num), str(self.charge), str(self.exp_neutral_mass), str(self.calc_neutral_mass), str(self.e_value), str(self.xcorr), str(self.delta_cn), str(self.sp_score), str(self.ions_matched), str(self.ions_total), str(self.plain_peptide), str(self.modified_peptide), str(self.prev_aa), str(self.next_aa), str(self.protein), str(self.protein_count), str(self.modifications), str(self.retention_time_sec), str(self.sp_rank), str(self.target), str(self.file_name), str(self.q_value)]
        return '\t'.join(res) + '\n'
        #return '\t'.join(('%s' % item for item in self.__dict__.values())) + '\n'

def find_rank1(lines):
    """选择rank1，并判断是否为decoy"""
    mins = float(lines[0].split('\t')[5])
    res = 0
    for i in range(1, len(lines)):
        if float(lines[i].split('\t')[5]) < mins:
            mins = float(lines[i].split('\t')[5])
            res = i
    proteins = lines[res].split('\t')[15].split(',')
    for protein in proteins:
        if protein.strip()[0 : 5] != "DECOY":
            return lines[res].strip() + '\ttarget'
    return lines[res].strip() + '\tdecoy'

def string_to_PSM(line):
    segs = line.strip().split('\t')
    return PSM(int(segs[0]), int(segs[1]), int(segs[2]), float(segs[3]), float(segs[4]), float(segs[5]), float(segs[6]), float(segs[7]), float(segs[8]), int(segs[9]), int(segs[10]), segs[11], segs[12], segs[13], segs[14], segs[15], int(segs[16]), segs[17], float(segs[18]), int(segs[19]), segs[20], segs[21])

def read_comet_lines(comet_file):
    """处理单个comet结果，得到target标记和删除非rank1
    返回由每行组成的list""" 
    import os.path
    filename = os.path.basename(comet_file).split('.')[0]
    with open(comet_file, 'r') as f:
        ff = f.readlines()
    res = []
    tmp = [ff[2]]  #存储每个scan的全部匹配，不断更新
    print(filename)
    scan_no_last = int(ff[2].strip().split('\t')[0])  #存储上一个scan号
    i = 3
    while(i < len(ff)):
        segs  = ff[i].strip().split('\t')
        #print(segs[17])
        scan_no_now = int(segs[0]) #当前考察的scan号
        if scan_no_last != scan_no_now: #tmp存完
            line = find_rank1(tmp) #只存rank1
            try:
                res.append(string_to_PSM(line + '\t' + filename))
            except:
                print(line)
            tmp = [ff[i]]
            scan_no_last = scan_no_now
        else: #tmp没有存完
            tmp.append(ff[i])
        i += 1
    line = find_rank1(tmp) #存储最后一组
    res.append(string_to_PSM(line + '\t' + filename))
    return res

def get_commet_result(folder_path, out_file, add1 = False):
    """将文件夹中的结果文件都合并（按照e-value从小到大排序）"""
    title = "scan\tnum\tcharge\texp_neutral_mass\tcalc_neutral_mass\te-value\txcorr\tdelta_cn\tsp_score\tions_matched\tions_total\tplain_peptide\tmodified_peptide\tprev_aa\tnext_aa\tprotein\tprotein_count\tmodifications\tretention_time_sec\tsp_rank\ttarget/decoy\tFile\tq_value\n"
    import os
    res = []
    file_names = os.listdir(folder_path)
    for file_ in file_names:
        if file_.split('.')[-1] == 'txt':
            res.extend(read_comet_lines(folder_path + file_))
    res.sort(key=lambda x: x.e_value)
    #重打分步骤
    fdr = []
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


folder_path = "C:\\Users\\kfwang\\results\\griffin\\two-step\\pFind\\0.1\\comet\\result\\"
out_file = "\\result.txt"
get_commet_result(folder_path, out_file)
        





