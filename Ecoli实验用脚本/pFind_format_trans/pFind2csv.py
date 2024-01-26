'''
Descripttion: 
version: 
Author: Kaifei
Date: 2023-09-19 11:48:54
LastEditors: Kaifei
LastEditTime: 2023-09-19 23:05:01
'''
# -*- coding: utf-8 -*-

"""脚本功能：将pFind-Filtered.spectra文件转化成MSGF+的csv格式"""

from tqdm import tqdm
from typing import List
from typing import Dict

class pFind_PSM:
    def __init__(self, File_Name="", Scan_No="", Exp_MHplus="", Charge="", Q_value="", Sequence="", Calc_MHplus="",
                 Mass_Shift="", Raw_Score="", Final_Score="", Modification="", Specificity="", Proteins="",
                 Positions="", Label="", Target="", Miss_Clv_Sites="", Avg_Frag_Mass_Shift="", Others=""):
        self.File_Name = File_Name
        self.Scan_No = Scan_No
        self.Exp_MHplus = Exp_MHplus
        self.Charge = Charge
        self.Q_value = Q_value
        self.Sequence = Sequence
        self.Calc_MHplus = Calc_MHplus
        self.Mass_Shift = Mass_Shift
        self.Raw_Score = Raw_Score
        self.Final_Score = Final_Score
        self.Modification = Modification
        self.Specificity = Specificity
        self.Proteins = Proteins
        self.Positions = Positions
        self.Label = Label
        self.Target = Target
        self.Miss_Clv_Sites = Miss_Clv_Sites
        self.Avg_Frag_Mass_Shift = Avg_Frag_Mass_Shift
        self.Others = Others

    def to_string(self):
        return str(self.File_Name) + '\t' + str(self.Scan_No) + '\t' + str(self.Exp_MHplus) + '\t' + str(
            self.Charge) + '\t' + str(self.Q_value) + '\t' + str(self.Sequence) + '\t' + str(
            self.Calc_MHplus) + '\t' + str(self.Mass_Shift) + '\t' + str(self.Raw_Score) + '\t' + str(
            self.Final_Score) + '\t' + str(self.Modification) + '\t' + str(self.Specificity) + '\t' + str(
            self.Proteins) + '\t' + str(self.Positions) + '\t' + str(self.Label) + '\t' + str(self.Target) + '\t' + str(
            self.Miss_Clv_Sites) + '\t' + str(self.Avg_Frag_Mass_Shift) + '\t' + str(self.Others) + '\n'
        # return '\t'.join(('%s' % item for item in self.__dict__.values())) + '\n'

"""读取pFind结果文件，存成PSM的形式"""
def ReadPSM(pfind_res_path):
    res = []
    with open (pfind_res_path, 'r') as file:
        lines = file.readlines()
    del lines[0]
    for line in lines:
        segs = line.split('\t')
        res.append(pFind_PSM(segs[0], segs[1], segs[2], segs[3], segs[4], segs[5], segs[6], segs[7], segs[8], segs[9], segs[10], segs[11], segs[12], segs[13], segs[14], segs[15], segs[16], segs[17]))
    return res

"""MSGF格式的结果和pFind的scan号是对不上的，其ScanNum是在文件中的index，需要做一次转换
title_dic_list: 存数据集文件名对应表的list，其中key = scan name, value = index
dic_index: key = mgf file name, value = the index for the dict in title_dic_list"""
def GetCsvSpec(title_dic_list, dic_index, pfind_spec, inst_mod):
    # Ecoli-1to1to1-un-C13-N15-1000mM-20150823.48516.48516.2.0.dta
    # Ecoli-1to1to1-un-C13-N15-1000mM-20150823_HCDFT.mgf
    segs = pfind_spec.split('.', 1)
    name = segs[0] + "_" + inst_mod + ".mgf"
    num = title_dic_list[dic_index[name]][segs[1]]
    return name, num

"""将pFind的肽段报告转化成MSGF的格式"""
def GetCsvPep(psm: pFind_PSM, mod2mass: Dict[str, str]):
    ##TODO: 没有添加MSGF肽段酶切前后的那两个氨基酸，如何确定呢？毕竟能回贴到很多蛋白
    res_pep = ""
    flag = [0] * len(psm.Sequence) # 每个位点都初始为没有修饰
    mods = psm.Modification.split(';')
    mods_name = []
    # 对修饰遍历
    for i in range(0, len(mods) - 1):
        segs = mods[i].split(',')
        mods_name.append(segs[1])
        flag[int(segs[0])] = len(mods_name)
    # 遍历sequence的每一个氨基酸
    for i in range(0, len(psm.Sequence)):
        res_pep += psm.Sequence[i]
        if flag[i] != 0: # 如果该氨基酸带有修饰
            res_pep += mod2mass[mods_name[flag[i]]]
    return res_pep
    
"""处理pFind格式的PSM，转换成csv格式的PSM
time_dic_list: 存数据集中谱图保留时间表的list，其中key = spectrum index, value = retain time"""
def PSM2CSV(pfind_psms: List[pFind_PSM], title_dic_list, dic_index, inst_mod, time_dic_list, mod2mass):
    csv_psm = []
    tmp_line = ""
    for psm in pfind_psms:
        name, num = GetCsvSpec(title_dic_list, dic_index, psm.File_Name, inst_mod)
        tmp_line = name + "index=" + str(num) + str(time_dic_list[dic_index[name]][num]) + inst_mod + psm.Exp_MHplus + "0" + psm.Mass_Shift + psm.Charge
        tmp_line += GetCsvPep(psm, mod2mass)
        tmp_line += psm.Proteins[0 : len(psm.Proteins) - 1].replace('/', ',') + "\t0\t0\t0\t0\tTRUE\t" + psm.Q_value
        csv_psm.append(tmp_line)
        tmp_line = ""
    return csv_psm
        
        
        
        
    
        
if __name__ == "__main__":
    in_path = "/Users/kaifeiwang/Desktop/trans/origin-two/msms.txt"
    out_path = "/Users/kaifeiwang/Desktop/trans/origin-two/test.spectra"
    mgf_folder = "/Users/kaifeiwang/Desktop/mgf/"
    source_type = "maxquant"  # comet, msgf, maxquant
    protein = to_pfind(source_type, in_path, out_path)
