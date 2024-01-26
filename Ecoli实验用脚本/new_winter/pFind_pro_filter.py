#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from functools import cmp_to_key

"""根据.spectra文件，对protein打分，增加protein FDR的过滤条件
"""


class PSM:
    def __init__(self, name='', scanno=0, emh=0.0, chg=0, qv=0.0, seq='', tmh=0.0, Mass_Shift = 0.0, raw_score=0.0,
                 final_score=0.0, modstr='', specificity=3, prostr='', sitestr='', labelstr='', target=True):
        self.name = name
        self.scanno = scanno
        self.emh = emh
        self.chg = chg
        self.qv = qv
        self.seq = seq
        self.tmh = tmh
        self.Mass_Shift = Mass_Shift
        self.raw_score = raw_score
        self.final_score = final_score
        self.modstr = modstr
        self.specificity = specificity
        self.prostr = prostr
        self.sitestr = sitestr
        self.labelstr = labelstr
        self.target = target
        self.pro_value = 0
    def change (self):
        """convert psm to string"""
        res = str(self.name) +'\t' + str(self.scanno) + '\t' + str(self.emh) + '\t' + str(self.chg) + '\t' + str(self.qv) + '\t' + str(self.seq) + '\t' + str(self.tmh) + '\t' +str(self.Mass_Shift) + '\t' + str(self.raw_score) + '\t' + str(self.final_score) + '\t' + str(self.modstr) + '\t' + str(self.specificity) + '\t' + str(self.prostr) + '\t' + str(self.sitestr) + '\t' + str(self.labelstr) + '\t' + str(self.target) + '\n'
        return res

   
def read_pfind_spectra(fpath, qv_limit=0.05):
    """read.spectra file"""
    ret = []
    with open(fpath, 'r') as f:
        line = f.readline()
        while True:
            line = f.readline().rstrip()
            if not line:
                break
            segs = line.split('\t')
            if float(segs[4]) > qv_limit:
                break
            ret.append(segs)
    return ret
     
def read_pfind_spectra_as_psms(fpath, qv_limit=0.05):
    """Convert each line in the .spectra file into a PSM object"""
    ret = []
    tmp = read_pfind_spectra(fpath, qv_limit)
    for segs in tmp:
        if len(segs) > 5:
            ret.append(PSM(segs[0], int(segs[1]), float(segs[2]), int(segs[3]), float(segs[4]), segs[5], float(segs[6]), float(segs[7]),
                           float(segs[8]), float(segs[9]), segs[10], int(segs[11]), segs[12], segs[13], segs[14], segs[15] == 'target'))
        else:
            ret.append(PSM(segs[0], int(segs[1]), float(segs[2]), int(segs[3]), float(segs[4])))
    return ret


def _make_key(psm, mod=True, chg=True, label=True):
    return psm.seq + (psm.modstr if mod else '') \
        + (psm.labelstr if label else '') + (str(psm.chg) if chg else '')#seq+mod+label+charge


def recalc_pro_value(psms):
    """Recalculate q-value
    input: psms sorted according to final_score"""
    t = 1e-15  #target count
    d = 0 #decoy count
    pro = set()
    wrong = set()
    fdr = [] #FDR computing process
    for psm in psms:
        proteins = psm.prostr.strip().split('/')
        for i in range(0, len(proteins) - 1):
            #if proteins[i] not in pro and proteins[i] not in wrong:
            if 'REV' == proteins[i][0 : 3]:
                pro.add(proteins[i])
                d += 1
            else:
                wrong.add(proteins[i])
                t += 1
        fdr.append(d / t)
        #print(d / t)
    for i in range(len(psms) - 2, -1, -1):
        if fdr[i] > fdr[i + 1]:
            fdr[i] = fdr[i + 1]
    for i in range(len(psms)):    #Recalculate q-value
        psms[i].pro_value = fdr[i]
    return psms

def get_protein(psms, threshold_pro = 0.01, threshold_pep = 0.01):
    """得到FDR过滤后的蛋白质"""
    res = set()
    wrong = set()
    count = 0
    for psm in psms:
        count += 1
        if psm.pro_value > threshold_pro or psm.qv > threshold_pep:
            print("psm.pro_value: ", psm.pro_value)
            print("psm.qv: ", psm.qv)
            print("psm: ", count)
            break
        proteins = psm.prostr.strip().split('/')
        for i in range(0, len(proteins) - 1):
            if 'REV' == proteins[i][0 : 3]:
                    continue
            res.add(proteins[i])
            if 'sp' == proteins[i][0 : 2]:
                wrong.add(proteins[i])
    print("all:", len(res))
    print("sp", len(wrong))
    return res

def write_PSM(fpath, psms, qv_limit = 0.01):
    title = "File_Name\tScan_No\tExp.MH+\tCharge\tQ-value\tSequence\tCalc.MH+\tMass_Shift(Exp.-Calc.)\tRaw_Score\tFinal_Score\tModification\tSpecificity\tProteins\tPositions\tLabel\tTarget/Decoy\n"
    filterd_spectra_file = os.path.join(fpath, 'new.Filtered.spectra')
    
    with open(filterd_spectra_file, "w") as f: #psms which under the qv limit
        f.write(title)
        for psm in psms:
            if psm.qv > qv_limit:  
                break
            if not psm.target:
                continue
            f.write(psm.change())
        f.close()

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
    specpath = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\Ecoli\\database2\\O2\\pFind.spectra"
    fasta = "Y:\\chihao_dong111\\dong_111_raw\\ecoli-plus-human-reviewed.fasta"
    out = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\pFind_database_pro\\new_0.01_0.01.fasta"
    psms = read_pfind_spectra_as_psms(specpath)
    psms = recalc_pro_value(psms)
    protein_names = get_protein(psms, threshold_pro = 0.01, threshold_pep = 0.01)
    get_sub (protein_names, fasta, out)
    
    
    

