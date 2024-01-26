# -*- coding: utf-8 -*-
# @author Wang kaifei

import sys
import os
import re

class Mod:
    def __init__(self, site = "", mod = ""):
        self.site = site
        self.mod = mod

class Mgfstart:
    def __init__(self):
        self.title = ""
        self.scans = ""
        self.charge = ""
        self.mass = ""
    def Output(self, i: int):
        return "BEGIN IONS\n" + "TITLE=" + self.title + "." + str(self.scans) + "." + str(self.scans) + "." + str(self.charge) + ".dta\nCHARGE=" + self.charge + "+\nPEPMASS=" + self.mass + "\n"

class Labelline:
    def __init__(self, name = "", pep1 = "0", seq = "", score = "", mod: list = []):
        self.name= name
        self.pep1 = pep1
        self.seq = seq
        self.score = score
        self.mod = mod
    def Output(self):
        res = "name=" + self.name + "\n" + "pep1=" + self.pep1 + " " + self.seq + " " + self.score + " "
        for mod in self.mod:
            res += mod.site + "," + mod_dic[mod.mod] + " "
        res += "\n"
        return res

spec_num = int(1)
mod_dic = {}
fix_mod = {} # key=AA; value=mod_name 

def GetcolDic(title: str, delimiter: str):
    res = {}
    segs = re.split("[" + delimiter + "]", title)
    for i in range (0, len(segs)):
        res[segs[i]] = i
    return res

def GetLabel(res_file: str):
    f = open(res_file)
    lines = f.readlines()
    labels = []
    title_col_dic = GetcolDic(lines[0], '\t') # key = col name; value = index
    rawfile_c = title_col_dic["Raw file"]
    scan_c = title_col_dic["Scan number"]
    modseq_c = title_col_dic["Modified sequence"]
    seq_c = title_col_dic["Sequence"]
    charge_c = title_col_dic["Charge"]
    score_c = title_col_dic["Score"]
    global fix_mod
    for i in range(1, len(lines)):
        segs = re.split('[\t]', lines[i])
        name = segs[rawfile_c] + "." + segs[scan_c] + "." + segs[scan_c] + "." + segs[charge_c] + ".dta"
        mods = [] #存储本条肽段的全部修饰
        mod_seq = segs[modseq_c]
        j = 0
        cnt_site = 0
        mod = Mod() #临时单个修饰存储
        while(j < len(mod_seq)):
            while(j < len(mod_seq) and mod_seq[j] != '('):
                if (mod_seq[j] >= 'A' and mod_seq[j] <= 'Z'):
                    cnt_site += 1
                    if mod_seq[j] in fix_mod.keys(): #如果该氨基酸发生了固定修饰
                        for v in fix_mod[mod_seq[j]]:
                            mod.site = str(cnt_site)
                            mod.mod += v
                            mods.append(mod)
                            mod = Mod()
                j += 1
                continue
            if (j == len(mod_seq)): #已无字符
                break
            #找到左括号
            mod.site = str(cnt_site)
            j += 1
            while(mod_seq[j] != ')'): #存储修饰信息直到右括号
                mod.mod += mod_seq[j]
                j += 1
            mod.mod += ")"
            mods.append(mod) #存储当前修饰，并为下一次修饰的存储做准备
            mod = Mod()
            j += 1
        label = Labelline(name = name, seq = segs[seq_c], score = segs[score_c], mod = mods)
        labels.append(label)
    return labels

def ApltoMgf(apl_file: str, out_file: str):
    mgfstart = Mgfstart()
    f = open(apl_file)
    lines = f.readlines()
    global spec_num
    with open(out_file, "a") as file:
        for line in lines:
            if("start" in line or line[0] == 'f'):
                continue
            if(line[0] == 'm'):
                mgfstart.mass = line.strip().split('=')[1]
            elif(line[0] == 'c'):
                mgfstart.charge = line.strip().split('=')[1]
            elif(line[0] == 'h'):
                segs = line.strip().split()
                mgfstart.title = segs[1]
                mgfstart.scans = segs[3]
                file.write(mgfstart.Output(spec_num))
                spec_num += 1 
                mgfstart = Mgfstart()
            elif("end" in line):
                file.write("END IONS\n")
            else:
                file.write(line)
    f.close()

def WriteHeader(mgf_file: str, cnt: int, out_file: str):
    res = "[FilePath]\nFile_Path=" + mgf_file + "\n" + "[Modification]\n"
    for key in mod_dic:
        res += mod_dic[key] + "=" + key + "\n"
    res += "[xlink]\nxlink=NULL\n[Total]\ntotal=" + str(cnt) + "\n"
    with open(out_file, "w") as file:
        file.write(res)

def OutputPlabel(labels: list, out_file: str):
    with open(out_file, "a") as file:
        i = 1
        for label in labels:
            file.write("[Spectrum" + str(i) + "]\n" + label.Output())
            i += 1

def ReadMod(sum_file: str):
    f = open(sum_file)
    lines = f.readlines()
    title_col_dic = GetcolDic(lines[0], '\t')
    segs = re.split('[\t]', lines[1])
    mods = re.split(';', segs[title_col_dic["Fixed modifications"]]) #此时存的是fix mod
    global fix_mod
    for mod in mods:
        aas = re.split('[()]', mod)
        for aa in aas[1]: #遍历所有发生该修饰的氨基酸
            if aa in fix_mod:
                fix_mod[aa].append(mod)
            else:
                fix_mod[aa] = [mod]
    mods += re.split(';', segs[title_col_dic["Variable modifications"]]) #添加可变修饰
    global mod_dic
    i = 1
    for mod in mods:
        mod_dic[mod] = str(i)
        i += 1
    return mods

if __name__ == "__main__":
    if (len(sys.argv) != 4):
        print("Error: Wrong param number!")
        exit(0)
            
    path = sys.argv[1] #apl文件夹目录
    mqres_path = sys.argv[2] #maxquant搜索结果路径 msms.txt文件
    mqsum_path = sys.argv[3] #maxquant搜索结果路径 summary.txt文件

    #生成.mgf文件
    mgf_file = path + "\\all.mgf"
    with open(mgf_file, "w") as file:
        file.write("MASS=Monoisotopic\n")
    files= os.listdir(path) #得到文件夹下的所有文件名称
    s = []
    for file in files: #遍历文件夹
        if not os.path.isdir(file): #判断是否是文件夹，不是文件夹才打开
            if (file.strip().split('.')[-1] == "apl"):
                ApltoMgf(path+"\\"+file, mgf_file)

    #生成.plabel文件
    plabel_file = path + "\\mq.plabel"
    ReadMod(mqsum_path)
    labels = GetLabel(mqres_path)
    WriteHeader(mgf_file, len(labels), plabel_file)
    OutputPlabel(labels, plabel_file)
