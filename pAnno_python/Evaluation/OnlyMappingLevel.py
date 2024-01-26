# -*- coding: utf-8 -*-

"""脚本功能：评测分为三个层次：肽段鉴定、基因组推断和全流程
如果单对基因组推断结果做评测，需要剔除肽段鉴定结果的影响。本脚本处理pAnno报告结果，如果覆盖到的肽段都来自于benchmark结果之外，则删除
"""
from tqdm import tqdm

"""从fasta文件中读取蛋白质序列"""
def ReadProSeq(fasta_path):
    seqs = []
    with open(fasta_path) as file:
        lines = file.readlines()
    seq = ""
    for line in lines:
        if line[0] == ">":
            if seq != "":
                seqs.append(seq)
            seq = ""
        else:
            seq += line.strip()
    if seq != "":
        seqs.append(seq)
    tqdm.write(f"Build end, all: {len(seqs)} protein sequences.")
    return seqs

def test(pep_str, bm_peps):
    peps = pep_str.split(';')
    for pep in peps:
        seq = pep.replace('.', '')
        if len(seq) <= 6:
            print(f"Error: pep length < 5 aa  {seq}")
        if seq in bm_peps:
            return True
    return False

"""过滤pAnno的报告结果"""
def FiterpAnnoRes(in_path, out_path, bm_peps):
    out_lines = []
    in_f = open(in_path)
    line = in_f.readline() # 过滤掉title行
    title = line
    while line:
        line = in_f.readline()
        segs = line.strip().split('\t')
        if len(segs) < 8:
            print(f"length < 8: {line}")
            break
        if test(segs[8], bm_peps):
             out_lines.append(line)
    with open (out_path, 'w') as fout:
        fout.write(title)
        for line in out_lines:
            fout.write(line)
    in_f.close()

"""从.spectra文件中读取肽段"""
def ReadPep(pfind_res_path):
    res = set()
    in_f = open(pfind_res_path)
    line = in_f.readline() # 过滤掉title行
    while line:
        line = in_f.readline()
        segs = line.strip().split('\t')
        if len(segs) < 10:
            break
        res.add(segs[5])
    return res

def GetPep(protein:str):
    """生成蛋白质酶切结果，最大遗漏酶切为2"""
    res = []
    pep = ""
    for aa in protein:
        if aa == 'K' or aa == 'R':
            pep += aa
            res.append(pep)
            pep = ""
        else:
            pep += aa
    if pep != "":
        res.append(pep)
    length = len(res)
    for i in range(1, length): # 遗漏一个酶切位点
        res.append(res[i - 1] + res[i])
    for i in range(2, length): # 遗漏两个酶切位点
        res.append(res[i - 2] + res[i - 1] + res[i])
    return res
    
def Digestion(proteins):
    """将proteins模拟酶切，返回生成的肽段列表"""
    peps = []
    for protein in proteins:
        peps.extend(GetPep(protein))
    return set(peps)

def TestOnePep (pep, pros):
    """判断肽段是否以子串形式存在于pros list"""
    for pro in pros:
        if pep in pro:
            return True
    return False

def TestPepConFast(pep, cleve_peps, pros):
    """判断肽段是否在fasta文件中存在"""
    if pep in cleve_peps: # 先判断肽段是否能匹配上全特异性酶切
        return True
    return TestOnePep (pep, pros) # 再考察以子串的形式匹配

"""过滤pFind鉴定结果，提取序列没有在被保留proteins中出现的PSM输出
作为完美引擎的鉴定结果
"""
def FilterpFindRes(in_path, out_path, fasta_path):
    db_pros = ReadProSeq(fasta_path) # 被保留的蛋白质
    cleve_peps = Digestion(db_pros) # 被保留蛋白质对应的肽段
    res = []
    covered_res = [] # 存储对应到已注释肽段的PSM
    with open(in_path) as file:
        lines = file.readlines()
    title = lines[0]
    del lines[0]
    for line in tqdm(lines):
        seq = line.split('\t')[5]
        if not TestPepConFast(seq, cleve_peps, db_pros):
            res.append(line)
        else:
            covered_res.append(line)
    with open (out_path, 'w') as fout:
        fout.write(title)
        for line in res:
            fout.write(line)
    with open (out_path + "covered", 'w') as fout:
        fout.write(title)
        for line in covered_res:
            fout.write(line)

def OutPepSeq(in_path, out_path):
    """把pfind.spectra文件的sequence输出"""
    with open(in_path) as file:
        lines = file.readlines()
    del lines[0]
    print(len(lines))
    seqs = set()
    for line in lines:
        seq = line.split('\t')[5]
        if (len(seq) > 6):
            seqs.add(seq)
    with open (out_path, 'w') as fout:
        for seq in seqs:
            fout.write(seq + "\n")

def Readpepline(in_path):
    with open(in_path) as file:
        lines = file.readlines()
    pep = set(lines)
    return pep

def Readpeptitle(in_path):
    pep = set()
    with open(in_path) as file:
        lines = file.readlines()
    for line in lines:
        if '>' != line[0]:
            pep.add(line)
    return pep

def TestGSSP(in_path1, in_path2):
    pep1 = Readpepline(in_path1)
    pep2 = Readpeptitle(in_path2)
    print(f"#new GSSP {len(pep1)}, #old GSSP {len(pep2)}.\n #Inter {len(pep1 & pep2)}")

if __name__ == "__main__":
    pAnno_res_path = r"G:\kfwang\human\PaceEva\2\output\dna_anno-filtered-nom.panno"
    out_res_path = r"G:\kfwang\human\PaceEva\2\output\pFind-Filtered.spectracovered"
    bm_pep_path = r"G:\kfwang\human\database\pFind_Full\pFind-Filtered.spectra"
    re_fasta = r"G:\kfwang\human\database\target_group_pace2.fasta"
    pep_path = r"G:\kfwang\human\PaceEva\2-ideal\pepold.txt"
    old_gssp_path = r"G:\kfwang\human\PaceEva\2\output\raw_res\target_GSSP.txt"
    #TestGSSP(pep_path, old_gssp_path)
    #FilterpFindRes(bm_pep_path, out_res_path, re_fasta)
    OutPepSeq(out_res_path, pep_path)
    #peps = ReadPep(bm_pep_path)
    #print(f"# Identified peptides: {len(peps)}")
    