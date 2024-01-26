# -*- coding: utf-8 -*-

"""脚本功能：高等生物结果判定脚本
1. 读取被删除蛋白质名称
2. 遍历gff文件，读取被删除蛋白相关行CDS【目标就是把这些表达片段找回来】
3. 遍历基因报告文件，读取除RNA和剪接外的行，需要有anno_pro信息
4. 第3步与第2步得到的片段比对，交集为TP（其实是第三步的覆盖第二步的）
对于precision: 如果覆盖到的这个gff行是正确的但没有质谱数据支持，也应该判定为对
对于recall：只考察那些有质谱数据覆盖到的gff行
"""
import re
import numpy as np
from tqdm import tqdm
import sys
import ahocorasick
import math
import multiprocessing
from functools import partial
from queue import Queue
from threading import Thread
import concurrent.futures
import time

genes = [] # 存储class Gene的容器
anno_res = []

title_csv = "Is_rna\tIs_mutation\tIs_splice\tIs_subset\tChr\tID\tFrame\tpep_q_value\tNovel_pep\tOld_pep\tStart_codon\tStart\tEnd\tStrand\tMutation_info\tSplice_info\tAnno_pro\tPro_seq\tGene_seq\tGene_q_value\tTarget/Decoy\n"

class AnnoRes:
    """.panno文件的一行"""
    def __init__(self, chr_name, pro_id, pn, min_pos, max_pos, split_info, line_num, legal, aa_seq, total_line, peps, anno_info):
        self.chr_name = chr_name # 染色体名
        self.pro_id = pro_id # 染色体名
        self.pn = pn # 正负链标记
        self.min_pos = min_pos # 左端位置，+start -end
        self.max_pos = max_pos # 右端位置
        self.split_info = split_info
        self.line_num = line_num # 存在panno文件的第几行
        self.legal = legal # 是否为合法报告
        self.aa_seq = aa_seq # 蛋白质序列
        self.total_line = total_line
        self.peps = peps # 对应到此回贴结果的肽段
        self.anno_info = anno_info # 重注释到的gff行
    def __eq__(self, other):
        # return self.chr_name == other.chr_name and self.pn == other.pn and self.min_pos == other.min_pos and self.max_pos == other.max_pos and self.split_info == other.split_info
        return self.pro_id == other.pro_id

class AAseg:
    """已注释基因的一段编码序列"""
    def __init__(self, start, end, length):
        self.start = start # seg左端位置
        self.end = end # seg右端位置
        self.length = length # seg长度 碱基长度
        self.is_recall = False # 此外显子片段是否被召回
        self.peps = [] # 这一段覆盖到的肽段
        
"""Gene from GFF file.
pro_id : segs(start, end, length)"""
class Gene: 
    def __init__(self, chr_name, pro_id, pn):
        self.chr_name = chr_name # 染色体名
        self.pro_id = pro_id # 蛋白质名【与fasta中的蛋白名匹配】
        self.segs = [] # 各段信息，存储aaseg
        self.filtered_segs = [] # 有质谱证据支持的段[待召回的段]
        self.pn = pn # 正负链标记
        self.legal = True # 是否与gff+DNA相对应
        self.is_del = False # 标记该基因是否对应的是被删除的蛋白（应该被召回）

    def Print(self):
        print(f"PN:  {self.pn}  CHR:  {self.chr_name}")
        for seg in self.segs:
            print(f"({seg.start}, {seg.end})")

"""构造name2seq dict"""
def BuildDic(fasta_path):
    name_seq = {}
    with open(fasta_path) as file:
        lines = file.readlines()
    name = ""
    seq = ""
    for line in lines:
        if line[0] == ">":
            if seq != "":
                name_seq[name] = seq
            line = line.strip()
            pos1 = line.find(" ")
            pos2 = line.find('|')
            pos3 = line.find('\t')
            if pos3 > 0 and (pos3 < pos1 or pos1 < 0):
                pos1 = pos3
            if pos2 < pos1 and pos2 > 15:
                pos1 = pos2
            if pos1 < 0:
                pos1 = len(line)
            name = line[1:pos1]
            name = name.split('/')[0]
            seq = ""
        else:
            seq += line.strip()
    # 存储最后一个
    if seq != "" and name != "":
        name_seq[name] = seq
    print("Build end, all:", len(name_seq), " proteins")
    return name_seq

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

"""
从gff文件中读取pro_names中蛋白所对应到的基因
gff_file: gff文件路径
right_pro_name: 报告的应该被判定为正确的蛋白质名
del_pro_name: 删掉的蛋白名（应该被召回的蛋白质名，是right_pro_name的子集）

return:
genes: 删库蛋白对应的基因
seg_num: 删除的基因包含的段数（大于待召回的数量，因为有一些无质谱证据）
"""
def DealGFF(gff_file, right_pro_name, del_pro_name):
    seg_num = 0
    cnt_del_pro = 0
    pre_gene = Gene("", "", "") # 上一个存储的基因
    for line in open(gff_file): 
        segs = line.split('\t')
        if len(segs) < 9 or segs[2] != "CDS":
            continue
        minisegs = segs[8].split("ID=")
        if len(minisegs) == 1: # Can't find Pro ID
            continue
        pro_id = re.split("cds-", re.split(";|\t|\n", minisegs[1])[0])[1]
        if pro_id not in right_pro_name:
            continue
        aaseg = AAseg(int(segs[3]), int(segs[4]), int(segs[4]) - int(segs[3]) + 1)
        if pro_id != pre_gene.pro_id: # 如果是新读入的基因
            if pre_gene.pro_id != "": # 存储前面已经读完的基因
                genes.append(pre_gene)
                if pre_gene.pro_id in del_pro_name:
                    pre_gene.is_del = True
                    cnt_del_pro += 1
            pre_gene = Gene(segs[0], pro_id, segs[6]) # 初始化
        pre_gene.segs.append(aaseg) # 存储段
        seg_num += 1 # 段数自增
    if pre_gene.pro_id != "": # 存储最后一个基因
        genes.append(pre_gene)
        if pre_gene.pro_id in del_pro_name:
            pre_gene.is_del = True
            cnt_del_pro += 1
    print(f"Read GFF file end. And the total seg number is {seg_num}")
    print(f"deleted pro: {cnt_del_pro}")
    return seg_num

"""根据段信息将pro切分，即根据编码段位置信息，把蛋白质分成各自编码段对应的氨基酸序列"""
def PartPro2Seg(gene, pro):
    res = []
    start = 0 # 切分起始位置
    remain = 0 # 上一段留下的碱基数
    for seg in gene.segs:
        aa_len = int(seg.length / 3)
        res.append(pro[start : start + aa_len])
        remain = (seg.length % 3 + remain) % 3
        start += aa_len
        if remain != 0:
            start += 1
    return res

def contains_string(input_str, trie):
    return input_str in trie

"""获取蛋白质上被质谱数据覆盖的位置"""
def GetProCov(pro, trie):
    cov_poses = [] # 存储pro上鉴定到的位置，闭区间
    cov_peps = []
    for end_index, keyword in trie.iter(pro):
        start_index = end_index - len(keyword) + 1
        cov_poses.append((start_index, end_index))
        cov_peps.append(keyword)
    return cov_poses, cov_peps

"""处理待召回基因，对于其中的segs，如果没有除已注释肽段外的质谱证据支持，则删除

pro_names: 待召回蛋白质名
name2seq: key = protin name; value = sequece
trie: 由全部可能的novel peptides组成的自动机

"""
def RefineGeneRec(pro_names, name2seq, trie, pep_uni_dic):
    rec_seg_num = 0
    for gene in genes: # 遍历读入的基因
        pro = name2seq[gene.pro_id] # sequence
        cov_poses, cov_peps = GetProCov(pro, trie) # 本蛋白质被覆盖的区域以及相应的肽段
        pre_base = 0 # 已考察核苷酸的长度
        new_segs = [] # 存储本基因中被保存下来的segs
        for seg in gene.segs: # 遍历本基因包含的所有segs
            covered = False # 标记此seg是否被cover到 
            end = int((pre_base + seg.length) / 3)
            for i in range(len(cov_poses)):
                if cov_poses[i][0] >= end or cov_poses[i][1] <= pre_base / 3: # 完全未覆盖
                    continue
                else:
                    if cov_poses[i][0] >= pre_base / 3 and cov_poses[i][1] <= end: # 完全覆盖
                        seg.peps.append(cov_peps[i])
                        covered = True
                    elif pep_uni_dic[cov_peps[i]]: # 部分覆盖，若为unique肽段，则认为是覆盖
                        seg.peps.append(cov_peps[i])
                        covered = True
                    else: # 部分覆盖，且为非unique肽段，至少需要覆盖3个氨基酸，且肽段长度>=9
                        if cov_poses[i][0] <= pre_base / 3 and cov_poses[i][1] - pre_base / 3 >= 2 and len(cov_peps[i]) >= 9:
                            seg.peps.append(cov_peps[i])
                            covered = True
                        elif cov_poses[i][0] >= pre_base / 3 and end - cov_poses[i][0] >= 2 and len(cov_peps[i]) >= 9:
                            seg.peps.append(cov_peps[i])
                            covered = True
            if covered:
                new_segs.append(seg)
            pre_base += seg.length # 更新考察核苷酸数
        gene.filtered_segs = new_segs # 有质谱数据覆盖的gff段
        if gene.is_del:
            rec_seg_num += len(gene.filtered_segs)
    return rec_seg_num

"""判断seq是否被数据库所包含"""
def IsSubPro(seq, pros):
    for pro in pros:
        if seq in pro:
            return True
    return False

"""根据目标基因是否与gff+DNA序列相对应，过滤目标基因
具体规则：在翻译库trans_pro中是否存在该子串

note：感觉这个有点问题，翻译库中不存在，也许是由于aa最短长度限制，并不是删掉的基因不合法"""
def GeneFilterbyGFF(name2seq, trans_pros):
    for gene in genes:
        pro = name2seq[gene.pro_id]
        seqs = PartPro2Seg(gene, pro) # 得到每一段的aa序列
        for seq in seqs:
            if not IsSubPro(seq, trans_pros):
                print(f"Not match: {gene.pro_id}")
                gene.legal = False
                break

def ReadAnno(anno_file):
    with open(anno_file,'r') as f:
        f.readline() # 过滤掉第一行
        line = f.readline()
        i = 1 # 该行index为1
        while line:
            segs = line.split('\t')
            if len(segs) < 21:
                continue
            min_pos = min(int(segs[11]), int(segs[12]))
            max_pos = max(int(segs[11]), int(segs[12]))
            if (min_pos > max_pos):
                print(line)
            assert min_pos <= max_pos
            chr_name = re.split("_", segs[4], 2)[2]
            anno_res.append(AnnoRes(chr_name, segs[5], segs[13], min_pos, max_pos, segs[15], i, False, segs[17], line, segs[8].split(";"), [] if segs[16] == "." else segs[16].split(";")))
            line = f.readline()
            i += 1
    print(f"Read .panno file end. And the total report number is {len(anno_res)}")

def ReadAnnoRet(anno_file):
    anno_res = []
    with open(anno_file,'r') as f:
        f.readline() # 过滤掉第一行
        line = f.readline()
        i = 1 # 该行index为1
        while line:
            segs = line.split('\t')
            if len(segs) < 21:
                continu
            min_pos = min(int(segs[11]), int(segs[12])) #if segs[13] == "+" else int(segs[12])
            max_pos = max(int(segs[12]), int(segs[11])) #if segs[13] == "+" else int(segs[11])
            assert min_pos <= max_pos
            chr_name = re.split("_", segs[4], 2)[2]
            anno_res.append(AnnoRes(chr_name, segs[5], segs[13], min_pos, max_pos, segs[15], i, False, segs[17], line, segs[8].split(";"), [] if segs[16] == "." else segs[16].split(";")))
            line = f.readline()
            i += 1
    print(f"Read .panno file end. And the total report number is {len(anno_res)}")
    return anno_res

"""读取蛋白质名"""
def ReadName (fasta_file):
    names = set()
    for line in open(fasta_file):
        if line[0] == '>':
            names.add(line.split()[0][1:])
    return names

def ReadpFindPeps(res_file):
    """读取pFind鉴定到的sequence"""
    peps = set()
    pep_uni_dic = {} # key = sequence, value = true(unique peptide) or false
    with open(res_file) as file:
        lines = file.readlines()
    del lines[0]
    for line in lines:
        segs = line.split('\t')
        pep = segs[5]
        if pep in peps or len(pep) <= 6:
            continue
        peps.add(pep)
        pros = segs[12].split('/')
        cnt_tpro = 0
        for pro in pros:
            if len(pro) >= 3 and pro[0 : 3] != "REV":
                cnt_tpro += 1
        if cnt_tpro == 1:
            pep_uni_dic[pep] = True
        else:
            pep_uni_dic[pep] = False
    return peps, pep_uni_dic

"""扩张含IL氨基酸的肽段"""
def ILExtend(pep):
    if "L" not in pep and "I" not in pep:
        return [pep]
    res = []
    pos = []
    for i in range (len(pep)): # 存储所有出现IL的位置
        if pep[i] == "L" or pep[i] == "I":
            pos.append(i)
    cnt = pow(2, len(pos))
    s = list(pep)
    for i in range(cnt):
        code = bin(i)[2:].zfill(len(pos))
        for j in range(len(code)):
            if code[j] == "1": # "L"
                s[pos[j]] = "L"
            else:
                s[pos[j]] = "I"
        new_pep = ''.join(s)
        res.append(new_pep)
    return res

"""判断肽段是否都能回贴到pro上"""
def IsPepMatch(pro, peps):
    for pep in peps:
        if pep in pro:
            continue
        is_match = False
        pep_ils = ILExtend(pep)
        for pep_il in pep_ils: # 只要有一个扩张肽段能match即可
            if pep_il in pro:
                is_match = True
                break
        if not is_match:
            return False
    return True

"""判断dna_seg是否被re_anno_pos所包含，即位置是否一致"""
def IsContainSeg(re_anno_pos, dna_seg):
    for pos in re_anno_pos:
        if pos[0] == dna_seg.start and pos[1] == dna_seg.end:
            return True
    return False

"""判断某pAnno报告是否合法
anno_res: 待考察报告
genes: 目标基因"""

"""计算外显子precision recall
name2seq: key = pro_id value = aa seqs
genes: 被删除的gene
anno_res: pAnno报告的行
index: 待考察pAnno报告行的index
"""
def TestAnnoRes(name2seq, genes, anno_res, index):
    rec_segs = [] # 存储被召回的gene seg的信息，存储的是pair (i,j)，标记第i个的第j个seg被召回
    re_anno_pos = [] # 本新基因报告所注释到的片段位置
    for re_anno_segs in anno_res[index].anno_info: # 遍历本新基因报告重注释到的gff行
        a = int(re_anno_segs.split(",")[0].split(":")[1].strip('"'))
        b = int(re_anno_segs.split(",")[1].split(":")[1].strip('"'))
        re_anno_pos.append([min(a, b), max(a, b)])
    for gene in genes: # Precision
        if not gene.legal: # 过滤掉目标基因中的脏数据
            continue
        if gene.chr_name == anno_res[index].chr_name and gene.pn == anno_res[index].pn:
            for seg in gene.segs: # Precision
                if IsContainSeg(re_anno_pos, seg):
                    anno_res[index].legal = True
                    # seg.is_recall = True
    i = 0
    for gene in genes: # Recall
        if not gene.legal or not gene.is_del: # 仅考虑被删除的基因
            i += 1
            continue
        j = 0 # 段下标
        if gene.chr_name == anno_res[index].chr_name and gene.pn == anno_res[index].pn:
            for seg in gene.filtered_segs: # Recall
                if IsContainSeg(re_anno_pos, seg):
                    seg.is_recall = True
                    rec_segs.append((i, j)) # 第i个基因的第j个seg被召回
                j += 1
        i += 1
    return anno_res[index], set(rec_segs)

def IsSupTP(anno_line, unp_pros, ncb_pros) -> bool:
    for pro in ncb_pros: # 是来自NCBI库的子串，说明是误匹配的其他位置
        if anno_line.aa_seq in pro:
            return False
    for pro in unp_pros: # 能匹配到其他库，猜测是由于评测实验的库考虑不全所致
        if anno_line.aa_seq in pro:
            anno_line.legal = True
            return True
    return False

"""考察pAnno报告的所有结果是否为TP报告,并计算precision recall
anno_res: 从.panno文件读取的全部报告结果
genes: 被删除的目标基因(希望能被.panno召回)
seg_num: 被删除的gff行数(recall分母)
"""
def pAnnoResPart (seg_num, unp_pros, ncb_pros, name2seq, anno_file, right_file, not_right_file):
    TP_p = 0
    TP_r = 0
    TP_sup = 0
    true_file = open(right_file, "w")
    true_file.write(title_csv)
    false_file = open(not_right_file, "w")
    false_file.write(title_csv)

    TestAnnoRes_1 = partial(TestAnnoRes, name2seq, genes, anno_res)
    pool = multiprocessing.Pool(processes=30)
    results = list(tqdm(pool.imap(TestAnnoRes_1, [i for i in range(len(anno_res))], chunksize=2000), total=len(anno_res), ascii=True))
    rec_segs = set()
    anno_tf = []
    for a, b in results:
        anno_tf.append(a)
        rec_segs.update(b)
    # for anno_line in tqdm(anno_res, ascii=True): # 考察每一个报告结果
    #     TestAnnoRes(anno_line, genes, name2seq)

    for anno_line in tqdm(anno_tf, ascii=True):
        if anno_line.legal:
            TP_p += 1
            true_file.write(anno_line.total_line)
        #else:
        #    if IsSupTP(anno_line, unp_pros, ncb_pros):
         #      TP_sup += 1
          #  false_file.write(anno_line.total_line)
    TP_r = len(rec_segs)
    
    print(f"TP_sup : {TP_sup}")
    print(f"true identified anno line {TP_p}")
    print(f"Recalled segs: {TP_r}")
    print(f"Recall: {TP_r/ seg_num}") #召回率
    print(f"recall denominator  {seg_num}")
    print(f"Pre: {TP_p / len(anno_res)}") #准确率
    print(f"Pre denominator  {len(anno_res)}")
    print(f"Pre add sup pros: {(TP_p + TP_sup) / len(anno_res)}") #准确率
    print(f"Pre denominator  {len(anno_res)}")
    return rec_segs

def AnnoLineIsSub (lines, test_line):
    for tmp_line in lines:
        if test_line == tmp_line:
            return -1
    return test_line.line_num

"""求两个anno结果的差集，写出文件anno_res1-anno_res2"""
def AnnoResSub(anno_res1, anno_res2, res_file, out_file):
    out_file = open(out_file, "w")
    out_file.write(title_csv)
    indexs = []
    AnnoLineIsSub_1 = partial(AnnoLineIsSub, anno_res2)
    pool = multiprocessing.Pool(processes=30)
    results = list(tqdm(pool.imap(AnnoLineIsSub_1, [test_line for test_line in anno_res1], chunksize=5000), total=len(anno_res1), ascii=True))
    # for anno_line in anno_res1:
    #     is_write = True
    #     for tmp_line in anno_res2:
    #         if anno_line == tmp_line:
    #             is_write = False
    #             break
    #     if is_write:
    #         indexs.append(anno_line.line_num)
    with open(res_file,'r') as f: # 将报告结果分两类输出,TP和FP
        panno_lines = f.readlines()
    for index in results:
        if index != -1:
            out_file.write(panno_lines[index])

"""处理spectra文件，得到已解析scan set，所以这里的输出.spectra文件应该是搜索全部已注释结果的"""
def GetIdScan (spec_file):
    scan_set = set()
    with open(spec_file) as file:
        lines = file.readlines()
    del lines[0]
    for line in lines:
        scan_set.add(line.split('\t', 1)[0].rsplit('.', 4)[0])
    return scan_set

"""处理搜索定制数据库的结果，删除灰色地带的谱图匹配，生成新的.spectra文件"""
def DelGrayPSM (spec_file, scan_set, out_file):
    of = open(out_file, "w")
    with open(spec_file) as file:
        lines = file.readlines()
    of.write(lines[0])
    del lines[0]
    for psm in tqdm(lines, ascii=True):
        scan = psm.split('\t', 1)[0].rsplit('.', 4)[0]
        if scan in scan_set:
            of.write(psm)
    of.close()

"""过滤肽段，仅保留不被蛋白质序列包含的那些肽段"""
def PepFilterbyPro (trie, peps, pros):
    matches = set() # 存储被覆盖到的肽段
    for pro in pros:
        for end_index, string in trie.iter(pro):
            matches.add(string)
    res = peps - matches # 返回差集，即没有被模拟已注释库覆盖的肽段
    print(f"matched pep: {len(matches)}")
    print(f"filtered pep: {len(res)}")
    return res

if __name__ == "__main__":    
    start_time = time.time()
    all_data = r"G:\kfwang\human\NCBI-download\protein.faa" # 全体正确的蛋白质序列
    target_data = r"G:\kfwang\human\database\target100%.fasta" # 被覆盖的库
    retain_data = r"G:\kfwang\human\database\target_group_pace2.fasta" # 模拟不完整的已注释库
    trans_db_file = r"G:\kfwang\human\PaceEva\2\trans_res\trans_DNA.fasta" # DNA翻译库
    panno_file = r"G:\kfwang\human\PaceEva\2-ideal\output\dna_anno-filtered-group.panno" # 待评测的pAnno结果
    not_right_file = r"G:\kfwang\human\PaceEva\2-ideal\output\test_anno+score_nr.panno"
    right_file = r"G:\kfwang\human\PaceEva\2-ideal\output\test_anno+score_r.panno"
    gff_file = r"G:\kfwang\human\NCBI-download\GCF_000001405.40_GRCh38.p14_genomic.gff" # gff文件
    sup_pro_file = r"G:\kfwang\human\NCBI-download\protein.faa" # 所有可能正确的蛋白质（解释为虽然结果不在目标库中，但能找到别的证据）
    pfind_res = r"G:\kfwang\human\database\pFind_Full\pFind-Filtered.spectra" # pFind_full .spectra 文件
    not_rec_gff = r"G:\kfwang\human\PaceEva\2-ideal\output\not_rec.panno"
    #second_stage_res = r"G:\kfwang\human\PaceEva\2\output\pFind_res2\pFind-Filtered-raw.spectra"
    #fir_stage_res = r"G:\kfwang\human\PaceEva\2\output\pFind_res1\pFind-Filtered.spectra"
    #fir_peps = ReadpFindPeps(fir_stage_res)
    anno_pros = set(ReadProSeq(retain_data)) # 模拟已注释库的蛋白质序列
    peps, pep_uni_dic = (ReadpFindPeps(pfind_res))# - ReadpFindPeps(fir_stage_res)) & ReadpFindPeps(second_stage_res) # 有质谱证据支持的全部肽段（完整目标库）
    print(f"$$$$ {len(peps)}")
    """得到有质谱证据支持的，且没有被模拟已注释库包含的肽段（这些肽段才是可能正确的novel peptides）"""
    trie = ahocorasick.Automaton() # Create a trie data structure
    for pep in peps:
        trie.add_word(pep, pep)
    trie.make_automaton()    # Build the trie data structure
    filtered_peps = PepFilterbyPro (trie, peps, anno_pros) # 有质谱证据支持的非模拟已注释库肽段
    print(f"$$$$ filtered pep {len(filtered_peps)}") # 49k vs. 41k
    trie_filt = ahocorasick.Automaton() # Create a trie data structure
    for pep in filtered_peps:
        trie_filt.add_word(pep, pep)
    trie_filt.make_automaton()    # Build the trie data structure
    
    
    sup_pros = set(ReadProSeq(sup_pro_file)) # 存在证据的所有蛋白质
    ncb_pros = sup_pros - anno_pros # 存在证据的所有蛋白质序列（删掉模拟已注释库的部分）
    ncb_name2seq = BuildDic(all_data) # key = protein name; value = sequence

    del_pro_name = ReadName (target_data) - ReadName (retain_data) # 需要被召回的目标蛋白名（是right_pro_name的子集）
    right_pro_name = ReadName (all_data) - ReadName (retain_data) # 报告的可以被判定为正确的蛋白质名
    seg_num_raw = DealGFF(gff_file, right_pro_name, del_pro_name) # 将待召回蛋白对应到gff文件，以gene形式存储，同时得到总段数
    
    print(f"seg_num_raw: {seg_num_raw}")
    rec_seg_num = RefineGeneRec(del_pro_name, ncb_name2seq, trie_filt, pep_uni_dic) # 精细化处理gene所对应的gff行，仅保留能被召回的
    print(f"rec_seg_num: {rec_seg_num}")
    # GeneFilterbyGFF(genes, ncb_name2seq, ReadProSeq(trans_db_file))
    ReadAnno(panno_file) # 读取panno报告的所有行
    rec_segs = pAnnoResPart (rec_seg_num, sup_pros, ncb_pros, ncb_name2seq, panno_file, right_file, not_right_file) # 得到被召回的seg信息
    not_rec_file = open(not_rec_gff, "w")
    for i in range(len(genes)):
        for j in range(len(genes[i].filtered_segs)):
            if (i, j) not in rec_segs and genes[i].is_del:
                pep_str = ';'.join([str(item) for item in genes[i].filtered_segs[j].peps])
                not_rec_file.write(genes[i].chr_name + "\t" + genes[i].pro_id + "\t" + genes[i].pn + "\t" + str(genes[i].filtered_segs[j].start) + "\t" + str(genes[i].filtered_segs[j].end) + "\t" + pep_str + "\n")
    end_time = time.time()
    print("脚本运行时间为：", end_time - start_time, "秒")
