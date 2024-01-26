'''
Descripttion: 
version: 
Author: Kaifei
Date: 2024-01-23 21:43:35
LastEditors: Kaifei
LastEditTime: 2024-01-24 23:38:26
'''
# -*- coding: utf-8 -*-

from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import seaborn as sns

def CometPep(res_path):
    seqs = set()
    with open(res_path, 'r') as f:
        ff = f.readlines()
    # 0:scan	1:num	2:charge	3:exp_neutral_mass	4:calc_neutral_mass	5:e-value	6:xcorr	7:delta_cn	8:sp_score	9:ions_matched	10:ions_total	11:plain_peptide	12:modified_peptide	13:prev_aa	14:next_aa	15:protein	16:protein_count	17:modifications	18:retention_time_sec	19:sp_rank	20:Target	21:File	22:q_value
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        seqs.add(segs[11])
    print(f"Comet: {len(seqs)}")
    return seqs, "Comet"

def msgf_mod_seq(raw_pep):
    mod_dic = {"+15.995":"Oxidation[M]", "+57.021":"Carbamidomethyl[C]", "-17.027":"Gln->pyro-Glu[AnyN-termQ]", "+42.011":"Acetyl[ProteinN-term]"}
    # 如果最左端出现（N端）修饰，则index是0
    mods = ""
    seq = "" # pFind格式的序列
    aa_index = 0 # 目前有多少个氨基酸出现
    temp_mod = "" # 暂时存储的修饰
    mid = raw_pep.split('.', 1)[1].rsplit('.', 1)[0]
    for cc in mid:
        if cc.isalpha():
            if temp_mod != "": # 如果前面的氨基酸有修饰，需要存储
                mods += str(aa_index) + "," + mod_dic[temp_mod] + ";"
                temp_mod = ""
            seq += cc
            aa_index += 1
        elif cc == "+" or cc == "-":
            if temp_mod != "": # 如果前面的氨基酸有修饰，需要存储
                mods += str(aa_index) + "," + mod_dic[temp_mod] + ";"
            temp_mod = cc
        else:
            temp_mod += cc
    return seq, mods
  
def MSGFPep(res_path):
    seqs = set()
    with open(res_path, 'r') as f:
        ff = f.readlines()
    # #Spec_file:0	Spec_ID:1	ScanNum:2	Scan_time:3	FragMethod:4	Precursor:5	IsotopeError:6	Precursor_Error:7	charge:8	Peptide:9	Protein:10	DeNovoScore:11	MSGFScore:12	SpecEvalue:13	Evalue:14	target:15	q_value:16
    for i in range(1, len(ff)):    
        segs = ff[i].strip().split('\t')
        seqence, _ = msgf_mod_seq(segs[9])
        seqs.add(seqence)
    print(f"MSGF+: {len(seqs)}")
    return seqs, "MSGF+"

def MaxquantPep(res_path):
    seqs = set()
    with open(res_path, 'r') as f:
        ff = f.readlines()
    # Raw file:0	Scan number:1	Scan index:2	Sequence:3	Length:4	Missed cleavages:5	Modifications:6	Modified sequence:7	Oxidation (M) Probabilities:8	Oxidation (M) Score diffs:9	Acetyl (Protein N-term):10	Gln->pyro-Glu:11	Oxidation (M):12	Proteins:13	Charge:14	Fragmentation:15	Mass analyzer:16	Type:17	Scan event number:18	Isotope index:19	m/z:20	Mass:21	Mass error [ppm]:22	Mass error [Da]:23	Simple mass error [ppm]:24	Retention time:25	PEP:26	Score:27	Delta score:28	Score diff:29	Localization prob:30	Combinatorics:31	PIF:32	Fraction of total spectrum:33	Base peak fraction:34	Precursor full scan number:35	Precursor Intensity:36	Precursor apex fraction:37	Precursor apex offset:38	Precursor apex offset time:39	Matches:40	Intensities:41	Mass deviations [Da]:42	Mass deviations [ppm]:43	Masses:44	Number of matches:45	Intensity coverage:46	Peak coverage:47	Neutral loss level:48	ETD identification type:49	Reverse:50	All scores:51	All sequences:52	All modified sequences:53	Reporter PIF:54	Reporter fraction:55	id:56	Protein group IDs:57	Peptide ID:58	Mod. peptide ID:59	Evidence ID:60	Oxidation (M) site IDs:61
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        seqs.add(segs[3])
    print(f"Maxquant: {len(seqs)}")
    return seqs, "Maxquant"

def pFindPep(res_path):
    seqs = set()
    with open(res_path, 'r') as f:
        ff = f.readlines()
    # Raw file:0	Scan number:1	Scan index:2	Sequence:3	Length:4	Missed cleavages:5	Modifications:6	Modified sequence:7	Oxidation (M) Probabilities:8	Oxidation (M) Score diffs:9	Acetyl (Protein N-term):10	Gln->pyro-Glu:11	Oxidation (M):12	Proteins:13	Charge:14	Fragmentation:15	Mass analyzer:16	Type:17	Scan event number:18	Isotope index:19	m/z:20	Mass:21	Mass error [ppm]:22	Mass error [Da]:23	Simple mass error [ppm]:24	Retention time:25	PEP:26	Score:27	Delta score:28	Score diff:29	Localization prob:30	Combinatorics:31	PIF:32	Fraction of total spectrum:33	Base peak fraction:34	Precursor full scan number:35	Precursor Intensity:36	Precursor apex fraction:37	Precursor apex offset:38	Precursor apex offset time:39	Matches:40	Intensities:41	Mass deviations [Da]:42	Mass deviations [ppm]:43	Masses:44	Number of matches:45	Intensity coverage:46	Peak coverage:47	Neutral loss level:48	ETD identification type:49	Reverse:50	All scores:51	All sequences:52	All modified sequences:53	Reporter PIF:54	Reporter fraction:55	id:56	Protein group IDs:57	Peptide ID:58	Mod. peptide ID:59	Evidence ID:60	Oxidation (M) site IDs:61
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        seqs.add(segs[5])
    print(f"pFind: {len(seqs)}")
    return seqs, "pFind"

def MSFraggerPep(res_path):
    seqs = set()
    with open(res_path, 'r') as f:
        ff = f.readlines()
    # Raw file:0	Scan number:1	Scan index:2	Sequence:3	Length:4	Missed cleavages:5	Modifications:6	Modified sequence:7	Oxidation (M) Probabilities:8	Oxidation (M) Score diffs:9	Acetyl (Protein N-term):10	Gln->pyro-Glu:11	Oxidation (M):12	Proteins:13	Charge:14	Fragmentation:15	Mass analyzer:16	Type:17	Scan event number:18	Isotope index:19	m/z:20	Mass:21	Mass error [ppm]:22	Mass error [Da]:23	Simple mass error [ppm]:24	Retention time:25	PEP:26	Score:27	Delta score:28	Score diff:29	Localization prob:30	Combinatorics:31	PIF:32	Fraction of total spectrum:33	Base peak fraction:34	Precursor full scan number:35	Precursor Intensity:36	Precursor apex fraction:37	Precursor apex offset:38	Precursor apex offset time:39	Matches:40	Intensities:41	Mass deviations [Da]:42	Mass deviations [ppm]:43	Masses:44	Number of matches:45	Intensity coverage:46	Peak coverage:47	Neutral loss level:48	ETD identification type:49	Reverse:50	All scores:51	All sequences:52	All modified sequences:53	Reporter PIF:54	Reporter fraction:55	id:56	Protein group IDs:57	Peptide ID:58	Mod. peptide ID:59	Evidence ID:60	Oxidation (M) site IDs:61
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        seqs.add(segs[0])
    print(f"MSFragger: {len(seqs)}")
    return seqs, "MSFragger"

def Get1idPep(pep1, pep2, pep3, pep4, pep5):
    """得到至少1个引擎鉴定的肽段"""
    union = pep1 | pep2 | pep3 | pep4 | pep5
    print(f"1id pep: {len(union)}")
    return union

def GetCredPep(pep1, pep2, pep3, pep4, pep5):
    """从5个肽段集合中得到可信肽段"""
    inter12 = pep1 & pep2
    inter13 = pep1 & pep3
    inter14 = pep1 & pep4
    inter15 = pep1 & pep5
    inter23 = pep2 & pep3
    inter24 = pep2 & pep4
    inter25 = pep2 & pep5
    inter34 = pep3 & pep4
    inter35 = pep3 & pep5
    inter45 = pep4 & pep5
    #求并集
    union = inter12 | inter13 | inter14 | inter15 | inter23 | inter24 | inter25 | inter34 | inter35 | inter45
    print(f"Credible pep: {len(union)}")
    # 将pep写出
    with open("credible_pep.txt", 'w') as f:
        for i in union:
            f.write(i + '\n')
    return union

def Get3idPep(pep1, pep2, pep3, pep4, pep5):
    """得到至少3个引擎鉴定的肽段"""
    inter123 = pep1 & pep2 & pep3
    inter124 = pep1 & pep2 & pep4
    inter125 = pep1 & pep2 & pep5
    inter134 = pep1 & pep3 & pep4
    inter135 = pep1 & pep3 & pep5
    inter145 = pep1 & pep4 & pep5
    inter234 = pep2 & pep3 & pep4
    inter235 = pep2 & pep3 & pep5
    inter245 = pep2 & pep4 & pep5
    inter345 = pep3 & pep4 & pep5
    #求并集
    union = inter123 | inter124 | inter125 | inter134 | inter135 | inter145 | inter234 | inter235 | inter245 | inter345
    print(f"3id pep: {len(union)}")
    return union

def Get4idPep(pep1, pep2, pep3, pep4, pep5):
    """得到至少4个引擎鉴定的肽段"""
    inter1234 = pep1 & pep2 & pep3 & pep4
    inter1235 = pep1 & pep2 & pep3 & pep5
    inter1245 = pep1 & pep2 & pep4 & pep5
    inter1345 = pep1 & pep3 & pep4 & pep5
    inter2345 = pep2 & pep3 & pep4 & pep5
    #求并集
    union = inter1234 | inter1235 | inter1245 | inter1345 | inter2345
    print(f"4id pep: {len(union)}")
    return union
    
def Get5idPep(pep1, pep2, pep3, pep4, pep5):
    """得到5个引擎鉴定的肽段"""
    inter12345 = pep1 & pep2 & pep3 & pep4 & pep5
    print(f"5id pep: {len(inter12345)}")
    return inter12345
    
def Draw(set1, set2, set3, set4, set5):
    # 创建子图对象和绘图的布局
    fig, ax = plt.subplots(nrows=2, ncols=5)
    fig.set_size_inches(18, 8)
    # 定义集合和标签的对应关系
    set_labels = {
        (0, 0): (set1[1], set2[1]),
        (0, 1): (set1[1], set3[1]),
        (0, 2): (set1[1], set4[1]),
        (0, 3): (set1[1], set5[1]),
        (0, 4): (set2[1], set3[1]),
        (1, 0): (set2[1], set4[1]),
        (1, 1): (set2[1], set5[1]),
        (1, 2): (set3[1], set4[1]),
        (1, 3): (set3[1], set5[1]),
        (1, 4): (set4[1], set5[1])
    }
    tmp_sets = [set1, set2, set3, set4, set5]
    row = 0
    col = 0
    for i in range(5):
        for j in range(i+1, 5):
            current_labels = set_labels[(row, col)]
            venn_diagram = venn2([tmp_sets[i][0], tmp_sets[j][0]], set_labels=current_labels, ax=ax[row][col])
            ax[row][col].set_title('{} vs. {}'.format(*current_labels), fontsize=16)
            for label in venn_diagram.set_labels:
                label.set_fontsize(14)  # 设置集合名字体大小
            for label in ax[row][col].texts:  # 遍历每个子图中的文本
                label.set_fontsize(14)  # 设置文本字体大小
            col += 1
            if col == 5:
                row += 1
                col = 0
    # 调整子图之间的间距
    plt.subplots_adjust(hspace=0.4, wspace=0.25)
    # 保存高分辨率图片
    plt.savefig('venn_diagrams.png', format='png', dpi=300, bbox_inches='tight')  # 保存为高分辨率的图片
    # 显示图形
    plt.show()

def DrawCol():
    """输入是一个4 * 5的矩阵"""
    # Ecoli, yeast, human, 拟南芥
    data = [[6892, 1428, 1468, 2818, 11680], [19982, 4153, 5300, 11306, 37173], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
    name = ['Dong-E.coli', 'Xu-Yeast', 'Gygi-Human', 'Arabidopsis thaliana']
    # 设置子图和柱子的数量
    num_subplots = 4
    num_bars = 5
    # 创建子图和 ax 对象
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
    # 设置使用 seaborn 的色卡
    colors = sns.color_palette('Set2')
    # 遍历子图，画柱状图
    for i in range(num_subplots):
        row = i // 2
        col = i % 2
        ax = axes[row][col]
        # 设置 x 轴标签和标题
        # ax.set_xlabel('X Label')
        ax.set_ylabel('#Peptides')
        ax.set_title(name[i])
        # 画柱状图
        sns.barplot(x=list(range(1, num_bars + 1)), y=data[i], ax=ax, palette=colors)
    # 调整子图之间的间距
    plt.tight_layout()
    # 保存高分辨率图片
    plt.savefig('venn_diagrams.png', format='png', dpi=300, bbox_inches='tight')  # 保存为高分辨率的图片
    # 显示图形
    plt.show()    
    
if __name__ == "__main__":
    comet_path = r"E:\wkf\yeast\credible_pep\Comet\res\result.txtfilter"
    msgf_path = r"E:\wkf\yeast\credible_pep\MSGF+\tsv\all_result.tsvfilter"
    maxquant_path = r"E:\wkf\yeast\credible_pep\MaxQuant\combined\txt\msms.txt"
    pFind_path = r"E:\wkf\yeast\credible_pep\pFind\pFind-Filtered.spectra"
    msfragger_path = r"E:\wkf\yeast\credible_pep\MSFragger\peptide.tsv"
    comet_pep, comet_name = CometPep(comet_path)
    msgf_pep, msgf_name = MSGFPep(msgf_path)
    maxquant_pep, maxquant_name = MaxquantPep(maxquant_path)
    pfind_pep, pfind_name = pFindPep(pFind_path)
    msfragger_pep, msfragger_name = MSFraggerPep(msfragger_path)
    # 画韦恩图
    # Draw([comet_pep, comet_name], [msgf_pep, msgf_name], [maxquant_pep, maxquant_name], [pfind_pep, pfind_name], [msfragger_pep, msfragger_name])
    # Get1idPep(comet_pep, msgf_pep, maxquant_pep, pfind_pep, msfragger_pep)
    # 得到可信肽段
    # GetCredPep(comet_pep, msgf_pep, maxquant_pep, pfind_pep, msfragger_pep)
    # Get3idPep(comet_pep, msgf_pep, maxquant_pep, pfind_pep, msfragger_pep)
    # Get4idPep(comet_pep, msgf_pep, maxquant_pep, pfind_pep, msfragger_pep)
    # Get5idPep(comet_pep, msgf_pep, maxquant_pep, pfind_pep, msfragger_pep)
    DrawCol()