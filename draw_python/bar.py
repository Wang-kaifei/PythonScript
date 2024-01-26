import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import csv
plt.rcParams["font.sans-serif"]=["Microsoft YaHei"]
plt.rcParams['xtick.labelsize'] = 13 # x轴ticks的size
plt.rcParams['ytick.labelsize'] = 13  # y轴ticks的size

def GetBar2c(x1, y1, x2, y2, filepath):
    """条形图，两个子图"""
    fig, axes = plt.subplots(1, 2)
    sns.set_style("whitegrid")
    # 构建数据
    g = sns.barplot(x = x1, y = y1, palette="Blues_d", linewidth=1, ax = axes[0])
    plt.setp(g.patches, linewidth=0)
    for i in range(0, len(x1)):
        g.text(i, y1[i] + 2, round(y1[i] + 2), ha = "center")
    g = sns.barplot(x = x2, y = y2, palette="Blues_d", linewidth=1, ax = axes[1])
    plt.setp(g.patches, linewidth=0)
    for i in range(0, len(x2)):
        g.text(i, y2[i] + 2, round(y2[i] + 2), ha = "center")
    axes[0].set_title("Open")
    axes[1].set_title("Restricted")
    axes[0].set_ylabel("Number")
    axes[0].set_xlabel("Classification")
    axes[1].set_xlabel("Classification")
    plt.savefig(filepath, dpi=200)
    plt.show()

def GetBar1(x, y, out_path, is_rot = False):
    """条形图"""
    sns.set_style("whitegrid")
    if is_rot:
        plt.xticks(rotation=20)   # 设置横坐标显示的角度，角度是逆时针，自己看
    #x = [i/10 for i in x]
    g = sns.barplot(x = x, y = y, palette="Blues_d", linewidth=1)
    plt.setp(g.patches, linewidth=0)
    for i in range(0, len(x)):
        g.text(i, y[i], round(y[i]), ha = "center", fontsize=13)
    plt.savefig(out_path, dpi=200)
    plt.show()


def GetBar2(dic : dict, out_path):
    """条形图"""
    sns.set_style("whitegrid")
    x = list(dic.keys())
    y = list(dic.values())
    print(x)
    g = sns.barplot(x = x, y = y, palette="Blues_d", linewidth=1)
    plt.setp(g.patches, linewidth=0)
    for i in range(0, len(x)):
        g.text(i, y[i], round(y[i]), ha = "center", fontsize=13)
    plt.savefig(out_path, dpi=200)
    plt.show()


def GetAnnoQvalue(anno_file):
    """获取pAnno结果文件中q_value的分布情况，分成10个区间"""
    res = [0] * 10
    with open(anno_file) as file:
        lines = file.readlines()
    for i in range(1, len(lines)):
        segs = lines[i].split('\t')
        q_value = float(segs[6]) * 1000
        res[int(q_value)] = res[int(q_value)] + 1
    return res


def GetAnnoStart(anno_file):
    """获取pAnno结果文件中起始密码子的分布情况"""
    start_code_index = {'ATG': 0, 'GTG': 1, 'CTG': 2, 'TTG': 3, 'AAT': 4, 'ATA': 5, 'ZZZ': 6}
    res = [0] * 7
    with open(anno_file) as file:
        lines = file.readlines()
    for i in range(1, len(lines)):
        segs = lines[i].split('\t')
        code = segs[9].strip()
        if start_code_index[code] == 6:
            print(code)
        res[start_code_index[code]] = res[start_code_index[code]] + 1
    return res


def GetAnnoLength(anno_file):
    """获取pAnno结果文件中蛋白质序列长度分布
    <20 20-40 40-60 60-80 80-100 100-120 >120"""
    res = {'<20':0, '40':0, '60':0, '80':0, '100':0, '120':0, '140':0, '>140':0}
    with open(anno_file) as file:
        lines = file.readlines()
    for i in range(1, len(lines)):
        segs = lines[i].split('\t')
        length = len(segs[16].strip())
        if length < 20: 
            res['<20'] += 1
        elif length < 40:
            res['40'] += 1
        elif length < 60:
            res['60'] += 1
        elif length < 80:
            res['80'] += 1
        elif length < 100:
            res['100'] += 1
        elif length < 120:
            res['120'] += 1
        elif length < 140:
            res['140'] += 1
        else:
            res['>140'] += 1
    return res


def ReadFromFeature(feature_path):
    """从feature.csv文件中提取length、起始子、q_value
    分别对应着0，1，2列"""
    q_value = [0] * 10
    length = {'<20':0, '40':0, '60':0, '80':0, '100':0, '120':0, '140':0, '>140':0}
    code = [0] * 7
    with open(feature_path, "r") as f:
        reader = csv.reader(f)
        lines = list(reader)
        for line in lines:
            l = float(line[0])
            if l < 20: 
                length['<20'] += 1
            elif l < 40:
                length['40'] += 1
            elif l < 60:
                length['60'] += 1
            elif l < 80:
                length['80'] += 1
            elif l < 100:
                length['100'] += 1
            elif l < 120:
                length['120'] += 1
            elif l < 140:
                length['140'] += 1
            else:
                length['>140'] += 1
            q = float(line[2]) * 1000
            q_value[int(q)] = q_value[int(q)] + 1
            code[int(float(line[1]))] = code[int(float(line[1]))] + 1
    return q_value, length, code


if __name__ == "__main__":
    right_anno = "C:\\Users\\pFind\\Desktop\\dl_model_data\\1\\feature_1.csv"
    n_right_anno = "C:\\Users\\pFind\\Desktop\\dl_model_data\\1\\feature_0.csv"
    start_codes = ['ATG', 'GTG', 'CTG', 'TTG', 'AAT', 'ATA', 'ZZZ']
    x = [i / 1000.0 for i in range(1, 11)]
    # r_q_value = GetAnnoQvalue(right_anno)
    # nr_q_value = GetAnnoQvalue(n_right_anno)
    # r_start_code = GetAnnoStart(right_anno)
    # nr_start_code = GetAnnoStart(n_right_anno)
    # r_length = GetAnnoLength(right_anno)
    # nr_length = GetAnnoLength(n_right_anno)
    r_file_qvalue = "D:\\OneDrive - mails.ucas.ac.cn\\组内工作汇报\\2022.04.08\\ecoli\\r_qvalue.png"
    n_r_file_qvalue = "D:\\OneDrive - mails.ucas.ac.cn\\组内工作汇报\\2022.04.08\\ecoli\\nr_qvalue.png"
    r_file_code = "D:\\OneDrive - mails.ucas.ac.cn\\组内工作汇报\\2022.04.08\\ecoli\\r_code.png"
    n_r_file_code = "D:\\OneDrive - mails.ucas.ac.cn\\组内工作汇报\\2022.04.08\\ecoli\\nr_code.png"
    r_file_length = "D:\\OneDrive - mails.ucas.ac.cn\\组内工作汇报\\2022.04.08\\ecoli\\r_length.png"
    n_r_file_length = "D:\\OneDrive - mails.ucas.ac.cn\\组内工作汇报\\2022.04.08\\ecoli\\nr_length.png"
    # GetBar1(x, r_q_value, r_file_qvalue)
    # GetBar1(x, nr_q_value, n_r_file_qvalue)
    # GetBar1(start_codes, r_start_code, r_file_code)
    # GetBar1(start_codes, nr_start_code, n_r_file_code)
    # GetBar2(r_length, r_file_length)
    # GetBar2(nr_length, n_r_file_length)
    
    # feature_path = "C:\\Users\\pFind\\Desktop\\dl_model_data\\yeast\\neg.csv"
    fe_q_value, fe_length, fe_codes = ReadFromFeature(right_anno)
    GetBar1(x, fe_q_value, r_file_qvalue, True)
    GetBar1(start_codes, fe_codes, r_file_code)
    GetBar2(fe_length, r_file_length)
    fe_q_value, fe_length, fe_codes = ReadFromFeature(n_right_anno)
    GetBar1(x, fe_q_value, n_r_file_qvalue, True)
    GetBar1(start_codes, fe_codes, n_r_file_code)
    GetBar2(fe_length, n_r_file_length)
    
    
    