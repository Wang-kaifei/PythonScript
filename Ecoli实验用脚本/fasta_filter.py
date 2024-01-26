# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""脚本功能:去除重复蛋白质，建立新fasta"""
def get_seq(filename):
    res = []
    seq = ""
    with open(filename) as f:
        f = f.readlines()
    for i in range(0, len(f)):
        if f[i][0] == '>':
            res.append(seq)
            seq = ""
            continue
        seq += f[i].strip()
    print (len(res))
    return res

def pro_filter(filename, seq):
    with open(filename) as f:
        f = f.readlines()
    res = ""
    seq_now = ""  #存储当前考察的蛋白质序列
    pro = ""      #存储当前考察的蛋白质信息（包含名称和序列）
    same = 0
    for i in range(0, len(f)):
        if f[i][0] == '>':
            if seq_now not in seq:
                res += pro
            else:
                same += 1
            pro = f[i]
            seq_now = ""
            continue
        pro += f[i]  
        seq_now += f[i].strip()
    print (same)
    print (len(res))
    return res

if __name__ == "__main__":
    fasta1 = "Y:\\chihao_dong111\\dong_111_raw\\ecoli_MG1655_20151014_NC_000913_con.fasta"
    fasta2 = "D:\\users\kfwang\\two-step\\Ecoli\\Ecoli-unreviewed\\database\\uniprot_unreviewed.fasta"
    out = "D:\\users\kfwang\\two-step\\Ecoli\\Ecoli-unreviewed\\database\\new.fasta"
    seq = get_seq(fasta1)
    res = pro_filter(fasta2, seq)
    with open(out, 'w') as f:
        f.write(res)
        f.close()


