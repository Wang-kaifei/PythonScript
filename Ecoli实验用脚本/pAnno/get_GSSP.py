# -*- coding: utf-8 -*-
"""从鉴定结果中提取GSSP，先提取sequece，然后再与之前已经存好的蛋白质库seq比对"""
def ProteinDatabase(Filename):
    #提取蛋白库中的seq
    res = set()
    with open(Filename) as f1:
        f11 = f1.readlines()
    i = 0
    while i < len(f11):
        if f11[i][0] == '>':
            protein_ = ""
            i += 1
            while i < len(f11) and f11[i][0] != ">":
                protein_ += f11[i].strip()
                i += 1
            res.add(protein_)
    return res

def PfindSeq(Filename):
    #从pFind鉴定结果提取seq
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(1, len(f11)):
        setA.add(f11[i].split('\t', 6)[5])
    print("Total num  ", len(setA))
    return setA

def GapeSeq(Filename):
    #从Gape鉴定结果中提取seq
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(0, len(f11)):
        if f11[i][0] == '>':
            continue
        setA.add(f11[i])
    print("Total num  ", len(setA))
    return setA

def GetGSSPPfind(spectra_path, database_path, output_path):
    database = ProteinDatabase(database_path)
    seqs = PfindSeq(spectra_path)
    res = set()
    for seq in seqs:
        for pro in database:
            if seq in pro:
                res.add(seq)
                break
    sub =  seqs - res
    print("GSSP num:  ", len(sub))
    with open(output_path, "w") as f:
        #title = ">pFind_GSSP\n"
        for r in sub:
            #f.write(title)
            f.write(r + '\n')

def GetGSSPGape(spectra_path, database_path, output_path):
    database = ProteinDatabase(database_path)
    seqs = GapeSeq(spectra_path)
    res = set()
    for seq in seqs:
        for pro in database:
            if seq in pro:
                res.add(seq)
                break
    sub = seqs - res
    print("GSSP num:  ", len(sub))
    with open(output_path, "w") as f:
        title = ">GAPE_GSSP\n"
        for r in sub:
            f.write(title)
            f.write(r + '\n')

if __name__ == "__main__":
    spectra_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\trypsin\\pFind-Filtered.spectra"
    #spectra_path = "C:\\Users\\kfwang\\Desktop\\pt1_GAPE_seq.txt"
    database_path = "C:\\Users\\kfwang\\Desktop\\Phatr2_mitochondrion_chloroplast_Proteomics.fasta"
    output_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\trypsin\\pFind_r_gssp.txt"
    #output_path = "C:\\Users\\kfwang\\Desktop\\pt1_GAPE_GSSP.txt"
    GetGSSPPfind(spectra_path, database_path, output_path)
    #GetGSSPGape(spectra_path, database_path, output_path)