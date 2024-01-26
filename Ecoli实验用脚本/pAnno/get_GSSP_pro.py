# -*- coding: utf-8 -*-
"""提取GSSP，规则变为比对蛋白质名称"""
def GetProteinName(Filename):
    #提取蛋白库中的名称
    res = set()
    with open(Filename) as f1:
        f11 = f1.readlines()
    i = 0
    while i < len(f11):
        if f11[i][0] == '>':
            res.add(f11[i].strip().split()[0].split('|')[0][1 : ])
        i += 1
    print("protein num : ", len(res))
    return res

def PfindPro(psm):
    #psm(一行)信息中提取蛋白质名
    res = set()
    proteins = psm.split('\t')[12].split('/')
    for i in range(0, len(proteins) - 1):
        res.add(proteins[i].strip())
    return res

def TestIn(protein_name, proteins):
    for pro in proteins:
        for name in protein_name:
            #print(name)
            #print(pro)
            if name in pro:
                return True
    return False

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
    protein_name = GetProteinName(database_path)
    res = set()
    with open(spectra_path) as f1:
        f11 = f1.readlines()
    for i in range(1, len(f11)):
        if i % 100 == 0:
            print(i)
        proteins = PfindPro(f11[i])
        if not TestIn(protein_name, proteins):
            res.add(f11[i].split('\t')[5])
    print("GSSP num:  ", len(res))
    with open(output_path, "w") as f:
        title = ">pFind_GSSP\n"
        for r in res:
            f.write(title)
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
    output_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\open\\pFind_gssp.txt"
    #output_path = "C:\\Users\\kfwang\\Desktop\\pt1_GAPE_GSSP.txt"
    GetGSSPPfind(spectra_path, database_path, output_path)
    #GetGSSPGape(spectra_path, database_path, output_path)