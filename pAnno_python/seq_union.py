# -*- coding: utf-8 -*-

def Unite(setA, setB):
    unite_num = len(setA.intersection(setB))
    return unite_num, float(unite_num) / min(len(setA), len(setB))

def GapeSeq(Filename):
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(0, len(f11)):
        if f11[i][0] == '>':
            continue
        setA.add(f11[i])
    return setA

def PfindSeq(Filename):
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(1, len(f11)):
        setA.add(f11[i].split('\t', 6)[5])
    return setA



if __name__ == "__main__":
    pFind_file = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\pFind-Filtered.spectra"
    gape_file = "C:\\Users\\kfwang\\Desktop\\pt1_All_seq.txt"
    gape_seq = GapeSeq(gape_file)
    pFind_seq = PfindSeq(pFind_file)
    print("gape_seq ::  ", len(gape_seq))
    print("pFind_seq ::  ", len(pFind_seq))
    print(len(Unite(gape_seq, pFind_seq)))