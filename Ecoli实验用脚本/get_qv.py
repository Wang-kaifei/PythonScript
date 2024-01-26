# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""





def get_qv(nanfile, spectra):
    qv = []
    with open(nanfile) as f1:
        f1 = f1.readlines()
    with open(spectra) as f2:
        f2 = f2.readlines()
    for i in range(0, len(f1)):
        for j in range(1, len(f2)):
            if f1[i].split('\t')[0] == f2[j].split('\t')[0]:
                qv.append(float(f2[j].split('\t')[4]))
                break
    qv.sort()
    print len(qv)
    return qv

        
if __name__ == "__main__":
    nanfile = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_0.05\\self\\pQuant_spectra.list_rnan"
    spectra = "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_0.05\\pFind-Filtered.spectra"
    qv = get_qv(nanfile, spectra)
    out = "C:\\Users\\kfwang\\Desktop\\qv2.txt"
    with open(out, 'w') as f:
        for q in qv:
            f.write(str(q) + '\n')
        f.close()
                  
