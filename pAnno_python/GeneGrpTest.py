# -*- coding: utf-8 -*-

def BuildDict(res_file):
    res = {}
    with open(res_file) as file:
        lines = file.readlines()
    del lines[0]
    for line in lines: 
        segs = line.split('\t', 6)
        if len(segs) < 5:
            continue
        res[segs[5]] = segs[3]
    print(f"Total number of gene: {len(res)}")
    return res

if __name__ == "__main__":
    res_file1 = r"G:\kfwang\human\PaceEva\2-ideal\output\rna_anno.panno"
    res_file2 = r"G:\kfwang\human\PaceEva\2-ideal\output\rna_anno1.panno"
    dic1 = BuildDict(res_file1)
    dic2 = BuildDict(res_file2)
    for key, value in dic2.items():
        if dic1[key] != value:
            print(key, value)