# -*- coding: utf-8 -*-
"""将GSSP分成多个文件以供定位问题"""

def PartGSSP(gssp_file, out_path, cnt):
    with open(gssp_file) as f1:
        lines = f1.readlines()
    res = []
    id = 0
    i = 0
    for line in lines:
        if i == cnt:
            out_file = out_path + str(id) + ".txt"
            with open(out_file,'w') as f:
                f.write("".join(res))
            id += 1
            i = 0
            res = []
        res.append(line)
        i += 1
    out_file = out_path + str(id) + ".txt"
    id += 1
    with open(out_file,'w') as f:
        f.write("\n".join(res))


if __name__ == "__main__":
    gssp_file = r"D:\users\kfwang\pAnno\mini\output\first_GSSP.txt"
    out_path = "D:\\users\\kfwang\\pAnno\\mini\\output\\part_gssp\\"
    cnt = 10000
    PartGSSP(gssp_file, out_path, cnt)