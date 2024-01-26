# -*- coding: utf-8 -*-
"""提取各GSSP文件的交集，画venn图"""

def ReadSeq(GSSP_path):
    with open(GSSP_path) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(0, len(f11)):
        if f11[i][0] == '>':
            continue
        setA.add(f11[i])
    return setA

def WriteFile(output_path, out_set):
    with open(output_path, "w") as f:
        for s in out_set:
            f.write(s)

if __name__ == "__main__":
    GAPE_path = "Z:\\kfwang\\pAnno\\pt1\\All_GSSP_only.txt"
    oepn_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\open\\pFind_open_gssp.txt"
    res_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\trypsin\\pFind_r_gssp.txt"

    open_set = ReadSeq(oepn_path)
    res_set = ReadSeq(res_path)
    gape_set = ReadSeq(GAPE_path)
    print("open_set : ", len(open_set))
    print("res_set : ", len(res_set))
    print("gape_set : ", len(gape_set))
    gape_open =  gape_set - open_set
    open_gape = open_set - gape_set
    gape_res = gape_set - res_set
    res_gape = res_set - gape_set
    res_open = res_set - open_set
    open_res = open_set - res_set
    print("gape - open : ", len(gape_open))
    print("open - gape : ", len(open_gape))
    print("gape - res : ", len(gape_res))
    print("res - gape : ", len(res_gape))
    print("res - open : ", len(res_open))
    print("open - res: ", len(open_res))
    
    print("gape & open : ", len(open_set & gape_set))
    print("gape & res : ", len(res_set & gape_set))
    print("res & open : ", len(res_set & open_set))
    print("res & open & gape : ", len(res_set & open_set & gape_set))
    gape_open_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\gape_open.txt"
    open_gape_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\open_gape.txt"
    gape_res_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\gape_res.txt"
    res_gape_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\res_gape.txt"
    open_res_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\open_res.txt"
    res_open_path = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\res_open.txt"
    
    WriteFile(gape_open_path, gape_open)
    WriteFile(open_gape_path, open_gape)
    WriteFile(gape_res_path, gape_res)
    WriteFile(res_gape_path, res_gape)
    WriteFile(open_res_path, open_res)
    WriteFile(res_open_path, res_open)

