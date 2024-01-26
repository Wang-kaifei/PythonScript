import os

def Unite(setA, setB):
    #求集合A和B的交集的元素个数
    unite_num = len(setA.intersection(setB))
    return unite_num, float(unite_num) / min(len(setA), len(setB))

def Union(setA, setB):
    #求集合A和B的并集的元素个数，同时返回交集占小集合的比例
    return len(setA.union(setB))

def Sub(setA, setB):
    #求差集，在A中但是不在B中的元素个数
    return len(setA.difference(setB))

def test_sequence(Filename):
    """处理后缀为.spectra的文件，提取第五列（即sequence集合）,放入set中返回"""
    Filename += "\\pFind-Filtered.spectra"
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(1, len(f11)):
        setA.add(f11[i].split('\t', 6)[5])
    return setA

def test_peptides(Filename):
    """处理后缀为.spectra的文件，提取第五列+第十列（即peptide集合）,放入set中返回"""
    Filename += "\\pFind-Filtered.spectra"
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(1, len(f11)):
        setA.add(f11[i].split('\t', 12)[5] + f11[i].split('\t', 12)[10])
    return setA

def protein_group(Filename):
    """处理后缀为.protein的文件，提取符合要求的protein_group名，放入set中返回"""
    Filename += "\\pFind.protein"
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(2, len(f11)):
        if f11[i].split('\t', 3)[0].isdigit():
            if str.lower(f11[i].split('\t', 3)[1][0 : 3]) != "rev":
                setA.add(f11[i].split('\t', 3)[1])
                # print(f11[i].split('\t', 3)[1])
        elif "--" in f11[i].split('\t', 3)[0]:
            break
    # print(len(setA))
    return setA

def test_PSM(filename, merge):
    """函数功能：从.spectra文件中按照nan文件的格式提取psm"""
    if merge:
        filename += "\\result5.spectra"
    else:
        filename += "\\pFind-Filtered.spectra"
    with open(filename) as f:
        f = f.readlines()
    psms = set()
    for i in range(1, len(f)):
        segs = f[i].split('\t')
        psm = segs[0] + segs[5]
        mods = segs[10].split(';')
        for j in range(0, len(mods) -  1):
            psm += mods[j]
        psms.add(psm)
    return psms

def protein(Filename):
    """函数功能：处理后缀为.protein的文件，提取非诱饵protein名，放入set中返回"""
    Filename += "\\pFind.protein"
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(2, len(f11)):
        if f11[i].split('\t', 3)[0].isdigit():   #找到一个protein_group
            if str.lower(f11[i].split('\t', 3)[1][0 : 3]) != "rev":
                setA.add(f11[i].split('\t', 3)[1])
            j = i + 1   #开始考察下一行
            while f11[j].split('\t', 3)[1] == "SubSet":  #考察这一group的subset
                if str.lower(f11[j].split('\t', 3)[2][0 : 3]) != "rev":
                    setA.add(f11[j].split('\t', 3)[2])
                j += 1
            i = j
        elif "--" in f11[i].split('\t', 3)[0]:
            break
    return setA
    
    
    
    
    
def compare(dir1, dir2):
    """对于dir1和dir2文件夹下的spectra文件，求第五列（文件夹之间）sequence的交并差个数
    和第五列+第十列 peptide 的交并差个数
    protein_group 的交并差个数"""
    result1 = dir1.rsplit("result", 1)[1] + "@" + dir2.rsplit("result", 1)[1] + '\t'
    result2 = dir1.rsplit("result", 1)[1] + "@" + dir2.rsplit("result", 1)[1] + '\t'
    result3 = dir1.rsplit("result", 1)[1] + "@" + dir2.rsplit("result", 1)[1] + '\t'
    setA1 = set()
    setB1 = set()
    setA2 = set()
    setB2 = set()
    setA3 = set()
    setB3 = set()
    for dirpath, dirnames, filenames in os.walk(dir1):
        for file in filenames:
            if file == "pFind-Filtered.spectra":
                setA1 = setA1.union(test_sequence(os.path.join(dirpath, file)))
                setA2 = setA2.union(test_peptides(os.path.join(dirpath, file)))
            elif file == "pFind.protein":
                setA3 = setA3.union(protein_group(os.path.join(dirpath, file)))
    # print("The sequence number of ", dir1, " : ", len(setA1))
    # print("The peptides number of ", dir1, " : ", len(setA2))
    for dirpath, dirnames, filenames in os.walk(dir2):
        for file in filenames:
            if file == "pFind-Filtered.spectra":
                setB1 = setB1.union(test_sequence(os.path.join(dirpath, file)))
                setB2 = setB2.union(test_peptides(os.path.join(dirpath, file)))
            elif file == "pFind.protein":
                setB3 = setB3.union(protein_group(os.path.join(dirpath, file)))
    # print("The sequence number of ", dir2, " : ", len(setB1))
    # print("The peptides number of ", dir2, " : ", len(setB2))
    result1 += str(Unite(setA1, setB1)) + "\t" + str(Union(setA1, setB1)) + "\t" + str(Sub(setA1, setB1)) + "\t" + str(Sub(setB1, setA1)) + "\n"
    result2 += str(Unite(setA2, setB2)) + "\t" + str(Union(setA2, setB2)) + "\t" + str(Sub(setA2, setB2)) + "\t" + str(Sub(setB2, setA2)) + "\n"
    result3 += str(Unite(setA3, setB3)) + "\t" + str(Union(setA3, setB3)) + "\t" + str(Sub(setA3, setB3)) + "\t" + str(Sub(setB3, setA3)) + "\n"
    # print(result1, result2)
    return result1, result2, result3

def sum_sequences_peptides(dirname):
    """为周虎数据集处理特有
    求文件夹dirname下sequence和peptide的总数和protein_group总数"""
    setA = set()
    setB = set()
    setC = set()
    for dirpath, dirnames, filenames in os.walk(dirname):
        for file in filenames:
            if file == "pFind-Filtered.spectra":
                setA = setA.union(test_sequence(os.path.join(dirpath, file)))
                setB = setB.union(test_peptides(os.path.join(dirpath, file)))
            elif file == "pFind.protein":
                setC = setC.union(protein_group(os.path.join(dirpath, file)))
    return len(setA), len(setB), len(setC)

def analyze (list_info, dir_list):
    """函数功能：将list_info中各元素集合相互求交并差
    dir_list表示list_info中各元素的来源"""
    res = ""
    for i in range(0, len(list_info) - 1):
        for j in range(i + 1, len(list_info)):
            res += dir_list[i].rsplit("\\", 1)[1] + " VS " + dir_list[j].rsplit("\\", 1)[1] + "\t"
            temp1, temp2 = Unite(list_info[i], list_info[j])
            res += str(temp1) + "\t"
            res += str(Union(list_info[i], list_info[j])) + "\t"
            res += str(Sub(list_info[i], list_info[j])) + "\t"
            res += str(Sub(list_info[j], list_info[i])) + "\t"
            res += str(temp2) + "\n"
    return res
            

if __name__ == "__main__":
    dir_list = [
            "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\complete1",
                  "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_0.01",
                  "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_0.05",
                  "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database1\\VS_C\\OR1_1",
            "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\complete2",
                  "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database2\\VS_C\\OR2_0.01",
                  "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database2\\VS_C\\OR2_0.05",
                  "C:\\test_wkf\\two_step_test\\new\\Ecoli\\database2\\VS_C\\OR2_1"]                
    #flag = [False, False, False, False, False, False, False, False]
    flag = [False, True, True, True, False, True, True, True]
    res = "A VS B name\tUnite\tUnion\tSub A-B\tSub B-A\tscale\n"
    list_psm = []
    list_seq = []  #按照dir_list中的顺序存储各PSM、序列、肽段、蛋白质、蛋白质群集合
    list_pep = []
    list_pro = []
    list_pro_gro = []
    #获取序列、肽段、蛋白质、蛋白质集合
    for i in range(0, len(dir_list)):
        list_psm.append(test_PSM(dir_list[i], flag[i]))
        print("psm get end!\n")
        list_seq.append(test_sequence(dir_list[i]))
        print("sequences get end!\n")
        list_pep.append(test_peptides(dir_list[i]))
        print("peptides get end!\n")
        list_pro.append(protein(dir_list[i]))
        print("proteins get end!\n")
        list_pro_gro.append(protein_group(dir_list[i]))
        print("protein_groups get end!\n")
    res += "psm\n"
    res += analyze(list_psm, dir_list)
    res += "sequences\n"
    res += analyze(list_seq, dir_list)
    res += "peptides\n"
    res += analyze(list_pep, dir_list)
    res += "proteins\n"
    res += analyze(list_pro, dir_list)
    res += "protein_groups\n"
    res += analyze(list_pro_gro, dir_list)
    with open("C:\\Users\\kfwang\\Desktop\\analysis.txt", "w") as f:
        f.write(str(res))
        f.close()

"""if __name__ == "__main__":
    #查看一次搜索和二次搜索的交并差集
    dir_lists = [["C:\\test_wkf\\two_step_test\\Ecoli\\database1\\O1",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\OO1_0.01",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\OO1_0.05",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\OO1_1",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\RO1_0.01",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\RO1_0.05",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\RO1_1"
                 ],
                ["C:\\test_wkf\\two_step_test\\Ecoli\\database1\\R1",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\OR1_0.01",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\OR1_0.05",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\OR1_1",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\RR1_0.01",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\RR1_0.05",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\RR1_1"
                 ],
                 ["C:\\test_wkf\\two_step_test\\Ecoli\\database2\\O2",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\OO2_0.01",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\OO2_0.05",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\OO2_1",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\RO2_0.01",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\RO2_0.05",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\RO2_1"
                 ],
                ["C:\\test_wkf\\two_step_test\\Ecoli\\database2\\R2",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\OR2_0.01",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\OR2_0.05",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\OR2_1",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\RR2_0.01",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\RR2_0.05",
                  "C:\\test_wkf\\two_step_test\\Ecoli\\database2\\RR2_1"
                 ]
                ]
    res = "A VS B name\tUnite\tUnion\tSub A-B\tSub B-A\tscale\n"
    for dir_list in dir_lists:
        list_seq = []  #按照dir_list中的顺序存储各序列、肽段、蛋白质、蛋白质群集合
        list_pep = []
        list_pro = []
        list_pro_gro = []
        #获取序列、肽段、蛋白质、蛋白质集合
        for path in dir_list:
            list_seq.append(test_sequence(path))
            print("sequences get end!\n")
            list_pep.append(test_peptides(path))
            print("peptides get end!\n")
            list_pro.append(protein(path))
            print("proteins get end!\n")
            list_pro_gro.append(protein_group(path))
            print("protein_groups get end!\n")
        res += "sequences\n"
        res += analyze(list_seq, dir_list)
        res += "peptides\n"
        res += analyze(list_pep, dir_list)
        res += "proteins\n"
        res += analyze(list_pro, dir_list)
        res += "protein_groups\n"
        res += analyze(list_pro_gro, dir_list)
        with open("C:\\Users\\kfwang\\Desktop\\analysis.txt", "a") as f:
            f.write(str(res))
            f.close()
        res = ""

"""



    



