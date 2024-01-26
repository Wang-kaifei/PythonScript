import os

def Unite(setA, setB):
    #�󼯺�A��B�Ľ�����Ԫ�ظ���
    unite_num = len(setA.intersection(setB))
    return unite_num, float(unite_num) / min(len(setA), len(setB))

def Union(setA, setB):
    #�󼯺�A��B�Ĳ�����Ԫ�ظ�����ͬʱ���ؽ���ռС���ϵı���
    return len(setA.union(setB))

def Sub(setA, setB):
    #������A�е��ǲ���B�е�Ԫ�ظ���
    return len(setA.difference(setB))

def GapeSeq(Filename):
    """��ȡgape��seq�������"""
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(0, len(f11)):
        if f11[i, 0] == '>':
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

def test_sequence(Filename):
    """������׺Ϊ.spectra���ļ�����ȡ�����У���sequence���ϣ�,����set�з���"""
    Filename += "\\pFind-Filtered.spectra"
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(1, len(f11)):
        setA.add(f11[i].split('\t', 6)[5])
    return setA

def test_peptides(Filename):
    """������׺Ϊ.spectra���ļ�����ȡ������+��ʮ�У���peptide���ϣ�,����set�з���"""
    Filename += "\\pFind-Filtered.spectra"
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(1, len(f11)):
        setA.add(f11[i].split('\t', 12)[5] + f11[i].split('\t', 12)[10])
    return setA

def protein_group(Filename):
    """������׺Ϊ.protein���ļ�����ȡ����Ҫ���protein_group��������set�з���"""
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
    """�������ܣ���.spectra�ļ��а���nan�ļ��ĸ�ʽ��ȡpsm"""
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
    """�������ܣ�������׺Ϊ.protein���ļ�����ȡ���ն�protein��������set�з���"""
    Filename += "\\pFind.protein"
    with open(Filename) as f1:
        f11 = f1.readlines()
    setA = set()
    for i in range(2, len(f11)):
        if f11[i].split('\t', 3)[0].isdigit():   #�ҵ�һ��protein_group
            if str.lower(f11[i].split('\t', 3)[1][0 : 3]) != "rev":
                setA.add(f11[i].split('\t', 3)[1])
            j = i + 1   #��ʼ������һ��
            while f11[j].split('\t', 3)[1] == "SubSet":  #������һgroup��subset
                if str.lower(f11[j].split('\t', 3)[2][0 : 3]) != "rev":
                    setA.add(f11[j].split('\t', 3)[2])
                j += 1
            i = j
        elif "--" in f11[i].split('\t', 3)[0]:
            break
    return setA
    
    
    
    
    
def compare(dir1, dir2):
    """����dir1��dir2�ļ����µ�spectra�ļ���������У��ļ���֮�䣩sequence�Ľ��������
    �͵�����+��ʮ�� peptide �Ľ��������
    protein_group �Ľ��������"""
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
    """Ϊ�ܻ����ݼ���������
    ���ļ���dirname��sequence��peptide��������protein_group����"""
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
    """�������ܣ���list_info�и�Ԫ�ؼ����໥�󽻲���
    dir_list��ʾlist_info�и�Ԫ�ص���Դ"""
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
    pFind_file = "C:\\Users\\kfwang\\Desktop\\pt1_restricted_gluC\\pFind-Filtered.spectra"
    gape_file = "C:\\Users\\kfwang\\Desktop\\pt1_All_seq.txt"
    gape_seq = GapeSeq(gape_file)
    pFind_seq = PfindSeq(pFind_file)
    print("gape_seq ::  ", len(gape_seq))
    print("pFind_seq ::  ", len(pFind_seq))
    print(len(Unite(gape_seq, pFind_seq))