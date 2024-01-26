import os
def get_basic_info(Filename):
    """�������ܣ���ñ���ʵ���������������
    Ϊ����д���ļ����㣬�����ַ���������Ϣ��\t��Ϊ�ָ���
    ����Filename��λ�ļ��У�����ȡ��·����1.aa�ļ����޸�ʱ�䣬�Դ�������ʱ���"""
    results = ""
    with open(Filename) as f1:
        f11 = f1.readlines()
    #��ͷ����Щ��Ϣ
    for i in range(1, 7):
        results += f11[i].split()[-1] + '\t'
    #��ȡIDrate
    for i in range(-1, -100, -1):
        if f11[i][0] == "I" and f11[i][1] == "D":
            j = i + 1
            while j <= -1 and f11[j][0] != "-":
                results += f11[j].split()[-1] + '\t'
                j += 1
            break
    #��ȡMixed MS2
    for i in range(-1, -100, -1):
        if f11[i][0] == "M" and f11[i][1] == "i" and f11[i][2] == "x":
            j = i + 1
            while j <= -1 and f11[j][0] != "-" and j - i < 6:
                results += f11[j].split()[-1] + '\t'
                j+= 1
            break
    #��������ʱ�䣬��ȷ������
    summary_time = os.path.getmtime(Filename)
    aa_name = Filename.rsplit("\\", 1)[0] + "\\1.aa"
    aa_time = os.path.getmtime(aa_name)
    results += str((summary_time - aa_time) / 60)
    return results

def file_name(file_dir):
    """�������ܣ���ȡfile_dir�ļ����µ�����pFind.summary�ļ�������"""
    L = []
    for dirpath, dirnames, filenames in os.walk(file_dir):
        for file in filenames:
            if file == 'pFind.summary':
                L.append(os.path.join(dirpath, file))
    return L

def write_info(infos, filename):
    """��������Ϣд��filename�ļ�"""
    title = "Name\tSpectra\tScans\tPeptides\tSequences\tProteins\tProtein Groups\tID Rate 1\tID Rate 2\tID Rate 3\tID Rate all\tMixed MS/MS 1\tMixed MS/MS 2\tMixed MS/MS 3\tMixed MS/MS 4\tMixed MS/MS 5\tTimes\n"
    with open(filename, "w") as f:
        f.write(title)
        for info in infos:
            f.write(info)
            
def get_sum(filename):
    """�˺���Ϊ�ܻ����ݴ�������
    �������ܣ��ǽ����ֱ�������������
    û��ȥ�أ���������  ����sequence_pro.py�еĺ������д���"""
    results = []
    with open(filename) as f1:
        f11 = f1.readlines()
    i = 1
    while i < len(f11):
        if "\\" in f11[i].split('\t', 1)[0]:  #����ǻ������ļ��е�����
            result = f11[i].split('\t', 1)[0].split("\\")[0]  #���ļ�����
            list1 = f11[i].split()
            list2 = f11[i + 1].split()
            list3 = f11[i + 2].split()
            for j in range(1, len(list1)):
                if "%" in list1[j]:
                    result += "\t" + str((
                        float(list1[j].split("%")[0]) + float(list2[j].split("%")[0]) + float(list3[j].split("%")[0])) / 3) + "%"
                else:
                    result += "\t" + str(float(list1[j]) + float(list2[j]) + float(list3[j]))
            result += "\n"
            results.append(result)
            i += 3
        i += 1
    with open(filename, "a") as f:
        for result in results:
            f.write(result)

"""if __name__ == "__main__":
    path = "C:\\test_wkf\\two_step_test\\after_two_step\\"
    path_list = file_name(path)  #��ȡ���ļ����µ�����summary�ļ�·��
    infos  = []
    for Filename in path_list:
        #������ȡ�ļ�������path�°������ļ�������
        info = str(Filename.split("\\")[-2]) + "\t"
#        if len(Filename.split("\\")) == 4:
#            info = str(Filename.split("\\")[2]) + "\t"
#        else:
#            info = str(Filename.split("\\")[2]) + "\\" + str(Filename.split("\\")[3]) + "\t"
        #��ȡ������Ϣ
        info += get_basic_info(Filename) + '\n'
        infos.append(info)

    File = "C:\\test_wkf\\two_step_test\\after_two_step\\tttt.txt"
    write_info(infos, File)
#    get_sum(File)
"""

if __name__ == "__main__":
    """path_list = ["C:\\test_wkf\\two_step_test\\TTE\\OO\\OO_0.01\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\OO\\OO_0.05\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\OO\\OO_1\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\OR\\OR_0.01\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\OR\\OR_0.05\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\OR\\OR_1\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\RR\\RR_0.01\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\RR\\RR_0.05\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\RR\\RR_1\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\RO\\RO_0.01\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\RO\\RO_0.05\\pFind.summary",
                 "C:\\test_wkf\\two_step_test\\TTE\\RO\\RO_1\\pFind.summary"
                 ]"""
    path_list = [
                "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\big\\O_mod_100\\OR_0.001\\pFind.summary",
        
            ]
    
    infos  = []
    for Filename in path_list:
        #������ȡ�ļ�������path�°������ļ�������
        info = str(Filename.split("\\")[-2]) + "\t"
        #��ȡ������Ϣ
        info += get_basic_info(Filename) + '\n'
        infos.append(info)

    File = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\big\\O_mod_100\\OR_0.001\\basic.txt"
    write_info(infos, File)
    
