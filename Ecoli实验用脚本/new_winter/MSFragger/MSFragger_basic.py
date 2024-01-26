# -*- coding: cp936 -*-
import os

"""�ű����ܣ�MSFragger������Ķ�����ʱ��"""

def get_pro_group(path):
    print(path)
    with open(path + "/protein.tsv") as f1:
        f11 = f1.readlines()
    print(len(f11))
    
def get_basic_info(path):
    """�������ܣ���ñ���ʵ���������������
    Ϊ����д���ļ����㣬�����ַ���������Ϣ��\t��Ϊ�ָ���
    ����Filename��λ�ļ��У�����ȡ��·����1.aa�ļ����޸�ʱ�䣬�Դ�������ʱ���"""

    #��������ʱ�䣬��ȷ������
    end_time = os.path.getmtime(path + "/modifications.tsv")
    start_time = os.path.getmtime(path + "/fragger.params")
    time = str((end_time - start_time) / 60)
    pep = 0
    pro = 0
    psm = 0
    wrong = set()
    protein = set()
    right = set()
    con = set()
    with open(path + "/peptide.tsv") as f1:
        f11 = f1.readlines()
        pep = len(f11) - 1
        for i in range(1, len(f11)):
            segs = f11[i].strip().split('\t')[-1]
            proteins = segs.split(',')
            if len(proteins) == 1:
                proteins = [segs]
            for pro in proteins:
                pro = pro.strip()
                if 'REV' == pro[0:3]:
                    continue
                protein.add(pro)
                if 'sp' == pro[0:2]:
                    wrong.add(pro)
                if 'gi' == pro[0:2]:
                    right.add(pro)
                if 'CON' == pro[0 : 3]:
                    con.add(pro)
    pro = len(protein)
    with open(path + "/psm.tsv") as f1:
        f11 = f1.readlines()
        psm = len(f11) - 1
        f1.close()
    print(path)
    print(len(con))
    print(len(wrong))
    print(len(right))
    print(time, pep, pro, psm)
    return time, pep, pro, psm

if __name__ == "__main__":
    path_list = [
       "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database1_comp",
                "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_comp",
       "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_comp",
       "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_comp_demo",
                "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_two_comp",
                "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2",
                "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\database2_0.05",
                "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\pFind_filter\\0.001",
          "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\pFind_filter\\0.01",
       "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\pFind_filter\\0.005",
        "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\pFind_filter_pro\\0.01_0.05",
       "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\pFind_database_pro\\spectra+pro\\0.01_0.03",
       "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\MSFragger\\pFind_database_pro\\pFind.protein\\0.2",
            ]
    
    for path in path_list:
        get_basic_info(path)

