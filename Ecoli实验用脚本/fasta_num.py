# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

def get_pronum(filepath):
    """函数功能：获取fasta文件中蛋白质的个数"""
    res = 0
    gi = 0
    con = 0
    with open(filepath) as file:
        file = file.readlines()
    for i in range(0, len(file)):  #遍历每一行
        if file[i][0 : 1] == ">":  #如果该行是蛋白质名
            res += 1
        if file[i][0:3].lower() == '>gi':
            gi += 1
        if file[i][0:4].lower() == '>con':
            con += 1
    print(res, gi, con)

def GetAAnum(filepath):
    """函数功能：获取fasta文件中氨基酸/碱基数量"""
    res = 0
    with open(filepath) as file:
        file = file.readlines()
    for i in range(0, len(file)):  #遍历每一行
        if file[i][0 : 1] == ">":  #如果该行是蛋白质名
            continue
        res += len(file[i].strip())
    print(res)



"""if __name__ == "__main__":
     fasta_path = [["C:\\test_wkf\\two_step_test\\limit1\\new_0.01.fasta", "C:\\test_wkf\\two_step_test\\limit1\\new_0.05.fasta", "C:\\test_wkf\\two_step_test\\limit1\\new_1.fasta"],
                    ["C:\\test_wkf\\two_step_test\\limit2\\new_0.01.fasta", "C:\\test_wkf\\two_step_test\\limit2\\new_0.05.fasta", "C:\\test_wkf\\two_step_test\\limit2\\new_1.fasta"],
                    ["C:\\test_wkf\\two_step_test\\open1\\new_0.01.fasta", "C:\\test_wkf\\two_step_test\\open1\\new_0.05.fasta", "C:\\test_wkf\\two_step_test\\open1\\new_1.fasta"],
                    ["C:\\test_wkf\\two_step_test\\open2\\new_0.01.fasta", "C:\\test_wkf\\two_step_test\\open2\\new_0.05.fasta", "C:\\test_wkf\\two_step_test\\open2\\new_1.fasta"]
                    ]
     for i in range(0, 4):
         for j in range(0, 3):
             get_pronum(fasta_path[i][j])"""
          
if __name__ == "__main__":
    fasta_path = ["C:\\Users\\pFind\\Desktop\\dl_model_data\\yeast\\trans_res\\database_td.fasta", "C:\\Users\\pFind\\Desktop\\dl_model_data\\yeast\\trans_res\\target.fasta", "C:\\Users\\pFind\\Desktop\\dl_model_data\\yeast\\trans_res\\decoy.fasta"]

    for path in fasta_path:
        print (path)
        get_pronum(path)
        GetAAnum(path)
                  
