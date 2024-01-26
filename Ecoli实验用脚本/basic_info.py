import os
def get_basic_info(Filename):
    """函数功能：获得本次实验结果的最基本数据
    为最终写入文件方便，按照字符串生成信息，\t作为分隔符
    根据Filename定位文件夹，并获取该路径下1.aa文件的修改时间，以此来计算时间差"""
    results = ""
    with open(Filename) as f1:
        f11 = f1.readlines()
    #开头的那些信息
    for i in range(1, 7):
        results += f11[i].split()[-1] + '\t'
    #提取IDrate
    for i in range(-1, -100, -1):
        if f11[i][0] == "I" and f11[i][1] == "D":
            j = i + 1
            while j <= -1 and f11[j][0] != "-":
                results += f11[j].split()[-1] + '\t'
                j += 1
            break
    #提取Mixed MS2
    for i in range(-1, -100, -1):
        if f11[i][0] == "M" and f11[i][1] == "i" and f11[i][2] == "x":
            j = i + 1
            while j <= -1 and f11[j][0] != "-" and j - i < 6:
                results += f11[j].split()[-1] + '\t'
                j+= 1
            break
    #计算所需时间，精确到分钟
    summary_time = os.path.getmtime(Filename)
    aa_name = Filename.rsplit("\\", 1)[0] + "\\1.aa"
    aa_time = os.path.getmtime(aa_name)
    results += str((summary_time - aa_time) / 60)
    return results

def file_name(file_dir):
    """函数功能：提取file_dir文件夹下的所有pFind.summary文件名返回"""
    L = []
    for dirpath, dirnames, filenames in os.walk(file_dir):
        for file in filenames:
            if file == 'pFind.summary':
                L.append(os.path.join(dirpath, file))
    return L

def write_info(infos, filename):
    """将基本信息写入filename文件"""
    title = "Name\tSpectra\tScans\tPeptides\tSequences\tProteins\tProtein Groups\tID Rate 1\tID Rate 2\tID Rate 3\tID Rate all\tMixed MS/MS 1\tMixed MS/MS 2\tMixed MS/MS 3\tMixed MS/MS 4\tMixed MS/MS 5\tTimes\n"
    with open(filename, "w") as f:
        f.write(title)
        for info in infos:
            f.write(info)
            
def get_sum(filename):
    """此函数为周虎数据处理特有
    函数功能：是将结果直接相加求总体结果
    没有去重，并不能用  可用sequence_pro.py中的函数进行处理"""
    results = []
    with open(filename) as f1:
        f11 = f1.readlines()
    i = 1
    while i < len(f11):
        if "\\" in f11[i].split('\t', 1)[0]:  #如果是还有子文件夹的数据
            result = f11[i].split('\t', 1)[0].split("\\")[0]  #父文件夹名
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
    path_list = file_name(path)  #提取该文件夹下的所有summary文件路径
    infos  = []
    for Filename in path_list:
        #首先提取文件名（即path下包含的文件夹名）
        info = str(Filename.split("\\")[-2]) + "\t"
#        if len(Filename.split("\\")) == 4:
#            info = str(Filename.split("\\")[2]) + "\t"
#        else:
#            info = str(Filename.split("\\")[2]) + "\\" + str(Filename.split("\\")[3]) + "\t"
        #提取基本信息
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
        #首先提取文件名（即path下包含的文件夹名）
        info = str(Filename.split("\\")[-2]) + "\t"
        #提取基本信息
        info += get_basic_info(Filename) + '\n'
        infos.append(info)

    File = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new\\Ecoli\\big\\O_mod_100\\OR_0.001\\basic.txt"
    write_info(infos, File)
    
