'''
Descripttion: 处理大mgf文件的情况  删除谱峰数 <= 6 的谱图，将谱图文件切分为多个新的子文件
version: 
Author: kfwang
Date: 2022-03-27 15:50:22
LastEditors: sueRimn
LastEditTime: 2022-04-02 15:22:21
'''
# -*- coding: utf-8 -*-
# @author Wang kaifei

def is_number(s):
    try:  # 如果能运行float(s)语句，返回True（字符串s是浮点数）
        float(s)
        return True
    except ValueError:  # ValueError为Python的一种标准异常，表示"传入无效的参数"
        pass  # 如果引发了ValueError这种异常，不做任何事情（pass：不做任何事情，一般用做占位语句）
    try:
        import unicodedata  # 处理ASCii码的包
        unicodedata.numeric(s)  # 把一个表示数字的字符串转换为浮点数返回的函数
        return True
    except (TypeError, ValueError):
        pass
    return False


def TitleLegal(mgf_path, out_path):
    with open(mgf_path) as f:
        lines = f.readlines()
    with open(out_path, 'w') as f:
        i = 0
        for line in lines:
            i += 1
            if (i % 1000 == 0):
                print(i)
            if line[0:5] == "TITLE":
                tmp = line.strip().split(" ")[0]
                tmp = tmp + ".dta\n"
                f.write(tmp)
            else:
                f.write(line)

            
def IfWrite(msms: list):
    cnt = 0
    for line in msms:
        if (is_number(line.split(' ')[0])):
            cnt += 1
    return cnt >= 7
            

def DelShortandDepart(mgf_path, out_path, num):
    with open(mgf_path) as f:
        lines = f.readlines()
    pace = len(lines) / num
    for i in range(10):
        f = open(out_path + str(num) + ".mgf", "w")
        msms = [] # 临时存储一张谱图
        for j in range(i * pace, min((i + 1) * pace, len(lines))):
            if lines[j][0:5] == "BEGIN":
                if IfWrite(msms): # 处理上一张谱图
                    for l in msms:
                        f.write(l)
                msms = [lines[j]]
            else:
                msms.append(lines[j])
        if IfWrite(msms): # 处理最后一张谱图
            for l in msms:
                f.write(l)
        f.close()

            
if __name__ == "__main__":
    mgf_path = "C:\\Users\\pFind\\Desktop\\ymk_mgf\\ZCQ_YMK-CR_Slot1-7_1_11663.mgf"
    out_path = "C:\\Users\\pFind\\Desktop\\ymk_mgf\\"
    DelShortandDepart(mgf_path, out_path, 10)
