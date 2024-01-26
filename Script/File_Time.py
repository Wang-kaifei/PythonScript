'''
Descripttion: 
version: 
Author: Kaifei
Date: 2022-06-07 19:27:36
LastEditors: Kaifei
LastEditTime: 2023-06-19 11:45:31
'''
# -*- coding: utf-8 -*-
import os
import time
import re
import sys
 
def get_FileSize(filePath):
    fsize = os.path.getsize(filePath)
    fsize = fsize/float(1024*1024)
    return round(fsize,2)

def GetFile(path, time_th):
    res = []
    times = time_th.split('-')
    for root, _, files in os.walk(path):
        for file in files:
            full_path = os.path.join(root, file)
            try:
                fsize = get_FileSize(full_path)
                mtime = os.stat(full_path).st_mtime
                file_modify_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(mtime))
                segs = re.split('-| ', file_modify_time)
                if int(segs[0]) >= int(times[0]) and int(segs[1]) >= int(times[1]) and int(segs[2]) >= int(times[2]):
                    res.append((full_path, file_modify_time, fsize))
            except:
                print(f"Find file {full_path} error!")
    return res

def CMP(a):
    return -a[2]

def WriteRes(out_path, res):
    res = sorted(res, key = CMP)
    with open(out_path, "w") as f:
        f.write("File path\tModify time\tSize(MB)\n")
        for line in res:
            f.write("\t".join(str(i) for i in line))
            f.write('\n')
            
import ahocorasick


# 定义一组字符串
strings = ['hello', 'world', 'python', 'programming']

# 创建AC自动机并将字符串添加到树中
ac = ahocorasick.Automaton()
for string in strings:
    ac.add_word(string, string)
ac.make_automaton()

# 要查找的文本串
text = 'hello pythonasdvvvvprogramming'

# 在AC自动机中查找所有匹配项
matches = set()
for end_index, string in ac.iter(text):
    matches.add(string)

# 判断是否存在匹配项
if matches:
    print('文本串中存在以下子串：', matches)
else:
    print('文本串中不存在任何子串。')
#     if (len(sys.argv) != 4):
#         print("Error: Wrong param number!")
#         exit(0) 
#     test_path = sys.argv[1] #要检查的目录
#     out_path = sys.argv[2] #输出目录
#     time_th = sys.argv[3] #时间阈值，格式为xxxx-xx-xx
#     WriteRes(out_path, GetFile(test_path, time_th))
