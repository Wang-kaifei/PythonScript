#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 12:15:00 2020

@author: kaifeiwang
"""
import os 
def get_build_time(path):
    """函数功能：获取进行第一次搜索的时间
    输入：.summary文件的路径"""
    summary_time = os.path.getmtime(path)
    aa_name = path.rsplit("\\", 1)[0] + "\\1.aa"
    aa_time = os.path.getmtime(aa_name)
    time = str((summary_time - aa_time) / 60)
    
    print(time)
    return time


if __name__ == "__main__":
    pathlist = ["C:\\test_wkf\\two_step_test\\limit1\\pFind.summary",
                "C:\\test_wkf\\two_step_test\\limit2\\pFind.summary",
                "C:\\test_wkf\\two_step_test\\open1\\pFind.summary",
                "C:\\test_wkf\\two_step_test\\open2\\pFind.summary"]
    
    for path in pathlist:
        get_build_time(path)
