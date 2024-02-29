'''
Descripttion: 
version: 
Author: Kaifei
Date: 2024-02-28 18:58:10
LastEditors: Kaifei
LastEditTime: 2024-02-28 19:53:27
'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""MSGF搜索大库时需要分割库分别搜索，得到的文件夹是按照库来组织的。但合并时，需要按照raw来重新组织"""

import os
import shutil

# 定义原始文件夹路径和新文件夹路径
original_folder_path = 'E:\\FTPtrans\\MSGF\\msgfres\\MSGF_res\\'
new_folder_path = 'E:\\FTPtrans\\MSGF\\msgfres\\merge\\'

# 创建一个字典用于存储文件名相同的文件
file_dict = {}

# 遍历所有文件夹
for folder in os.listdir(original_folder_path):
    folder_path = os.path.join(original_folder_path, folder)
    if os.path.isdir(folder_path):
        # 遍历文件夹中的文件
        for file in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file)
            if os.path.isfile(file_path):
                # 获取文件名（不包含扩展名）
                filename, file_extension = os.path.splitext(file)
                # 将相同文件名的文件加入到字典中
                if filename in file_dict:
                    file_dict[filename].append(file_path)
                else:
                    file_dict[filename] = [file_path]

# 创建新的文件夹存放重组后的文件
os.makedirs(new_folder_path, exist_ok=True)

# 遍历字典中的文件，将文件移动到新的文件夹中并进行重命名
for key, files_list in file_dict.items():
    new_folder = os.path.join(new_folder_path, key)
    os.makedirs(new_folder, exist_ok=True)
    # 移动文件并重命名
    for i, file_path in enumerate(files_list):
        new_file_name = f"{key}_{i}{os.path.splitext(file_path)[-1]}.mzid"
        new_file_path = os.path.join(new_folder, new_file_name)
        shutil.move(file_path, new_file_path)

print("File reorganization completed!")