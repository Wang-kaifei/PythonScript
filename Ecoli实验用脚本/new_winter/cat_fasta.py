#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang

将两个fasta文件拼接
"""




def add_Ecoli (fasta_path1, fasta_path2):
    with open(fasta_path2) as file:
        file = file.readlines()
    print(len(file))
    with open(fasta_path1) as file1:
        tmp = file1.readlines()
        print(len(tmp))
        file1.close()

    with open(fasta_path1, 'a') as file1:
        for line in file:
            file1.write(line)
        file1.close()

    with open(fasta_path1) as file1:
        tmp = file1.readlines()
        print(len(tmp))
        file1.close()
        


if __name__ == "__main__":
    fasta_path1 = "D:\\users\\kfwang\\two-step\\Gygi\\database\\uniprot_reviewed.fasta"
    fasta_path2 = "D:\\users\\kfwang\\two-step\\Gygi\\database\\uniprot-human-filtered-reviewed.fasta"
    add_Ecoli (fasta_path1, fasta_path2)





