#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""




def get_names (fasta_path):

    res = 0
    with open(fasta_path) as file:
        file = file.readlines()
    print(len(file))
    i = 0
    while(i < len(file)):
        if file[i][0] == ">":  #如果该行是蛋白质名
            if file[i][1] == "C":
                print(file[i])
            res += 1
        i += 1

    print(res)
        


if __name__ == "__main__":
    fasta_path = "C:\\Users\\kfwang\\gygi\\database\\uniprot_reviewed_human_+maxquant.fasta"
    #outfile = "C:\\Users\\kfwang\\test_wkf\\two_step_test\\new_10,000\\uniprot\\name.fasta"
    get_names (fasta_path)





