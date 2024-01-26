#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""



def delete (fasta_path, outfile):

    res = ""
    with open(fasta_path) as file:
        file = file.readlines()
    print(len(file))
    i = 0
    while i < len(file):
        if i% 100000 == 0:
            print(i)
        if file[i][0] == ">" and file[i][1 : 4].lower() == "con":  #如果该行是蛋白质
            i += 1
            while i < len(file) and file[i][0] != '>':
                i += 1
        else:
            res += file[i]
            i += 1

    with open(outfile, "w") as f:
        f.write(res)

def show(fasta_path):
    with open(fasta_path) as file:
        file = file.readlines()
    for i in range(len(file) - 100, len(file)):
        print(file[i])


if __name__ == "__main__":
    fasta_path = "C:\\Users\\kfwang\\gygi\\one-step\\pFind-no-coelute\\database\\0.1.fasta"
    outfile = "C:\\Users\\kfwang\\gygi\\one-step\\pFind-no-coelute\\database\\0.1_maxquant.fasta"
    delete (fasta_path, outfile)
    #show(fasta_path)





