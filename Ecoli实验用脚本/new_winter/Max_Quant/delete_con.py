#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""给污染蛋白重命名"""



def delete (fasta_path, outfile):

    res = ""
    with open(fasta_path) as file:
        file = file.readlines()
    print(len(file))
    i = 0
    cnt = 0
    while i < len(file):
        if i% 100000 == 0:
            print(i)
        if file[i][0 : 4] == ">CON":  #如果该行是蛋白质
            print(file[i])
            res += ">CON|" + str(cnt) + "\n"
            cnt += 1
            i += 1
        else:
            res += file[i]
            i += 1

    with open(outfile, "w") as f:
        f.write(res)


if __name__ == "__main__":
    fasta_path = "C:\\Users\\kfwang\\gygi\\one-step\\pFind-no-coelute\\database\\0.01.fasta"
    outfile = "C:\\Users\\kfwang\\gygi\\one-step\\pFind-no-coelute\\database\\0.01_maxquant.fasta"
    delete (fasta_path, outfile)





