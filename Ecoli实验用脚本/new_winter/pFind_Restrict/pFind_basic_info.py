#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""pFind结果基本信息
PSM数量、peptide数量、PSM与standard的交集数量"""
import os
def get_PSM_pep(SP_file):
    """读取SP_file中报告的PSM信息: file_name(0)、scan_num(1)、seq(5)、mod(10)""" 
    with open(SP_file, 'r') as f:
        ff = f.readlines()
    psms = set()
    peps = set()
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        psm = segs[0].strip() + segs[1].strip() + segs[5].strip() + segs[10].strip()
        psms.add(psm)
        peps.add(segs[5].strip() + segs[10].strip())

    return psms, peps

def get_time(path):
    summary_file = path + "\\pFind.summary"
    a_file = path + "\\1.aa"
    time = (os.path.getmtime(summary_file) - os.path.getmtime(a_file)) / 60 
    return time 


if __name__ == "__main__":
    standard_file = "D:\\users\\kfwang\\two-step\\Ecoli\\standard\\pFind_R\\pFind-Filtered.spectra"
    test_file = "D:\\users\\kfwang\\two-step\\Ecoli\\Ecoli-reviewed\\one-step\\pFind_old\\pFind-Filtered.spectra"
    psm_s, pep_s = get_PSM_pep(standard_file)
    psm_t, pep_t = get_PSM_pep(test_file)
    time = get_time(test_file.rsplit('\\', 1)[0])
    print("total_standard\n", len(psm_s), len(pep_s))
    print("total_test\n", len(psm_t), len(pep_t))
    print("union\n", len(psm_t & psm_s), len(pep_t & pep_s))
    print("time\n", time)

        





