#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""MSGF结果基本信息
PSM数量、peptide数量、PSM与standard的交集数量"""

def get_PSM_pep(comet_file):
    """读取MSGF中报告的PSM信息: file(0)、scan_num(1)、modified_pep(9)""" 
    with open(comet_file, 'r') as f:
        ff = f.readlines()
    psms = set()
    peps = set()
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        psm = segs[0].strip() + segs[1].strip() + segs[9].strip()
        psms.add(psm)
        peps.add(segs[9].strip())
    return psms, peps

if __name__ == "__main__":
    standard_file = "D:\\users\\kfwang\\two-step\\Ecoli\\standard\\MS-GF+\\tsv\\all_result_filter.tsv"
    test_file = "C:\\Users\\kfwang\\results\\griffin\\two-step\\pFind\\0.1\\MS-GF+\\tsv\\all_result_filter.tsv"
    #psm_s, pep_s = get_PSM_pep(standard_file)
    psm_t, pep_t = get_PSM_pep(test_file)
    #print("total_standard", len(psm_s), len(pep_s))
    print("total_test", len(psm_t), len(pep_t))
    #print("union", len(psm_t & psm_s), len(pep_t & pep_s))

        





