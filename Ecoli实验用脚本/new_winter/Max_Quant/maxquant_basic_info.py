#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:18:19 2021

@author: kaifeiwang
"""
"""max_quant结果基本信息
PSM数量、peptide数量、PSM与standard的交集数量"""

def get_PSM_pep(MQ_file):
    """读取MQ_file中报告的PSM信息: file_name(0)、scan_num(1)、modified_pep(7)""" 
    with open(MQ_file, 'r') as f:
        ff = f.readlines()
    psms = set()
    peps = set()
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        psm = segs[0].strip() + segs[1].strip() + segs[7].strip()
        psms.add(psm)
        peps.add(segs[7].strip())
    return psms, peps

if __name__ == "__main__":
    standard_file = "C:\\Users\\kfwang\\gygi\\standard\\MaxQuant\\combined\\txt\\msms.txt"
    test_file = "C:\\Users\\kfwang\\gygi\\one-step\\MaxQuant\\combined\\txt\\msms.txt"
    psm_s, pep_s = get_PSM_pep(standard_file)
    psm_t, pep_t = get_PSM_pep(test_file)
    print("total_standard\n", len(psm_s), len(pep_s))
    print("total_test\n", len(psm_t), len(pep_t))
    print("union\n", len(psm_t & psm_s), len(pep_t & pep_s))

        





