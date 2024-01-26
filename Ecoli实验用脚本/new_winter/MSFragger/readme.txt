MSFragger建库：
分别用纯sequence字符串匹配、
peptide.tsv文件查找蛋白质匹配、
psm.tsv文件查找蛋白质匹配、
protein.tsv文件直接匹配（舍弃，protein.tsv文件报告的是protein group
	最后选用MSFragger_build_fasta_peptide.py脚本，其结果与直接字符串匹配达到一致。