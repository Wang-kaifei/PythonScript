# -*- coding: utf-8 -*-


from tqdm import tqdm

class pFind_PSM:
    def __init__(self, File_Name="", Scan_No="", Exp_MHplus="", Charge="", Q_value="", Sequence="", Calc_MHplus="",
                 Mass_Shift="", Raw_Score="", Final_Score="", Modification="", Specificity="", Proteins="",
                 Positions="", Label="", Target="", Miss_Clv_Sites="", Avg_Frag_Mass_Shift="", Others=""):
        self.File_Name = File_Name
        self.Scan_No = Scan_No
        self.Exp_MHplus = Exp_MHplus
        self.Charge = Charge
        self.Q_value = Q_value
        self.Sequence = Sequence
        self.Calc_MHplus = Calc_MHplus
        self.Mass_Shift = Mass_Shift
        self.Raw_Score = Raw_Score
        self.Final_Score = Final_Score
        self.Modification = Modification
        self.Specificity = Specificity
        self.Proteins = Proteins
        self.Positions = Positions
        self.Label = Label
        self.Target = Target
        self.Miss_Clv_Sites = Miss_Clv_Sites
        self.Avg_Frag_Mass_Shift = Avg_Frag_Mass_Shift
        self.Others = Others

    def to_string(self):
        return str(self.File_Name) + '\t' + str(self.Scan_No) + '\t' + str(self.Exp_MHplus) + '\t' + str(
            self.Charge) + '\t' + str(self.Q_value) + '\t' + str(self.Sequence) + '\t' + str(
            self.Calc_MHplus) + '\t' + str(self.Mass_Shift) + '\t' + str(self.Raw_Score) + '\t' + str(
            self.Final_Score) + '\t' + str(self.Modification) + '\t' + str(self.Specificity) + '\t' + str(
            self.Proteins) + '\t' + str(self.Positions) + '\t' + str(self.Label) + '\t' + str(self.Target) + '\t' + str(
            self.Miss_Clv_Sites) + '\t' + str(self.Avg_Frag_Mass_Shift) + '\t' + str(self.Others) + '\n'
        # return '\t'.join(('%s' % item for item in self.__dict__.values())) + '\n'




def comet_mod_to_pfind(mod_str):
    if mod_str == "-":
        return ""
    mods = mod_str.split(',')
    mod_dic = {"15.994915":"Oxidation[M]", "57.021464":"Carbamidomethyl[C]", "-17.026549":"Gln->pyro-Glu[AnyN-termQ]", "42.010565":"Acetyl[ProteinN-term]"}
    res = ""
    for mod in mods:
        segs = mod.split('_')
        res += segs[0] + "," + mod_dic[segs[2]] + ";"
    return res

def comet_lines_to_pfind(in_path):
    res = []
    with open(in_path, 'r') as f:
        ff = f.readlines()

    # 0:scan	1:num	2:charge	3:exp_neutral_mass	4:calc_neutral_mass	5:e-value	6:xcorr	7:delta_cn	8:sp_score	9:ions_matched	10:ions_total	11:plain_peptide	12:modified_peptide	13:prev_aa	14:next_aa	15:protein	16:protein_count	17:modifications	18:retention_time_sec	19:sp_rank	20:Target	21:File	22:q_value
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        spec_name = segs[21].rsplit('_', 1)[0] + "." + segs[0] + "." + segs[0] + "." + segs[2] + ".0.dta"
        mods_pfind = comet_mod_to_pfind(segs[17])
        res.append(
            pFind_PSM(spec_name, segs[0], segs[3], segs[2], "0", segs[11], segs[4], "0", "0", "0", mods_pfind, "3", segs[15], "-",
                      "-",
                      "target" if segs[20].upper() == "TRUE" else "decoy", "0", "0", "-0"))
    return res

def msgf_mod_seq(raw_pep):
    mod_dic = {"+15.995":"Oxidation[M]", "+57.021":"Carbamidomethyl[C]", "-17.027":"Gln->pyro-Glu[AnyN-termQ]", "+42.011":"Acetyl[ProteinN-term]"}
    # 如果最左端出现（N端）修饰，则index是0
    mods = ""
    seq = "" # pFind格式的序列
    aa_index = 0 # 目前有多少个氨基酸出现
    temp_mod = "" # 暂时存储的修饰
    mid = raw_pep.split('.', 1)[1].rsplit('.', 1)[0]
    for cc in mid:
        if cc.isalpha():
            if temp_mod != "": # 如果前面的氨基酸有修饰，需要存储
                mods += str(aa_index) + "," + mod_dic[temp_mod] + ";"
                temp_mod = ""
            seq += cc
            aa_index += 1
        elif cc == "+" or cc == "-":
            if temp_mod != "": # 如果前面的氨基酸有修饰，需要存储
                mods += str(aa_index) + "," + mod_dic[temp_mod] + ";"
            temp_mod = cc
        else:
            temp_mod += cc
    return seq, mods

"""注意MGF格式的scan需要一步转换"""
def msgf_get_scan_name(mgf_file_name, line_index):
    # title_file是谱图名文件，里面的每一行都是一个谱图的title
    title_file = mgf_folder + mgf_file_name.rsplit('.', 1)[0] + ".title"
    with open(title_file, 'r') as f:
        lines = f.readlines()
    spec_name = lines[line_index].split("=")[1].strip()
    return spec_name

def msgf_lines_to_pfind(in_path):
    res = []
    with open(in_path, 'r') as f:
        ff = f.readlines()
    # #Spec_file:0	Spec_ID:1	ScanNum:2	Scan_time:3	FragMethod:4	Precursor:5	IsotopeError:6	Precursor_Error:7	charge:8	Peptide:9	Protein:10	DeNovoScore:11	MSGFScore:12	SpecEvalue:13	Evalue:14	target:15	q_value:16
    for i in tqdm(range(1, len(ff))):
        segs = ff[i].strip().split('\t')
        spec_name = msgf_get_scan_name(segs[0], int(segs[2]))
        scan_num = spec_name.split(".")[-4]
        charge = spec_name.split(".")[-3]
        if charge != segs[8]:
            print(ff[i])
        seqence, mods = msgf_mod_seq(segs[9])
        res.append(
            pFind_PSM(spec_name, scan_num, "0", segs[8], "0", seqence, "0", "0", "0", "0", mods, "3", segs[10], "-",
                      "-", "target" if segs[15].upper() == "TRUE" else "decoy", "0", "0", "0"))

    return res

def maxquant_mod(sequence, raw_mods):
    mods_dict = {"Oxidation (M)":"Oxidation[M]", "Gln->pyro-Glu":"Gln->pyro-Glu[AnyN-termQ]", "Acetyl (Protein N-term)":"Acetyl[ProteinN-term]"}
    mods = ""
    mods_cnt = {}
    for i in range(len(sequence)):
        if sequence[i] == "C":
            mods += str(i + 1) + "," + "Carbamidomethyl[C];"
    if "Unmodified" in raw_mods:
        return mods
    for mod in raw_mods:
        cnt = 1 if mod[0].isalpha() else int(mod[0])
        mod_name = mod if cnt == 1 else mod.split(" ", 1)[1]
        if mod_name in mods_cnt:
            print(sequence, raw_mods)
        mods_cnt[mod_name] = cnt
    for key, value in mods_cnt.items():
        for i in range(value):
            if key not in mods_dict:
                print(f"not match: {key}")
            mods += str(1) + "," + mods_dict[key] + ";"
    return mods    

def maxquant_lines_to_pfind(in_path):
    res = []
    with open(in_path, 'r') as f:
        ff = f.readlines()

    # Raw file:0	Scan number:1	Scan index:2	Sequence:3	Length:4	Missed cleavages:5	Modifications:6	Modified sequence:7	Oxidation (M) Probabilities:8	Oxidation (M) Score diffs:9	Acetyl (Protein N-term):10	Gln->pyro-Glu:11	Oxidation (M):12	Proteins:13	Charge:14	Fragmentation:15	Mass analyzer:16	Type:17	Scan event number:18	Isotope index:19	m/z:20	Mass:21	Mass error [ppm]:22	Mass error [Da]:23	Simple mass error [ppm]:24	Retention time:25	PEP:26	Score:27	Delta score:28	Score diff:29	Localization prob:30	Combinatorics:31	PIF:32	Fraction of total spectrum:33	Base peak fraction:34	Precursor full scan number:35	Precursor Intensity:36	Precursor apex fraction:37	Precursor apex offset:38	Precursor apex offset time:39	Matches:40	Intensities:41	Mass deviations [Da]:42	Mass deviations [ppm]:43	Masses:44	Number of matches:45	Intensity coverage:46	Peak coverage:47	Neutral loss level:48	ETD identification type:49	Reverse:50	All scores:51	All sequences:52	All modified sequences:53	Reporter PIF:54	Reporter fraction:55	id:56	Protein group IDs:57	Peptide ID:58	Mod. peptide ID:59	Evidence ID:60	Oxidation (M) site IDs:61
    for i in range(1, len(ff)):
        segs = ff[i].strip().split('\t')
        spec_name = segs[0] + "." + segs[1] + "." + segs[1] + "." + segs[14] + ".0.dta"
        mods = maxquant_mod(segs[3], segs[6].split(","))
        res.append(
            pFind_PSM(spec_name, segs[1], "0", segs[14], "0", segs[3], "0", "0", "0", "0", mods, "3", segs[13], "-",
                      "-", "target", "0", "0", "0"))

    return res


trans_dict = {
    "comet": comet_lines_to_pfind,
    "msgf": msgf_lines_to_pfind,
    "maxquant": maxquant_lines_to_pfind
}


def to_pfind(source_type, in_path, out_path):
    title = "File_Name\tScan_No\tExp_MHplus+\tCharge\tQ_value\tSequence\tCalc_MHplus+\tMass_Shift(Exp.-Calc.)\tRaw_Score\tFinal_Score\tModification\tSpecificity\tProteins\tPositions\tLabel\tTarget\ttMiss_Clv_Sites\tAvg_Frag_Mass_Shift\tOthers\n"
    res = []
    fun = trans_dict.get(source_type, None)
    if fun is None:
        print("Invalid Type, pls use: comet, msgf or maxquant")
        return -1
    res.extend(fun(in_path))

    # 将结果写入文件
    with open(out_path, 'w') as f:
        f.write(title)
        for line in res:
            # f.write(line.strip())
            f.write(line.to_string())

    print("Transfered " +
          str(len(res)) +
          " line(s) successed.")
    return 0

if __name__ == "__main__":
    in_path = "/Users/kaifeiwang/Desktop/trans/origin-two/msms.txt"
    out_path = "/Users/kaifeiwang/Desktop/trans/origin-two/test.spectra"
    mgf_folder = "/Users/kaifeiwang/Desktop/mgf/"
    source_type = "maxquant"  # comet, msgf, maxquant
    protein = to_pfind(source_type, in_path, out_path)
