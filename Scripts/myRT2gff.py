#!/usr/bin/env python
# Last update March 4 2021 by Fatemeh Sharifi

import sys
import re
import os

if len(sys.argv) > 2:
    DIR = sys.argv[2]
    print(DIR)
    ACC = sys.argv[1]
file_info_dictionary = {
    "gn_dom_txt": "-gn-dom.txt",
    "gn_list": "-gn.list",
    "RVT_gff": "-RVT.gff",
    "gn_gff": "-gn.gff",
    "gn_dom_list": "-gn-dom.list",
    "len_file": ".fna.len",
    "fgs_gff": "-FGS.gff",
}
file_paths_dictionary = dict()
for file_type in file_info_dictionary:
    file_paths_dictionary[file_type] = os.path.join(DIR, ACC + file_info_dictionary[file_type])

# gn: genomic neighborhood
# dom : domain
# RT/RVT: reverse transcriptase

gn_dom_list_file = open(file_paths_dictionary["gn_dom_list"], "w")
gn_list_file = open(file_paths_dictionary["gn_list"], "r")

#domains
dom = {}
#Genomic neighboorhood of CRISPR RTs
cas_gn = {}
# description of genes
des={}
RTs = list()
gn_list = list()
seq_seen = list()
gene_seen = list()
# CRISPR-Cas associated domains found in the neighborhood of CRISPR-RTs
matches = ["cas", "cmr", "csx", "Cmr", "Cas", "Csx"]
for aline in gn_list_file:
    casdom = 0
    if aline[0] == '#':
        continue
    subs = aline.strip().split(",")
    domlist = []
    RT_line = []
    for i in range(len(subs)):
        gn_dom_txt = open(file_paths_dictionary["gn_dom_txt"], "r")
        for domline in gn_dom_txt:
            par = domline.strip().split()
            if subs[i] == par[0] and par[1] != "-":
                domlist.append(par[1])
                if par[0] not in gn_list:
                    gn_list.append(par[0])
                if any(x in par[1] for x in matches):
                    casdom = casdom + 1
                if "RVT-" in par[1]:
                    # RT_line captures all RTs in one line of domfile
                    RT_line.append(par[0])
                    # RT captures all RTs
                    RTs.append(par[0])
    for j in RT_line:
        idx = RTs.index(j)
        cas_gn[idx] = casdom
    print >> gn_dom_list_file, ','.join(map(str, domlist))
    gn_dom_txt.close()
gn_list_file.close()
gn_dom_list_file.close()
# print RTs
des = {}
fgs_gff = open(file_paths_dictionary["fgs_gff"], "r")
for aline in fgs_gff:
    subs = aline.strip().split("\t")
    if aline[0] == '#':
        continue
    elif subs[2] == "region":
        continue
    else:
        subs2 = subs[-1].split(";")
        flagID = 0
        for ades in subs2:
            if (ades[:3] == 'ID=') and ades[3:] in gn_list:
                flagID = 1
                ourgn = ades[3:]
                idx = gn_list.index(ourgn)
                descr = ""
            elif flagID == 1 and (ades[:5] == 'gene=') and ades[5:] != "unk" and ades[5:] != "unknown":
                descr = descr + ades + ";"
            elif flagID == 1 and (ades[:4] == "Ori=") and ades[4:] != "unk" and ades[4:] != "unknown":
                descr = descr + ades
            if flagID == 1:
                des[idx] = descr
fgs_gff.close()
gff_file = open(file_paths_dictionary["gn_gff"], "w")
rt_file = open(file_paths_dictionary["RVT_gff"], "w")
gn_dom_txt = open(file_paths_dictionary["gn_dom_txt"], "r")
for aline in gn_dom_txt:
    if aline[0] == '#':
        continue
    subs = aline.strip().split()
    seqid = re.sub(r"_[0-9]*_[0-9]*_[-|+]", "", subs[0])
    par = subs[0][len(seqid) + 1:].split("_")
    start_gene = int(par[0])
    stop_gene = int(par[1])
    strand = par[2]
    if seqid not in seq_seen:
        seq_seen.append(seqid)
        len_file = open(file_paths_dictionary["len_file"], "r")
        for l in len_file:
            col = l.strip().split()
            if seqid == col[1]:
                len_seq = col[0]
                output = seqid + "\tGenbank\tregion\t1\t" + len_seq + "\t.\t" + strand + "\t.\tID=" + seqid
                print >> gff_file, output
        len_file.close()
    if subs[0] not in gene_seen:
        gene_seen.append(subs[0])
        if subs[0] in gn_list:
            print gn_list
            print des
            idx = gn_list.index(subs[0])
            descr = ";" + des[idx]
        else:
            descr = ""
        output = seqid + "\tFGS\tCDS\t" + str(start_gene) + "\t" + str(stop_gene) + "\t.\t" + strand + "\t0\tID=" + subs[
            0] + ";what=gene" + descr
        print >> gff_file, output
    dom_start = start_gene + (int(subs[2]) * 3)
    dom_end = start_gene + (int(subs[3]) * 3)
    if "RVT-" in subs[1]:
        if "RVT-CRISPR" == subs[1]:
            if cas_gn[RTs.index(subs[0])] == 0:
                subs[1] = "RVT-CRISPR-like2"
        # elif "RVT-CRISPR-like" in subs[1]:
        # if cas_gn[RTs.index(subs[0])] > 0:
        # subs[1]="RVT-CRISPR"
        output = seqid + "\tmyRT\tdom\t" + str(dom_start) + "\t" + str(
            dom_end) + "\t.\t" + strand + "\t0\tID=" + seqid + "_" + str(dom_start) + "_" + str(
            dom_end) + "_" + strand + ";what=RT;des=" + subs[1]
        print >> gff_file, output
        out = DIR + "/" + ACC + "\t" + seqid + "\tmyRT\t" + seqid + "_" + str(start_gene) + "_" + str(
            stop_gene) + "_" + strand + "\t" + subs[1]
        print >> rt_file, out
    elif subs[1] != "-":
        output = seqid + "\tmyRT\tdom\t" + str(dom_start) + "\t" + str(
            dom_end) + "\t.\t" + strand + "\t0\tID=" + seqid + "_" + str(dom_start) + "_" + str(
            dom_end) + "_" + strand + ";what=dom;des=" + subs[1]
        print >> gff_file, output
rt_file.close()
gff_file.close()
gn_dom_txt.close()

