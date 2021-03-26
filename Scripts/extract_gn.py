#!/usr/bin/env python

# Last update March 4 2021 by Fatemeh Sharifi

import sys
import re
import os

if len(sys.argv) > 2:
    DIR = sys.argv[2]
    print(DIR)
    ACC = sys.argv[1]
eps = 0.00001
intergenic_distance = 2000
file_info_dictionary = {
    "RVT_tmp": "-RVT.tmp",
    "RVT_cdd_pfam": "-RVT-cdd-pfamA.domtblout",
    "RVT_domtblout": "-RVT.domtblout",
    "RVT_faa": "-RVT.faa",
    "RT_faa": "-RT.faa",
    "QUERY_faa": "-query.faa",
    "Genomic_neighborhood_file": "-gn.list",
    "Genomic_neighborhood_faa": "-gn.faa",
    "gff_file": "-FGS.gff",
    "faa_file": "-FGS.faa",
}
file_paths_dictionary = dict()
for file_type in file_info_dictionary:
    file_paths_dictionary[file_type] = os.path.join(DIR, ACC + file_info_dictionary[file_type])

# print file_paths_dictionary["RVT_tmp"]
RTlist = list()
RT = list()
OUT = list()
query = list()
verified = list()
# Verifies the RTs based on cdd/pfam hits
inf = open(file_paths_dictionary["RVT_cdd_pfam"], "r")
matches = ["RVT_1", "RVT_N", "Cas1", "Abi", "RT_", "cas1", "GIIM", "Cas6", "cas6", "Cas2", "cas2"]
for aline in inf:
    if aline[0] == '#':
        continue
    subs = aline.strip().split()
    if subs[3] not in RTlist:
        RTlist.append(subs[3])
        ln = 1
        flagV = 0
        if "RVT_1" in subs[0] or "RT_" in subs[0] or "RVT_N" in subs[0]:
            verified.append(subs[3])
            flagV = 2
        elif any(x in subs[0] for x in matches):
            flagV = 1
    elif flagV < 2 and ln <= 10:
        ln += 1
        if flagV < 2 and (ln <= 10) and any(x in subs[0] for x in matches):
            flagV += 1
            if flagV >= 2:
                verified.append(subs[3])
    else:
        continue
inf.close()
# print (verified)
if len(verified) > 0:
    reversetranscriptase_faa = open(file_paths_dictionary["RT_faa"], "w")
    for v in verified:
        inf = open(file_paths_dictionary["RVT_faa"], "r")
        for aline in inf:
            if aline[0] == '>' and v in aline:
                print >> reversetranscriptase_faa, aline, inf.next(),
                break
        inf.close()
    reversetranscriptase_faa.close()
# Since we don't know how many hits we have for each RT gene and we only need a maximum of three, we continue reading until we either get the three top hits or get to the next RT gene
# Then we check the flagS to see if it's 1. flagS==1 means we have processed the hits for the previous RT. 
flagS = 2
inf = open(file_paths_dictionary["RVT_tmp"], "r")
for aline in inf:
    # Since aline is a line from a file, len(aline) > 0
    if aline[0] == '#':
        continue
    subs = aline.strip().split()
    subs[22] = re.sub(r"CRISPR-G[0-9]", "CRISPR", subs[22])
    subs[22] = re.sub(r"GII-I[I]*", "GII", subs[22])
    # We only process verified RTs.
    if subs[3] in verified:
        # This means that we are processing the top hit for a new RT, so let's make sure we are done with the previous RT before moving to this one.
        if subs[3] not in RT:
            # flagS==0 means that previous RT has had less than three top hits.
            if flagS == 0:
                # If the previous RT has has only one hit just output it.
                if ln == 1:
                    firsthit[0] = first
                    firsthit.append(firsthit[6])
                # If the previous RT has had two top hits lets see if they belong to one class (e.g. GII-I and GII-II that both belong to GII) or two different classes.
                elif ln == 2 and float(float(firsthit[6]) / float(secondhit[6])) > eps:
                    if (first == second):
                        firsthit[0] = first
                    else:
                        firsthit[0] = first + "/" + second  # two different classes
                        query.append(firsthit[3])
                    firsthit.append(float(float(firsthit[6]) / float(secondhit[6])))
                else:
                    firsthit[0] = first
                    firsthit.append(float(float(firsthit[6]) / float(secondhit[6])))
                OUT.append(firsthit)
                # flagS=1
            # Let's process this new RT.
            RT.append(subs[3])
            flagS = 0
            ln = 1  # We are processing the 1st hit for this RT.
            # Let's get rid of subclass information for example let's write CRISPR-G1 as CRISPR.
            first = re.sub(r"CRISPR-G[0-9]", "CRISPR", subs[0])
            first = re.sub(r"GII-I[I]*", "GII", first)
            firsthit = subs
        elif (ln < 2):
            ln += 1
            second = re.sub(r"CRISPR-G[0-9]", "CRISPR",
                            subs[0])  # CRISPR has 6 subclasses, we don't need the subclass info.
            second = re.sub(r"GII-I[I]*", "GII", second)  # GII has 2 subclasses, we don't need the subclass info.
            if subs[0] == firsthit[0]:  # This means the RT has two RVT_1 domains for example two RVT-CRISPR-G2.
                firsthit[0] = first
                firsthit.append(firsthit[6])
                OUT.append(firsthit)
                secondhit = subs
                secondhit[0] = second
                secondhit.append(secondhit[6])
                OUT.append(secondhit)
                flagS = 1  # Done
            elif first == second:  # This means that first and second hits are from the same class (both are CRISPR, or both are GII)
                firsthit.append(float(float(firsthit[6]) / float(subs[6])))
                firsthit[0] = first
                OUT.append(firsthit)
                flagS = 1  # Done
            else:
                secondhit = subs
        elif (ln < 3) and flagS != 1 and float(float(firsthit[6]) / float(secondhit[6])) > eps:
            ln += 1
            third = re.sub(r"CRISPR-G[0-9]", "CRISPR", subs[0])
            third = re.sub(r"GII-I[I]*", "GII", third)

            if (first == third) or (second == third):
                firsthit[0] = first + "/" + second  # firsthit & secondhit
                firsthit.append(float(float(firsthit[6]) / float(
                    secondhit[6])))  # we used e-value2 to be more strict, because e-value2 < e-value3
            else:
                firsthit[0] = first + "/" + second + "/" + third  # firsthit & secondhit & thirdhit
                firsthit.append(float(float(firsthit[6]) / float(subs[6])))
            query.append(firsthit[3])
            OUT.append(firsthit)
            flagS = 1
        else:
            continue
    else:
        continue
if flagS == 0:
    if ln == 1:
        firsthit[0] = first
        firsthit.append(firsthit[6])  # This means that there was only one tophit for previous RT
    elif ln == 2:
        if float(float(firsthit[6]) / float(secondhit[6])) > eps:
            if (first == second):
                firsthit[0] = first
            else:
                firsthit[0] = first + "/" + second
                query.append(firsthit[3])
        else:
            firsthit[0] = first
        firsthit.append(float(float(firsthit[6]) / float(secondhit[6])))
    OUT.append(firsthit)
    flagS = 1
if len(OUT) > 0:
    out = open(file_paths_dictionary["RVT_domtblout"], "w")
    for idx in range(len(OUT)):
        print >> out, ' '.join(map(str, OUT[idx]))

if (len(query) > 0):
    query_faa = open(file_paths_dictionary["QUERY_faa"], "w")
    for idx in range(len(query)):
        inf = open(file_paths_dictionary["faa_file"], "r")
        for aline in inf:
            if aline[0] == '#':
                continue
            elif aline[0] == '>' and query[idx] in aline:
                print >> query_faa, aline, inf.next(),
                break
        inf.close()
    query_faa.close()
if len(verified) > 0:
    genomicneighborhood_file = open(file_paths_dictionary["Genomic_neighborhood_file"], "w")
    genomicneighborhood_faa = open(file_paths_dictionary["Genomic_neighborhood_faa"], "w")
    genomicneighborhood_seen = list()
    for reversetranscriptase in verified:
        genomicneighborhood = list()
        inf = open(file_paths_dictionary["gff_file"], "r")
        features = []
        seqid = re.sub(r"_[0-9]*_[0-9]*_[-|+]", "", reversetranscriptase)
        for aline in inf:
            if aline[0] == '#':
                continue
            else:
                subs = aline.strip().split("\t")
                if len(subs) != 9:
                    print aline
                    sys.exit("wrong gff input")
                if subs[0] != seqid:
                    continue
                if subs[2] == 'region':
                    seqlen = int(subs[4]) - int(subs[3]) + 1
                    continue
                notes = subs[-1].split(";")
                ID = ""
                for anote in notes:
                    if "ID=" in anote:
                        ID = anote[3:]
                features.append([subs[0], int(subs[3]), int(subs[4]), subs[6], ID])
        inf.close()
        for idx in range(len(features)):
            if features[idx][4] == reversetranscriptase:
                if (idx > 0 and abs(features[idx][1] - features[idx - 1][2]) < intergenic_distance):
                    if (idx > 1 and abs(features[idx - 1][1] - features[idx - 2][2]) < intergenic_distance):
                        genomicneighborhood.append(str(features[idx - 2][4]))
                        genomicneighborhood.append(str(features[idx - 1][4]))
                    else:
                        genomicneighborhood.append(str(features[idx - 1][4]))
                genomicneighborhood.append(features[idx][4])
                if (idx < len(features) - 1 and abs(features[idx + 1][1] - features[idx][2]) < intergenic_distance):
                    if (idx < len(features) - 2 and abs(
                            features[idx + 2][1] - features[idx + 1][2]) < intergenic_distance):
                        genomicneighborhood.append(str(features[idx + 1][4]))
                        genomicneighborhood.append(str(features[idx + 2][4]))
                    else:
                        genomicneighborhood.append(str(features[idx + 1][4]))
                break
        print >> genomicneighborhood_file, ','.join(map(str, genomicneighborhood))
        for gene in genomicneighborhood:
            if gene not in genomicneighborhood_seen:
                genomicneighborhood_seen.append(gene)
                inf = open(file_paths_dictionary["faa_file"], "r")
                for aline in inf:
                    if aline[0] == '#':
                        continue
                    elif aline[0] == '>' and gene in aline:
                        print >> genomicneighborhood_faa, aline, inf.next(),
                        break
                inf.close()
    genomicneighborhood_file.close()
    genomicneighborhood_faa.close()
