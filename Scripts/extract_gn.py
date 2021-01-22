#!/usr/bin/env python

# Last update Dec 14 2020 by Fatemeh Sharifi

import sys
import re
DIR=sys.argv[2]
print(DIR)
ACC=sys.argv[1]
RVT_tmp=DIR+"/"+ACC+"-RVT.tmp"
RVT_cdd_pfam=DIR+"/"+ACC+"-RVT-cdd-pfamA.domtblout"
RVT_domtblout=DIR+"/"+ACC+"-RVT.domtblout"
RVT_faa=DIR+"/"+ACC+"-RVT.faa"
RT_faa=DIR+"/"+ACC+"-RT.faa"
QUERY_faa=DIR+"/"+ACC+"-query.faa"
GN_file=DIR+"/"+ACC+"-gn.list"
GN_faa=DIR+"/"+ACC+"-gn.faa"
gfffile=DIR+"/"+ACC+"-FGS.gff"
faafile=DIR+"/"+ACC+"-FGS.faa"
RTlist=[]
verified=[]
#Verifies the RTs based on cdd/pfam hits
inf = open(RVT_cdd_pfam, "r")
matches=["RVT_1","RVT_N","Cas1","Abi","RT_","cas1","GIIM","Cas6","cas6","Cas2","cas2"]
for aline in inf:
        if aline[0] == '#':
                continue
        subs = aline.strip().split()
        if subs[3] not in RTlist:
            RTlist.append(subs[3])
            ln=1
            flagV=0
            if "RVT_1" in subs[0] or "RT_" in subs[0] or "RVT_N" in subs[0]:
            	verified.append(subs[3])
                flagV=2
            elif any(x in subs[0] for x in matches):
                flagV=flagV+1
        elif flagV <2 and ln <=10:
            ln=ln+1
            if flagV <2 and (ln <=10) and any(x in subs[0] for x in matches):
            	flagV=flagV+1
            	if flagV >=2:
                	verified.append(subs[3])
        else:
            continue
inf.close()
#print (verified)
if len(verified)> 0:
	rt_faa=open(RT_faa,"w")
	for v in verified: 
		inf = open(RVT_faa,"r")
        	for aline in inf:
        		if aline[0] == '>' and v in aline:
                        	print >>rt_faa, aline,inf.next(),
                               	break
		inf.close()
        rt_faa.close()
RT=[]
OUT=[]
query=[]
flagS=2
inf = open(RVT_tmp, "r")
for aline in inf:
        if aline[0] == '#':
                continue
        subs = aline.strip().split()
        subs[22]=re.sub(r"CRISPR-G[0-9]","CRISPR",subs[22])
        subs[22]=re.sub(r"GII-I[I]*","GII",subs[22]) 
        if subs[3] in verified:
        	if subs[3] not in RT:
                        if flagS==0:
                        	if ln==1:
                                        firsthit[0]=first
                        		firsthit.append(firsthit[6]) # This means that there was only one tophit for previous RT
                                elif ln==2 and float(float(firsthit[6])/float(secondhit[6])) > 0.00001 :
                                        if (first==second):
                                                firsthit[0]=first
                                        else:
                                                firsthit[0]=first+"/"+second #first & second
                                        	query.append(firsthit[3])
                                        firsthit.append(float(float(firsthit[6])/float(secondhit[6])))
                                else:
                                	firsthit[0]=first
                                	firsthit.append(float(float(firsthit[6])/float(secondhit[6])))
                                OUT.append(firsthit)
                                flagS=1
            		RT.append(subs[3])
                        flagS=0
            		ln=1
                        first=re.sub(r"CRISPR-G[0-9]","CRISPR",subs[0])
                        first=re.sub(r"GII-I[I]*","GII",first)
                        firsthit=subs
        	elif (ln < 2):
                        ln=ln+1
                        second=re.sub(r"CRISPR-G[0-9]","CRISPR",subs[0])           # CRISPR has 6 subclasses, we don't need the subclass info 
                        second=re.sub(r"GII-I[I]*","GII",second)                    # GII has 2 subclasses, we don't need the subclass info
                        if subs[0]==firsthit[0]:        # This means the RT has two RVT_1 domains for example two RVT-CRISPR-G2 
                                firsthit[0]=first
                                firsthit.append(firsthit(6))
                        	OUT.append(firsthit)
                                secondhit=subs
                                secondhit[0]=second
                                secondhit.append(secondhit(6))
                                OUT.append(secondhit)
                                flagS=1  #Solved
                        elif first==second:                #This means that first and second hits are from the same class (both are CRISPR, or both are GII)
				firsthit.append(float(float(firsthit[6])/float(subs[6])))
                                firsthit[0]=first
                                OUT.append(firsthit)
                                flagS=1
                        else:
                        	secondhit=subs
                elif (ln < 3) and flagS!=1 and float(float(firsthit[6])/float(secondhit[6])) > 0.00001:
                        ln=ln+1
                        third=re.sub(r"CRISPR-G[0-9]","CRISPR",subs[0])
                        third=re.sub(r"GII-I[I]*","GII",third)
                        
                        if (second==third):                       
                                firsthit[0]=first+"/"+second             #first & second
                                firsthit.append(float(float(firsthit[6])/float(secondhit[6])))   # we used e-value2 to be more strict, because e-value2 < e-value3
                        elif (first==third):
                                firsthit[0]=first+"/"+second             #first & second
                                firsthit.append(float(float(firsthit[6])/float(secondhit[6])))   # we used e-value2 to be more strict, because e-value2 < e-value3
                        else:
                                firsthit[0]=first+"/"+second+"/"+third # first & second & third 
                                firsthit.append(float(float(firsthit[6])/float(subs[6])))
                        query.append(firsthit[3])
                        OUT.append(firsthit)
                        flagS=1
                else:   	
                	continue
	else:
        	continue
if flagS==0:
                                if ln==1:
                                        firsthit[0]=first
                                        firsthit.append(firsthit[6]) # This means that there was only one tophit for previous RT
                                elif ln==2:
                                        if float(float(firsthit[6])/float(secondhit[6])) > 0.00001 :
                                        	if (first==second):
                                                	firsthit[0]=first
                                        	else:
                                                	firsthit[0]=first+"/"+second #first & second
                                			query.append(firsthit[3])
                                        else:
                                        	firsthit[0]=first
                                	firsthit.append(float(float(firsthit[6])/float(secondhit[6])))
                                OUT.append(firsthit)
                                flagS=1
if len(OUT)> 0:
	out = open(RVT_domtblout, "w")
	for idx in range(len(OUT)):
        	#print (' '.join(map(str, OUT[idx])))
		print >>out,' '.join(map(str, OUT[idx]))

if (len(query) > 0):
	query_faa= open(QUERY_faa,"w")
	for idx in range(len(query)):
        	inf = open(faafile, "r")
        	for aline in inf:
        		if aline[0] == '#':
                		continue
                	elif aline[0] == '>' and query[idx] in aline:
                        	print >>query_faa, aline,inf.next(),
				break
		inf.close()
	query_faa.close()
dist=2000
if len(verified)> 0:
	gn_file = open(GN_file, "w")
	gn_faa = open(GN_faa, "w")
	gn_seen=[]
	for rt in verified:
        	gn=[]
		inf = open(gfffile, "r")
		#print "now process", gfffile, "id", seqid
        	features = []
        	seqid=re.sub(r"_[0-9]*_[0-9]*_[-|+]","",rt)
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
                        	features.append([subs[0], int(subs[3]), int(subs[4]), subs[6],ID])                       
        	inf.close()
		for idx in range(len(features)):
                        if features[idx][4]==rt:
                        	if ( idx > 0 and abs(features[idx][1]-features[idx-1][2]) < dist):
                                	if ( idx > 1 and abs(features[idx-1][1]-features[idx-2][2]) < dist):
 						gn.append(str(features[idx-2][4]))
                                                gn.append(str(features[idx-1][4]))
                                        else:
						gn.append(str(features[idx-1][4]))
                                gn.append(features[idx][4])
                                if (idx < len(features)-1 and abs(features[idx+1][1]-features[idx][2]) < dist):
                                        if ( idx < len(features)-2  and abs(features[idx+2][1]-features[idx+1][2]) < dist):
                                                gn.append(str(features[idx+1][4]))
                                                gn.append(str(features[idx+2][4]))
                                        else:
                                                gn.append(str(features[idx+1][4]))
				break
      		print >>gn_file,','.join(map(str, gn))
		for g in gn:
                	if g not in gn_seen:
                		gn_seen.append(g)	
				inf = open(faafile, "r")
        			for aline in inf:
        				if aline[0] == '#':
						continue
					elif aline[0] == '>' and g in aline:
                        			print >>gn_faa, aline,inf.next(),
						break	
           			inf.close()
	gn_file.close()
	gn_faa.close()       
