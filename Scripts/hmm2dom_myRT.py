#!/usr/bin/env python

# Last update Nov 5 2020 by Fatemeh Sharifi

import sys

#myRT/Scripts/hmm2dom_myRT.py NC_002950.2-gn.list  NC_002950.2-gn-cdd-PfamA.domtblout NC_002950.2-RVT.domtblout  NC_002950.2-gn-dom.txt

if len(sys.argv) < 4:
        sys.exit(sys.argv[0] + " gn.list vs-cdd-PfamA vs-RVT out-file")

# adding all the genes to idlist
idlist=[]
inf = open(sys.argv[3], "r")
for aline in inf:
        if aline[0] == '#':
                continue
        subs = aline.strip().split()
        if subs[3] not in idlist and subs[11] > 0.001:
            idlist.append(subs[3])
inf.close()
inf = open(sys.argv[2], "r")
for aline in inf:
        if aline[0] == '#':
            continue
        subs = aline.strip().split()
        if subs[3] not in idlist:
            idlist.append(subs[3])
inf.close()
print idlist


# To keep the sortings of the genes
gn = open(sys.argv[1], "r")
seqid, seq2dom = [], []
for aline in gn:
        line=aline.strip().split(",")
        for i in range(len(line)):
        	if line[i] not in seqid:
            		seqid.append(line[i])
            		tmp=[]
            		seq2dom.append(tmp)
        	if line[i] not in idlist:
            		idx=seqid.index(line[i])
            		seq2dom[idx].append(['-', 0, 0, 0])    
                
gn.close()
print seqid

correctLabel={}
correctevalue={}
score={}
if len(sys.argv) > 5:
	Correction=open(sys.argv[5], "r")
	for aline in Correction:
		subs = aline.strip().split()
                if subs[1] =="RVT-CL":
                	subs[1]="RVT-CRISPR-like"
                if subs[1] == "RVT-CLB":
                        subs[1]="RVT-CRISPR-like"
                if subs[1] == "RVT-CLC":
                        subs[1]="RVT-CRISPR-like"
                if subs[1] == "RVT-Retron":
                        subs[1]="RVT-Retrons"
        	correctLabel[subs[0]]=subs[1]
        	correctevalue[subs[0]]=1/(2**(float(subs[2])*100)) #score to e-value
      		score[subs[0]]=float(subs[2])
	Correction.close()

# assign RVT
RVT= open(sys.argv[3], "r")
added=[]
for aline in RVT:
        if aline[0] == '#':
            continue
        flag=0
        subs = aline.strip().split()
        if subs[3] in added:
            continue
        if subs[3] not in seqid: 
            seqid.append(subs[3])
            tmp=[]
            seq2dom.append(tmp)
        x=range(int(subs[17]), int(subs[18]))
        idx = seqid.index(subs[3])
        dom=seq2dom[idx]
        for d in dom:
            y=range(d[2],d[3])
            l=list(set(x) & set(y))
            if l:                   #Check for overlaps
                 flag=1
                 break
        if flag==0:
                label=subs[0].split("/")
                if (subs[3] in correctLabel.keys()):
			subs[11]=correctevalue[subs[3]]
                        correctlab=correctLabel[subs[3]].replace('-like','')
                	if correctLabel[subs[3]] in subs[0] or correctlab in subs[0] :
                		subs[0]=correctLabel[subs[3]]
                	elif (float(score[subs[3]]) > 0.60 ):
                        	subs[0]=correctLabel[subs[3]]+"/"+subs[0]
                	else:
                        	subs[0]=subs[0]+"/"+correctLabel[subs[3]]
        	seq2dom[idx].append([subs[0], subs[11], int(subs[17]), int(subs[18])])

RVT.close()

#assign pfam
inf = open(sys.argv[2], "r")
for aline in inf:
        if aline[0] == '#':
                continue
        flag=0
        subs = aline.strip().split()
        if subs[3] not in seqid:
            seqid.append(subs[3])
            tmp=[]
            seq2dom.append(tmp)
        idx = seqid.index(subs[3])
        dom=seq2dom[idx]
        x=range(int(subs[17]), int(subs[18]))
        for d in dom:
                y=range(d[2],d[3])
                l=list(set(x) & set(y))
                if l:                   #Check for overlaps
                        flag=1
                        break
        if(flag==0):
            seq2dom[idx].append([subs[0], subs[11], int(subs[17]), int(subs[18])])
inf.close()
print seq2dom
out = open(sys.argv[4], "w")
for idx in range(len(seqid)):
	dom = seq2dom[idx]
        domsorted = sorted(dom, key=lambda tmp: float(tmp[2]))  #sort based on start
	for d in domsorted:
		print >>out, seqid[idx],d[0], d[2], d[3], d[1]
out.close()        
