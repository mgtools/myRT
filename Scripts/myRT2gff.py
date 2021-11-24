#!/usr/bin/env python3
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
# Last update Nov 5 2020 by Fatemeh Sharifi

import sys
import re
DIR=sys.argv[2]
ACC=sys.argv[1]
domfile=DIR+"/"+ACC+"-gn-dom.txt"
listfile=DIR+"/"+ACC+"-gn.list"
RTfile=DIR+"/"+ACC+"-RT.gff"
GFFfile=DIR+"/"+ACC+"-gn.gff"
GN_dom=DIR+"/"+ACC+"-gn-dom.list"
Lenfile=DIR+"/"+ACC+".fna.len"
fgsgff=DIR+"/"+ACC+"-FGS.gff"

gn_dom=open(GN_dom,"w")
inf = open(listfile, "r")
dom={}
casgn={}
RTs=[]
gn_list=[]
matches=["cas","cmr","csx","Cmr","Cas","Csx"]
for aline in inf:
        casdom=0
        if aline[0] == '#':
                continue
        subs = aline.strip().split(",")
        domlist=[]
        RT_line=[]
        rt_curr=""
        for i in range(len(subs)):
                domf=open(domfile,"r")
                for domline in domf:
                        par = domline.strip().split()
                        if subs[i] == par[0] and par[1]!="-":
                                domlist.append(par[1])
                                if par[0] not in gn_list:
                                        gn_list.append(par[0])
                                if any(x in par[1] for x in matches) and par[1]!="SAM_DCase_Bsu":
                                        casdom=casdom+1
                                if "RVT-" in par[1]:
                                        #RT_line captures all RTs in one line of domf
                                        RT_line.append(par[0])
                                        #RT captures all RTs
                                        RTs.append(par[0])
                                        rt_curr=par[0]
        for j in RT_line:
                idx = RTs.index(j)
                casgn[idx]=casdom
        domslist=','.join(map(str, domlist))
        thisline=rt_curr + " " + domslist 
        print (thisline,file=gn_dom)
        domf.close()     
inf.close()
gn_dom.close()
#print RTs
des={}
inf = open(fgsgff, "r")
for aline in inf:
        subs = aline.strip().split("\t")
        if aline[0] == '#':
                continue
        elif subs[2] == "region":
                continue
        else:
                subs2 = subs[-1].split(";")
                flagID=0
                for ades in subs2:
                        if (ades[:3] == 'ID=') and ades[3:] in gn_list:
                                flagID=1
                                ourgn=ades[3:]
                                idx=gn_list.index(ourgn)
                                descr=""
                        elif flagID==1 and (ades[:5] == 'gene=') and ades[5:]!="unk" and ades[5:]!="unknown":
                                        descr=descr+ades+";"
                        elif flagID==1 and (ades[:4] == "Ori=") and ades[4:]!="unk" and ades[4:]!="unknown":
                                        descr=descr+ades+";"
                        if flagID==1:
                                        des[idx]=descr
inf.close()
seqseen=[]
geneseen=[]
gfffile=open(GFFfile,"w")
rtfile=open(RTfile,"w")
inf = open(domfile, "r")
for aline in inf:
        if aline[0] == '#':
                continue
        subs = aline.strip().split()
        seqid=re.sub(r"_[0-9]*_[0-9]*_[-|+]","",subs[0])
        par=subs[0][len(seqid)+1:].split("_")
        startgene=int(par[0])
        stopgene=int(par[1])
        strand=par[2]
        if seqid not in seqseen:
                seqseen.append(seqid)
                lenfile=open(Lenfile,"r")
                for l in lenfile:
                        col = l.strip().split()
                        if seqid==col[1]:
                                lenseq=col[0]
                                output=seqid+"\tGenbank\tregion\t1\t"+lenseq+"\t.\t"+strand+"\t.\tID="+seqid
                                print (output,file=gfffile)
                lenfile.close()
        if subs[0] not in geneseen:
                geneseen.append(subs[0])
                if subs[0] in gn_list:
                        print (gn_list)
                        print (des)
                        print (subs[0])
                        idx=gn_list.index(subs[0])
                        if idx in des:
                                descr=";" + des[idx]
                        else:
                                descr=";"
                else:
                        descr=""
                output= seqid + "\tFGS\tCDS\t" + str(startgene) + "\t" + str(stopgene) + "\t.\t"+ strand+ "\t0\tID="+subs[0] + ";what=gene"+descr
                print (output,file=gfffile)
        domstart=startgene+(int(subs[2])*3)
        domend=startgene+(int(subs[3])*3)
        if "RVT-" in subs[1]:
                if "RVT-CRISPR"==subs[1]:
                        if casgn[RTs.index(subs[0])]==0:
                                subs[1]="RVT-CRISPR-like"
                output=seqid + "\tmyRT\tdom\t" + str(domstart)+ "\t" + str(domend) + "\t.\t" +strand+ "\t0\tID=" + seqid + "_" + str(domstart) + "_" + str(domend) + "_"+ strand + ";what=RT;des="+subs[1]
                print (output,file=gfffile)
                out=DIR+ "/" + ACC + "\t" + seqid +"\tmyRT\t" + seqid + "_" + str(startgene) + "_" + str(stopgene) + "_" +strand+"\t" +subs[1]
                print (out,file=rtfile)
        elif subs[1]!="-":
                output=seqid + "\tmyRT\tdom\t" + str(domstart) + "\t" + str(domend) + "\t.\t" +strand+ "\t0\tID=" + seqid + "_" + str(domstart) + "_" + str(domend) + "_"+ strand + ";what=dom;des="+subs[1]
                print (output,file=gfffile)
rtfile.close()
gfffile.close()
inf.close()

