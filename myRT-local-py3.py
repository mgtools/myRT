#!/usr/bin/env python3
# Last update Nov 2021 by Fatemeh Sharifi
# myRT-local-py3.py Query BA000019.2 QUERY/BA000019.2/BA000019.2.fna
# myRT-local-py3.py Query BA000019.2 QUERY/BA000019.2/BA000019.2.fna QUERY/BA000019.2/BA000019.2.gff

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
import sys
import os
import subprocess
python="python3"


if "check_output" not in dir( subprocess ): #old python doesn't support check_output
    def f(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = f

def file_put_contents(filename, content):
    subprocess.check_output("echo \"" + content + "\" > " + filename, shell=True)   

cleanup=True

if len(sys.argv) < 4:
   sys.exit("usage: myRT-local-py3.py output-base output-prefix input-fna <input-gff>\n")


cmd=sys.argv[0]

#Local Directory
#myRT="/home/team/myRT"
myRT=os.getcwd()

# This works with Pfam-A.hmm and cdd.hmm too, but we recommend cdd-PfamA.hmm
# Download from here https://omics.informatics.indiana.edu/myRT/share/download.php 
pfammodel = os.path.join( myRT , "Models/cdd-pfamA/cdd-pfamA.hmm")
# =====================

# Customized hmm model built from RVT_1 motif sequences 
RVTAllmodel =os.path.join(myRT , "Models/HMM/RVT-All.hmm")
# =====================

# For gene prediction (if you use ncbi gff files you can skip using FragGeneScan )
FragGeneScan =os.path.join(myRT , "bin/FragGeneScan1.31/run_FragGeneScan.pl")
hmmscan = os.path.join (myRT , "bin/hmmer-3.2/src/hmmscan")
# =====================

# FINDING SIMIALR RTS
blastp = os.path.join(myRT , "bin/blast+/bin/blastp")
RTs_collection = os.path.join(myRT, "Models/RTs-collection.faa")
# If you make changes to the collection, use makeblastdb -in RTs-collection.faa -dbtype 'prot' to update it.
# =====================

# PHYLOGENETIC ANALAYSES
hmmalign = os.path.join(myRT , "bin/hmmer-3.2/src/hmmalign")
pplacer = os.path.join(myRT, "bin/pplacer-Linux-v1.1.alpha19/pplacer")
# Please check the script to download required R packages (treeio, castor, phangorn); tar files are in bin folder.
JplaceParser = os.path.join(myRT , "Scripts/JplaceParser.R")
# Profiles for the reference tree built using taxtastic 0.8.3 based on pplacer manual
Ref_tree_package = os.path.join (myRT , "Models/myRT-FastTree2.refpkg/" )
Ref_tree_sto = os.path.join (myRT , "Models/myRT-FastTree2.refpkg/RVT-ref.sto")
Ref_tree_hmm = os.path.join (myRT , "Models/myRT-FastTree2.refpkg/RVT-ref.hmm")
# =====================

# myRT SCRIPTS: 
extract_genomic_neighborhood=os.path.join(myRT , "Scripts/extract_gn.py")
extract_domains=os.path.join(myRT , "Scripts/hmm2dom_myRT.py")
build_gff_file=os.path.join(myRT , "Scripts/myRT2gff.py")
# =====================
# To generate protein sequences if gff is provided; otherwise uses fraggenescan)
gff2seq=os.path.join(myRT , "/Scripts/gff2seq.pl")
# To check the fasta file
check_input=os.path.join(myRT , "/Scripts/check_file.pl")
# =====================

# Required folders
if os.path.exists(myRT + "/Models") == False:
   sys.exit("myRT directory not set correctly: Models folder not found")
if os.path.exists(myRT + "/Scripts") == False:
   sys.exit("myRT directory not set correctly: Scripts folder not found")
if os.path.exists(myRT + "/bin") == False:
   sys.exit("myRT directory not set correctly: bin folder not found")
# =====================

# Output Directory 
# Example: QUERY
output_base=sys.argv[1]
# Example: BA000019.2
output_prefix=sys.argv[2]
# Example: QUERY/BA000019.2
output_dir=os.path.join(output_base,output_prefix)
# Example: QUERY/BA000019.2/BA000019.2
output=os.path.join(output_base,output_prefix,output_prefix)
orign_input_name=output_prefix
#Example QUERY/BA000019.2/BA000019.2.fna
old_input_file=sys.argv[3]

old_gff_file="doesnotexist"
old_gff_valid=False

# If gff is added skips FGS 
if len(sys.argv) >= 5:
   old_gff_file=sys.argv[4]
   old_gff_valid=True 

if os.path.exists(output_dir) == False:
    subprocess.check_output("mkdir -p " + output_dir, shell=True)


# RVT_faa is generated based on searching all genes against RVT-All.hmm, and RT_faa is generated based on searching RVT_faa against cdd-pfamA.hmm (See Figure 1 in our manuscipt). 

file_info_dictionary = {
        "new_input_file" :  ".fna",
        "RT_gn_list" : "-gn.list",
        "RT_gn_faa" : "-gn.faa",
        "fgs_faa" : "-FGS.gff",
        "RVT_faa" : "-RVT.faa",
        "RT_faa" : "-RT.faa",
        "RVT_tmp" : "-RVT.tmp",
        "RVT_cdd_pfamA" : "-RVT-cdd-pfamA.domtblout",
        "RVT_domtblout" : "-RVT.domtblout",
        "RT_gn_cdd_pfamA" : "-gn-cdd-pfamA.domtblout",
        "RT_gn_dom" : "-gn-dom.txt",
        "RT_vs_ALL" : "-RVT-vs-ALL.m8",
        "RT_gn_gff" : "-gn.gff",
        "RT_gn_dom_list" : "-gn-dom.list",
        "fgs_out" : "-FGS",
        "fgs_faa" : "-FGS.faa",
        "fgs_faa_tmp" : "-FGS.faa.tmp",
        "fgs_ffn" : "-FGS.ffn",
        "gff_fgs_file" : "-FGS.gff",
        "seq_note" :".note",
        "error_log" : ".log",
        "Done_file" : ".done",
        "RT_query" : "-query.faa",
        "RT_combo" : "-combo.sto",
        "RT_jplace": "-combo.jplace",
        "RT_tree" : "-combo.tre",
        "RT_pred" : "-pred.txt",
        "len_file" : ".fna.len",
        "R_log" : "-R.log",
        "fna_log": "-fna.log",
        "sm_output" : "-sm.fna"
}

file_paths_dictionary = dict()
for file_type in file_info_dictionary:
    file_paths_dictionary[file_type] = os.path.join(output_dir,output_prefix + file_info_dictionary[file_type])

cmd="rm -f " + file_paths_dictionary["error_log"]     #important for successful rerun
subprocess.check_output(cmd, shell=True)

print (output)


if os.path.exists(file_paths_dictionary["Done_file"]):
    log = open(file_paths_dictionary["error_log"], "w")
    print ("Results found in the folder, if you need to redo, delete the results then rerun this script",file=log)
    sys.exit("Results found in the folder " + output_dir + "--if you need to redo the annotation, delete the results then rerun this script\n")

#check input file
'''
check_result = subprocess.check_output("perl " + check_input +  " + old_input_file + " " + file_paths_dictionary["new_input_file"] + " " + file_paths_dictionary["seq_note"] + ">" + file_paths_dictionary["fna_log"] , shell=True)

if os.path.exists(file_paths_dictionary["fna_log"]) and os.stat(file_paths_dictionary["fna_log"]).st_size > 0:
    cmd="mv " + file_paths_dictionary["fna_log"] + " " + file_paths_dictionary["error_log "]
    subprocess.check_output(cmd, shell=True)
    sys.exit("Input file is not in FASTA format")
       
'''
if os.path.exists(old_input_file) and not os.path.exists(file_paths_dictionary["new_input_file"]):
        cmd="mv " + old_input_file + " " + file_paths_dictionary["new_input_file"]
        subprocess.check_output(cmd, shell=True)

#prepare protein file when gff file is given
if os.path.exists(old_gff_file) and os.stat(old_gff_file).st_size > 0 and ( not os.path.exists(file_paths_dictionary["fgs_faa"]) or os.stat(file_paths_dictionary["fgs_faa"]).st_size == 0 ):
        #if filecmp.cmp(old_gff_file,gff_fgs_file)==False:
        cmd="cp" + " " +  old_gff_file + " " + file_paths_dictionary["gff_fgs_file"]
        subprocess.check_output(cmd, shell=True)
        print ("Preparing protein sequences according to the gff file...\n")
        result=subprocess.check_output("perl " + myRT + "/Scripts/gff2seq.pl " + output, shell=True)
        if os.path.getsize(file_paths_dictionary["fgs_faa"]) <= 10:
                old_gff_valid=False
#predict proteins using FragGeneScan if gff was not provided, or the given gff did not contain information about proteins
if not old_gff_valid and ( not os.path.exists(file_paths_dictionary["fgs_faa"]) or os.stat(file_paths_dictionary["fgs_faa"]).st_size == 0 ):

        print ("Running FragGeneScan...\n")
        fraggenescan=FragGeneScan +" -genome=" + file_paths_dictionary["new_input_file"] + " -out=" + file_paths_dictionary["fgs_out"] + " -complete=1 -train=complete thread=9"
        try:
                subprocess.check_output(fraggenescan, shell=True)
        except subprocess.CalledProcessError as error:
                print (error)

if not os.path.exists(file_paths_dictionary["fgs_faa"]):
        log = open(file_paths_dictionary["error_log"], "w")
        print ("Generation of protein sequences failed!",file=log)
        sys.exit("Generation of protein sequences failed!Please make sure your sequence file is not empty")

# Searches all the genes against RVT-All.hmm (customized hmm model)
if not os.path.exists(file_paths_dictionary["RVT_tmp"]) or os.stat(file_paths_dictionary["RVT_tmp"]).st_size == 0:
        print ("Searching the genes against RVT-All.hmm (e-value < 0.0001) ...\n")
        cmd=hmmscan + " --noali -E 0.0001 --domE 0.0001   --domtblout " +  file_paths_dictionary["RVT_tmp"] + " " + RVTAllmodel + " " + file_paths_dictionary["fgs_faa"]
        try:
                subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as error:
                print (error)

if not os.path.exists(file_paths_dictionary["RVT_tmp"]) or os.stat(file_paths_dictionary["RVT_tmp"]).st_size == 0:
        log = open(file_paths_dictionary["error_log"], "w")
        print ("No RT domain was found!",file=log)
        sys.exit("No RT domain was found!")

# Extracts putative RTs
if not os.path.exists(file_paths_dictionary["RVT_faa"]) or os.stat(file_paths_dictionary["RVT_faa"]).st_size == 0 :
        print ("Extracting putative RTs...\n")
        cmd="cat " + file_paths_dictionary["RVT_tmp"] + """| awk '($1!~"#") {print $4}' | uniq | while read line; do grep -A 1 ${line} """ + file_paths_dictionary["fgs_faa"] + " >> " +  file_paths_dictionary["RVT_faa"] + ";done"
        try:
                subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as error:
            print (error)
if not os.path.exists(file_paths_dictionary["RVT_faa"]) or os.stat(file_paths_dictionary["RVT_faa"]).st_size == 0:
        log = open(file_paths_dictionary["error_log"], "w")
        print ("No RT domain was found!",file=log)
        sys.exit("No RT domain was found!")
# Searches putative RTs against PfamA.hmm (to remove False Positives if any)      
if not os.path.exists(file_paths_dictionary["RVT_cdd_pfamA"]) or os.stat(file_paths_dictionary["RVT_cdd_pfamA"]).st_size == 0:
        print ("Searching the putative RTs against CDD domains (e-value < 0.01)...\n")
        cmd=hmmscan + " --noali -E 0.01 --domE 0.01 --domtblout " +  file_paths_dictionary["RVT_cdd_pfamA"] + " " + pfammodel + " " + file_paths_dictionary["RVT_faa"]
        try:
                subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as error:
                print (error)
# Verifying RTs, and extracting genomic neighborhood of verified RTs
if not os.path.exists(file_paths_dictionary["RT_gn_faa"]) or os.stat(file_paths_dictionary["RT_gn_faa"]).st_size == 0 :
        print ("Verifying RTs, and extracting genomic neighborhood of verified RTs...\n")
        cmd=python + " " + extract_genomic_neighborhood + " " + output_prefix + " " + output_dir
        try:
                subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as error:
                print (error)

if not os.path.exists(file_paths_dictionary["RVT_domtblout"]) or os.stat(file_paths_dictionary["RVT_domtblout"]).st_size == 0:
        log = open(file_paths_dictionary["error_log"], "w")
        print ("No RT domain was found!",file=log)
        sys.exit("No RT domain was found!")

# seqlen
if not os.path.exists(file_paths_dictionary["len_file"]) or os.stat(file_paths_dictionary["len_file"]).st_size == 0:
        cmd = """awk '$0 ~ ">" {seqid=substr($1,2); seqlen[seqid]=0;} $0 !~ ">" {seqlen[seqid]=seqlen[seqid]+length($0);} END { for (seq in seqlen) print seqlen[seq] " " seq; }' """ + file_paths_dictionary["new_input_file"] + " > " + file_paths_dictionary["len_file"]
        subprocess.check_output(cmd, shell=True)

# blastp to search for similar RTs in our collection
if not os.path.exists(file_paths_dictionary["RT_vs_ALL"]) or os.stat(file_paths_dictionary["RT_vs_ALL"]).st_size == 0 :
        print ("Running Blastp...\n")
        cmd="nohup " + blastp + "  -query " + file_paths_dictionary["RT_faa"] + " -db "  + RTs_collection + " -num_alignments 1 -evalue 1e-3 -outfmt 6 -out " +  file_paths_dictionary["RT_vs_ALL"]
        try:
               subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as error:
                print (error)

# pplacer to incorprate the phylogenetic information to improve the accuary (this step is only done for low confidence/unclassified predictions)
# Aligning the queries and the reference 
if not os.path.exists(file_paths_dictionary["RT_combo"]) and os.path.exists(file_paths_dictionary["RT_query"]):
       # If there are Us in the sequence that cause pplacer to fail uncomment the following block to  remove them
       # sed function
       '''
       cmd="sed -i '/^[A-Z]/s/U//g'" + " " + file_paths_dictionary["RT_query"]
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print error
       '''
       print ("Incorporating phylogenetic information ..\n")
       cmd=hmmalign + " -o " + file_paths_dictionary["RT_combo"] + " --mapali " + Ref_tree_sto + " " + Ref_tree_hmm + " " + file_paths_dictionary["RT_query"] 
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
              print (error)
              print ("If the error is because of the Us in the sequence remove them using the sed function provided in the myRT script")
       
if not os.path.exists(file_paths_dictionary["RT_jplace"]) and os.path.exists(file_paths_dictionary["RT_combo"]):
       cmd=pplacer + " -c " + Ref_tree_package + " " + file_paths_dictionary["RT_combo"] + " -o " + file_paths_dictionary["RT_jplace"] 
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
            # Usually this error is because of U in the sequence, use 
            print (error)

# If you want to have a tree with the placed queries uncomment this block
'''
if not os.path.exists(RT_tree) and os.path.exists(RT_jplace):
       cmd=guppy + " tog " + file_paths_dictionary["RT_jplace"] + " -o " file_paths_dictionary["RT_tree"]
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print (error) 
'''
# Parsing pplacer output (.jplace) (prediction based on placement in the phylogenetic tree)
if (not os.path.exists(file_paths_dictionary["RT_pred"]) or os.stat(file_paths_dictionary["RT_pred"]).st_size == 0 ) and os.path.exists(file_paths_dictionary["RT_jplace"]):
       cmd="Rscript " + JplaceParser + " " + file_paths_dictionary["RT_jplace"] + " > " + file_paths_dictionary["RT_pred"] + " 2>>" + file_paths_dictionary["R_log"]
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print (error)
# Searching for other domains in the putative RT, and domains in the flanking genes
print ("Annotating the genomic neighborhood of RTs ...\n")
if os.path.exists(file_paths_dictionary["RT_gn_faa"]) and os.stat(file_paths_dictionary["RT_gn_faa"]).st_size > 0 and (not os.path.exists(file_paths_dictionary["RT_gn_cdd_pfamA"]) or os.stat(file_paths_dictionary["RT_gn_cdd_pfamA"]).st_size == 0):
       cmd=hmmscan + " --noali -E 0.001 --domE 0.001 --domtblout " +  file_paths_dictionary["RT_gn_cdd_pfamA"] + " " + pfammodel + " " + file_paths_dictionary["RT_gn_faa"]
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print (error)
# Combining hmmscan results, and pplacer (if used)
if not os.path.exists(file_paths_dictionary["RT_gn_dom"]): 
       print ("Preparing the domain files ...\n")     
       if os.path.exists(file_paths_dictionary["RT_pred"]): 
        cmd=python + " " + extract_domains + " " + file_paths_dictionary["RT_gn_list"] + " " +  file_paths_dictionary["RT_gn_cdd_pfamA"] + " " + file_paths_dictionary["RVT_domtblout"] + " " + file_paths_dictionary["RT_gn_dom"] + " " + file_paths_dictionary["RT_pred"]
       else:
        cmd=python + " " + extract_domains + " " + file_paths_dictionary["RT_gn_list"] + " " +  file_paths_dictionary["RT_gn_cdd_pfamA"] + " " + file_paths_dictionary["RVT_domtblout"] + " " + file_paths_dictionary["RT_gn_dom"]
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print (error)
# Preparing a gff file (RTs and flanking genes and domains) for visualization, and one containing the RTs and their predicted class  
# Also makes a list of domains in the genomic neighborhood of each RT
if not os.path.exists(file_paths_dictionary["RT_gn_gff"]):
       print ("Preparing the gff files ...\n")
       cmd=python + " " + build_gff_file + " " + output_prefix + " " + output_dir
       try:
                subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print (error)

# For METAGENOMES uncomment this block
'''
cmd="""awk -F $' ' 'FNR==NR{ split($9,array,";"); what=substr(array[2],6); if(what=="RT") idlist[$1]; next} {s=substr($1,1,1); if (s==">") { id=substr($1,2); if (id in idlist) { valid = 1; print ">" id; } else { valid = 0; }} else if (valid) print $0; }' """ + file_paths_dictionary["RT_gn_gff"] + " " + file_paths_dictionary["new_input_file"] + " > " + file_paths_dictionary["sm_output"]
subprocess.check_output(cmd, shell=True)
'''

if os.path.exists(file_paths_dictionary["error_log"]):
    log.close()
cmd="touch " + file_paths_dictionary["Done_file"]
subprocess.check_output(cmd, shell=True)
print ("Done!\n")

