#!/usr/bin/env python
import sys
import os
import subprocess
import time
from time import gmtime, strftime
import traceback
import filecmp
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
   sys.exit("usage: mydgr-local.py output-base output-prefix input-fna <input-gff>\n")

cmd=sys.argv[0]
#Local Directory
#myRT="/home/team/myRT"
myRT=os.getcwd()

# This works with Pfam-A.hmm and cdd.hmm too, but we recommend cdd-PfamA.hmm
pfammodel= myRT + "/Models/cdd-pfamA/cdd-pfamA.hmm"

FragGeneScan=myRT + "/bin/FragGeneScan1.31/run_FragGeneScan.pl"
hmmalign = myRT + "/bin/hmmer-3.2/src/hmmalign"
hmmscan = myRT + "/bin/hmmer-3.2/src/hmmscan"
blastp= myRT + "/bin/blast+/bin/blastp"
JplaceParser= myRT + "/Scripts/JplaceParser.R" 
pplacer= myRT+ "/bin/pplacer-Linux-v1.1.alpha19/pplacer"


if os.path.exists(myRT + "/Models") == False:
   sys.exit("myRT directory not set correctly: Models folder not found")
if os.path.exists(myRT + "/Scripts") == False:
   sys.exit("myRT directory not set correctly: Scripts folder not found")
if os.path.exists(myRT + "/bin") == False:
   sys.exit("myRT directory not set correctly: bin folder not found")

output_base=sys.argv[1]
output_prefix=sys.argv[2]

output_dir=output_base + "/" + output_prefix
orign_input_name=output_prefix
output = output_base + "/" + output_prefix + "/" + output_prefix
old_input_file=sys.argv[3]
old_gff_file="doesnotexist"
old_gff_valid=False
if len(sys.argv) >= 5:
   old_gff_file=sys.argv[4]
   old_gff_valid=True 

if os.path.exists(output_dir) == False:
	subprocess.check_output("mkdir -p " + output_dir, shell=True)
flagE=0
new_input_file=output + ".fna"
myRT_gn_list=output+"-gn.list"
myRT_gn_faa=output+"-gn.faa"
fgs_faa=output + "-FGS.gff"
myRT_RVT_faa=output+"-RVT.faa"
myRT_RT_faa=output+"-RT.faa"
myRT_RVT_tmp=output+"-RVT.tmp"
myRT_RVT_cdd_pfamA=output+"-RVT-cdd-pfamA.domtblout"
myRT_RVT=output+"-RVT.domtblout"
myRT_gn_cdd_pfamA=output+"-gn-cdd-pfamA.domtblout"
myRT_gn_dom=output+"-gn-dom.txt"
myRT_vs_ALL=output+"-RVT-vs-ALL.m8"
myRT_gn_gff=output+"-gn.gff"
myRT_gn_dom_list=output+"-gn-dom.list"
fgs_out=output+"-FGS"
fgs_faa=output + "-FGS.faa"
fgs_faa_tmp=output + "-FGS.faa.tmp"
fgs_ffn=output+"-FGS.ffn"
gff_fgs_file=output+"-FGS.gff"
seq_note=output+".note"
error_log=output+".log"
Done_file=output+".done"
myRT_query=output+"-query.faa"
myRT_combo=output+"-combo.sto"
myRT_jplace=output+"-combo.jplace"
myRT_tree=output+"-combo.tre"
myRT_pred=output+"-pred.txt"
len_file=output+".fna.len"

cmd="rm -f " + error_log     #important for successful rerun
subprocess.check_output(cmd, shell=True)

print output


if os.path.exists(Done_file):
        log = open(error_log, "w")
        print >> log, "Results found in the folder, if you need to redo, delete the results then rerun this script"
	sys.exit("Results found in the folder " + output_dir + "--if you need to redo the annotation, delete the results then rerun this script\n")

#check input file
check_result = subprocess.check_output("perl " + myRT + "/Scripts/check_file.pl " + old_input_file + " " + new_input_file + " " + seq_note, shell=True)


if os.path.exists(old_input_file) and not os.path.exists(new_input_file):
        cmd="mv " + old_input_file + " " + new_input_file
        subprocess.check_output(cmd, shell=True)

#prepare protein file when gff file is given
if os.path.exists(old_gff_file) and os.stat(old_gff_file).st_size > 0 and ( not os.path.exists(fgs_faa) or os.stat(fgs_faa).st_size == 0 ):
        #if filecmp.cmp(old_gff_file,gff_fgs_file)==False:
        cmd="cp" + " " +  old_gff_file + " " + gff_fgs_file
        subprocess.check_output(cmd, shell=True)
        print "Preparing protein sequences according to the gff file...\n"
        result=subprocess.check_output("perl " + myRT + "/Scripts/gff2seq.pl " + output, shell=True)
        if os.path.getsize(fgs_faa) <= 10:
                old_gff_valid=False
#predict proteins using FragGeneScan if gff was not provided, or the given gff did not contain information about proteins
if not old_gff_valid or ( not os.path.exists(fgs_faa) or os.stat(fgs_faa).st_size == 0 ):

        print "Running FragGeneScan...\n"
        fraggenescan=FragGeneScan +" -genome=" + new_input_file + " -out=" + fgs_out + " -complete=1 -train=complete thread=9"
        try:
                subprocess.check_output(fraggenescan, shell=True)
        except subprocess.CalledProcessError as error:
                print error

if not os.path.exists(fgs_faa):
        log = open(error_log, "w")
        print >> log, "Generation of protein sequences failed!"
        sys.exit("Generation of protein sequences failed!Please make sure your sequence file is not empty")

# Searches all the genes against RVT-All.hmm (customized hmm model)
if not os.path.exists(myRT_RVT_tmp) or os.stat(myRT_RVT_tmp).st_size == 0:
        print "Searching the genes against RVT-All.hmm (e-value < 0.0001) ...\n"
        cmd=hmmscan + " --noali -E 0.0001 --domE 0.0001   --domtblout " +  myRT_RVT_tmp + " " + myRT +  "/Models/HMM/RVT-All.hmm " + fgs_faa
        try:
                subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as error:
                print error

if not os.path.exists(myRT_RVT_tmp) or os.stat(myRT_RVT_tmp).st_size == 0:
        log = open(error_log, "w")
        print >> log, "No RT domain was found!"
        sys.exit("No RT domain was found!")

# Extracts putative RTs
if not os.path.exists(myRT_RVT_faa) or os.stat(myRT_RVT_faa).st_size == 0 :
        print "Extracting putative RTs...\n"
        cmd="cat " + myRT_RVT_tmp + """| awk '($1!~"#") {print $4}' | uniq | while read line; do grep -A 1 ${line} """ + fgs_faa + " >> " +  myRT_RVT_faa + ";done"
        try:
                subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as error:
        	print error
if not os.path.exists(myRT_RVT_faa) or os.stat(myRT_RVT_faa).st_size == 0:
        log = open(error_log, "w")
        print >> log, "No RT domain was found!"
        sys.exit("No RT domain was found!")
# Searches putative RTs against PfamA.hmm (to remove False Positives if any)      
if not os.path.exists(myRT_RVT_cdd_pfamA) or os.stat(myRT_RVT_cdd_pfamA).st_size == 0:
        print "Searching the putative RTs against CDD domains (e-value < 0.01)...\n"
        cmd=hmmscan + " --noali -E 0.01 --domE 0.01 --domtblout " +  myRT_RVT_cdd_pfamA + " " +pfammodel + " " + myRT_RVT_faa
        try:
                subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as error:
                print error
# Verifying RTs, and extracting genomic neighborhood of verified RTs
if not os.path.exists(myRT_gn_faa) or os.stat(myRT_gn_faa).st_size == 0 :
        print "Verifying RTs, and extracting genomic neighborhood of verified RTs...\n"
        cmd="python " + myRT + "/Scripts/extract_gn.py " + output_prefix + " " + output_dir
        try:
                subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as error:
                print error

if not os.path.exists(myRT_RVT) or os.stat(myRT_RVT).st_size == 0:
        log = open(error_log, "w")
        print >> log, "No RT domain was found!"
        sys.exit("No RT domain was found!")

# seqlen
if not os.path.exists(len_file):
        cmd = """awk '$0 ~ ">" {seqid=substr($1,2); seqlen[seqid]=0;} $0 !~ ">" {seqlen[seqid]=seqlen[seqid]+length($0);} END { for (seq in seqlen) print seqlen[seq] " " seq; }' """ + new_input_file + " > " + new_input_file + ".len"
        subprocess.check_output(cmd, shell=True)

# blastp to search for similar RTs in our collection 	
if not os.path.exists(myRT_vs_ALL) or os.stat(myRT_vs_ALL).st_size == 0 :
        print "Running Blastp...\n"
        cmd="nohup " + blastp + "  -query " + myRT_RT_faa + " -db "  + myRT + "/Models/RTs-collection.faa -num_alignments 1 -evalue 1e-3 -outfmt 6 -out " +  myRT_vs_ALL 
        try:
               subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as error:
                print error

#pplacer to incorprate the phylogenetic information to improve the accuary (this step is only done for low confidence/unclassified predictions)
# Aligning the queries and the reference 
if not os.path.exists(myRT_combo) and os.path.exists(myRT_query):
       print "Incorporating phylogenetic information ..\n"
       cmd=hmmalign + " -o " + myRT_combo + " --mapali " + myRT + "/Models/myRT-FastTree2.refpkg/RVT-ref.sto " + myRT + "/Models/myRT-FastTree2.refpkg/RVT-ref.hmm " + myRT_query 
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
              flagE=1
	      # Usually this error is because of U in the sequence, so let's fix it in the next step
# Pplacer 
if not os.path.exists(myRT_jplace) and os.path.exists(myRT_combo):
       cmd=pplacer + " -c " + myRT + "/Models/myRT-FastTree2.refpkg/ " + myRT_combo + " -o " + myRT_jplace 
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                # Usually this error is because of U in the sequence, so let's fix it in the next step
                flagE=1
# removing the error-causing U in the sequence, and running pplacer again
if flagE==1:
       cmd="sed -i '/^[A-Z]/s/U//g'" + " " + myRT_query
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print error
       cmd=hmmalign + " -o " + myRT_combo + " --mapali " + myRT + "/Models/myRT-FastTree2.refpkg/RVT-ref.sto " + myRT + "/Models/myRT-FastTree2.refpkg/RVT-ref.hmm " + myRT_query
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print error
       cmd=pplacer + " -c " + myRT + "/Models/myRT-FastTree2.refpkg/ " + myRT_combo + " -o " + myRT_jplace
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print error

# If you wanna have a tree with the placed queries uncomment this 
#if not os.path.exists(myRT_tree) and os.path.exists(myRT_jplace):
#       cmd=guppy + " tog " + myRT_jplace + " -o " myRT_tree
#       try:
#           subprocess.check_output(cmd, shell=True)
#       except subprocess.CalledProcessError as error:
#                print error 

# Parsing pplacer output (.jplace) (prediction based on placement in the phylogenetic tree)
if (not os.path.exists(myRT_pred) or os.stat(myRT_pred).st_size == 0 ) and os.path.exists(myRT_jplace):
       cmd="Rscript " + JplaceParser + " " + myRT_jplace + " > " + myRT_pred + " 2>>myRT-R.log" 
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print error
# Searching for other domains in the putative RT, and domains in the flanking genes
print "Annotating the genomic neighborhood of RTs ...\n"
if not os.path.exists(myRT_gn_cdd_pfamA) or os.stat(myRT_gn_cdd_pfamA).st_size == 0:
       cmd=hmmscan + " --noali -E 0.001 --domE 0.001 --domtblout " +  myRT_gn_cdd_pfamA + " " + pfammodel + " " + myRT_gn_faa
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print error
# Combining hmmscan results, and pplacer (if used)
if not os.path.exists(myRT_gn_dom): 
       print "Preparing the domain files ...\n"     
       if os.path.exists(myRT_pred): 
       	cmd="python " + myRT + "/Scripts/hmm2dom_myRT.py " + myRT_gn_list + " " +  myRT_gn_cdd_pfamA + " " + myRT_RVT + " " + myRT_gn_dom + " " + myRT_pred
       else:
       	cmd="python " + myRT + "/Scripts/hmm2dom_myRT.py " + myRT_gn_list + " " +  myRT_gn_cdd_pfamA + " " + myRT_RVT + " " + myRT_gn_dom
       try:
           subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print error
# Preparing a gff file (RTs and flanking genes and domains) for visualization, and one containing the RTs and their predicted class  
# Also makes a list of domains in the genomic neighborhood of each RT
if not os.path.exists(myRT_gn_gff):
       print "Preparing the gff files ...\n"
       cmd="python " + myRT + "/Scripts/myRT2gff.py " + output_prefix + " " + output_dir
       try:
                subprocess.check_output(cmd, shell=True)
       except subprocess.CalledProcessError as error:
                print error
if os.path.exists(error_log):
	log.close()
cmd="touch " + Done_file
subprocess.check_output(cmd, shell=True)
print "Done!\n"
