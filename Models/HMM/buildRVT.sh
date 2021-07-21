
# First search the RT sequences against pfam-A.hmm or cdd-PfamA.hmm 
# You can also check against RVT_1.hmm but you need to make sure that you are only searching RT sequencing (it's kind of tricky to search against only one model because DNA binding domains can also be similar to RVT_1, so it's best to search against all and grab RVT_1 hits)

# You can do this for RTs from each group (e.g. DGRs) seperately, or you can do it for All and then separatd based on the class

#myRT/bin/hmmer-3.2/src/hmmscan --noali -E 0.0001 --domE 0.0001   ALL_vs_pfam ../../models/Pfam/Pfam-A.hmm RTs-collection.faa

# Then create the RVT_1 motif sequences

#cat All-vs-pfam | grep "RVT_1" | awk '($1!~"#"){print $4,$18,$19}' | awk '!seen[$1]++' | while read id start stop ; do echo ">"$id"_"$start"_"$stop >> RVT-seq.faa ; grep -A 1 ${id} RVT-seq.faa | tail -n +2 | cut -c $start-$stop >> RVT-seq.faa ; done

# Then you need to check the sequences and filter out the ones that are too short, ....


# Then use Muscle to create one multiple sequence alignment for All RVT_1 motif sequences

#Muscle -in RVT-seq-filtered.faa -out RVT-seq-filtered.fst

# Then use Fasttree to create the phylogenetic information
#FastTreeMP RVT-seq-filtered.fst > RVT-seq-filtered.tre
# This gives you one replication
# To create the bootstrap values 
#Use python seqIO or other tools to convert fasta to phylip

## Python script (open python and paste the script, or create a .py script and run)
#from Bio import SeqIO
#records = SeqIO.parse("RVT-seq-filtered.fst", "fasta")
#count = SeqIO.write(records, "RVT-seq-filtered.phy", "phylip")
#print("Converted %i records" % count)
# End of python script


# Then use seqboot to create boostrap alignments
# When you get the seqboot results
#FastTreeMP -n 100 -log RVT-seq-filtered-bp100.log RVT-seq-filtered-bp100.phylip > RVT-seq-filtered-bp100.tre

# To use the tree for pplacer: 

# Use the filtered RVT_1 motif sequences as reference

# cp RVT-seq-filtered.faa RVT-ref.faa


#cd taxtastic/
#virtualenv taxtastic-env
#source taxtastic-env/bin/activate
#docker pull nghoffman/taxtastic:0.8.3
#docker run -it --rm nghoffman/taxtastic:0.8.3 taxit --version
#cd docker
#taxit create -l RT -P myRT-FastTree.refpkg --aln-fasta RVT-ref.fst --tree-file RVT-ref.tre --tree-stats  RVT-ref.log --author fsharifi --description myRT 

#To build the reference hmm
# first create RVT-ref.sto using reformat.pl or fasta2sto.py

#myRT/bin/hmmer-3.2/src/hmmbuild RVT-ref.hmm RVT-ref.sto

#Check the myRT-local.py to see how to use the reference tree and reference profiles (RVT-ref.hmm)


