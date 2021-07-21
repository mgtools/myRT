# Run as ./buildhmm.sh RVT-DGRs
# No need to run if you haven't updated the models
# Make sure the hmmer directory is correct, you may need to add the path to myRT to the hmmer path
file=$1
muscle -in ${file} -out ${file}.fst
./reformat.pl  fas sto ${file}.fst .sto
# You can use fasta2sto.py instead of reformat.pl
myRT/bin/hmmer-3.2/src/hmmbuild -n ${file} ${file}.hmm ${file}.sto
rm -f RVT-All.hmm*
cat *.hmm > RVT-All.hmm
myRT/bin/hmmer-3.2/src/hmmpress RVT-All.hmm
chmod a+rx RVT-All.hmm*
