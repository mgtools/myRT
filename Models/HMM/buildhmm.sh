# to update one model:  bash buildhmm.sh RVT-DGRs (don't add .fst)
# to update all models: cat list | while read line ; do bash buildhmm.sh ${line} ; done 
file=$1
muscle -in ${file} -out ${file}.fst
./reformat.pl  fas sto ${file}.fst .sto
../../bin/hmmer-3.2/src/hmmbuild -n ${file} ${file}.hmm ${file}.sto
rm -f RVT-All.hmm*
cat *.hmm > RVT-All.hmm
../../bin/hmmer-3.2/src/hmmpress RVT-All.hmm
chmod a+rx RVT-All.hmm*
