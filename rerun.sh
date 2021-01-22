str=$1
col=$2
rm -f ${col}/${str}/${str}.done

# If you already have FragGeneScan results, lets not redo it
if test -f "${col}/${str}/${str}-FGS.gff"; then
./myRT-local.py ${col} ${str} ${col}/${str}/${str}.fna ${col}/${str}/${str}-FGS.gff

# If ncbi gff file is already in the folder let's use it and skip FragGeneScan:
elif test -f "${col}/${str}/${str}.gff"; then
./myRT-local.py ${col} ${str} ${col}/${str}/${str}.fna ${col}/${str}/${str}.gff
else

# If no gff file was found, let's use FragGeneScan for gene prediction  (FASTEST)
./myRT-local.py ${col} ${str} ${col}/${str}/${str}.fna
fi
