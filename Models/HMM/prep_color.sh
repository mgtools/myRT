rm -f color.txt
rm -f colors.txt
echo "TREE_COLORS" >> colors.txt
echo "SEPARATOR SPACE" >> colors.txt
echo "DATA" >> colors.txt
grep ">" June3-cp.fst | sed 's/>//g' > color.txt 
cat color.txt | awk '($1~/_GII/){print $1,"range","#0033ff"}' >> colors.txt
cat color.txt | awk '($1~/_CAS/){print $1,"range","#ff6700"}' >> colors.txt
cat color.txt | awk '($1~/_DGR/){print $1,"range","#ff0033"}' >> colors.txt
cat color.txt | awk '($1~/_RTN/){print $1,"range","#33ff00"}' >> colors.txt
cat color.txt | awk '($1~/_AbiA/){print $1,"range","#660099"}' >> colors.txt
cat color.txt | awk '($1~/_AbiK/){print $1,"range","#CC66FF"}' >> colors.txt
cat color.txt | awk '($1~/_AbiP2/){print $1,"range","#9900ff"}' >> colors.txt
cat color.txt | awk '($1~/_UG/){print $1,"range","#330033"}' >> colors.txt
cat color.txt | awk '($1~/_G2L4/){print $1,"range","#66ffff"}' >> colors.txt
cat color.txt | awk '($1~/_CL/){print $1,"range","#ffcc66"}' >> colors.txt
