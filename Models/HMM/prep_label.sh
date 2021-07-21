rm -f color.txt
rm -f data.txt
echo "DATASET_TEXT" >> data.txt
echo "SEPARATOR COMMA" >> data.txt
echo "DATASET_LABEL,label" >> data.txt
echo "COLOR,#ff0000" >> data.txt
echo "DATA" >> data.txt
cat Idd4 > color.txt 
cat color.txt | awk '($1~/_GII/){print $1",-1,#0033ff,bold,3,0"}' >> data.txt
cat color.txt | awk '($1~/_CAS/){print $1",-1,#ff6700,bold,3,0"}' >> data.txt
cat color.txt | awk '($1~/_DGR/){print $1",-1,#ff0033,bold,3,0"}' >> data.txt
cat color.txt | awk '($1~/_RTN/){print $1",-1,#33ff00,bold,3,0"}' >> data.txt
cat color.txt | awk '($1~/_AbiA/){print $1",-1,#660099,bold,3,0"}' >> data.txt
cat color.txt | awk '($1~/_AbiK/){print $1",-1,#CC66FF,bold,3,0"}' >> data.txt
cat color.txt | awk '($1~/_AbiP2/){print $1",-1,#9900ff,bold,3,0"}' >> data.txt
cat color.txt | awk '($1~/_UG/){print $1",-1,#330033,bold,3,0"}' >> data.txt
cat color.txt | awk '($1~/_G2L4/){print $1",-1,#66ffff,bold,3,0"}' >> data.txt
cat color.txt | awk '($1~/_CL/){print $1",-1,#ffcc66,bold,3,0"}' >> data.txt
