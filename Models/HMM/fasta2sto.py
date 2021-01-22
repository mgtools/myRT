#!/usr/bin/env python
import sys

if len(sys.argv) < 3:
	sys.exit("fasta2sto.py fasta-input sto-output")

seqid, alignseq = [], []
inf = open(sys.argv[1], "r")
for aline in inf:
	subs = aline.strip().split()
	if aline[0] == '>':
		seqid.append(subs[0][1:])
		alignseq.append("")
	else:
		alignseq[-1] += aline.strip()
inf.close()

out = open(sys.argv[2], "w")
print >>out, "# STOCKHOLM 1.0"
print >>out, "#=GF SQ ", str(len(seqid))

for idx in range(len(seqid)):
	print >>out, seqid[idx], alignseq[idx]
print >>out, "//"
out.close()
