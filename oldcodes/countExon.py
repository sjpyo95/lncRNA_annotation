import sys
import module as mdl
infile = sys.argv[1]
outfile = open(sys.argv[2], 'w')

gtf = mdl.getGtf(infile)
trxs = sum(gtf.values(), [])
exonNdic = dict()
typedic = dict()
for trx in trxs:
	exons = trx.exons()
	exonN = len(exons)
	gtype = trx.genetype()
	if not typedic.has_key(gtype):
		typedic[gtype] = dict()
	if exonN < 5:
		if not exonNdic.has_key(exonN):
			exonNdic[exonN] = 0
		exonNdic[exonN] += 1
		if not typedic[gtype].has_key(exonN):
			typedic[gtype][exonN] = 0
		typedic[gtype][exonN] += 1
	else:
		if not exonNdic.has_key(5):
			exonNdic[5] = 0
		exonNdic[5] +=1
		if not typedic[gtype].has_key(5):
			typedic[gtype][5] = 0
		typedic[gtype][5] +=1

outfile.write('type\texonN\tnum\n')
for type in typedic.keys():
	exonNdic = typedic[type]
	for exonN in exonNdic.keys():
		num = exonNdic[exonN]
		outfile.write(type+'\t'+str(exonN)+'\t'+str(num)+'\n')
for exonN in exonNdic.keys():
	num = exonNdic[exonN]
	outfile.write('total\t'+str(exonN)+'\t'+str(num)+'\n')
outfile.close()
