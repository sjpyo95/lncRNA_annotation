import sys, os
import module as mdl

novelfile = sys.argv[1]
reffile = sys.argv[2]
outfile = sys.argv[3]
novelgtf = mdl.getGtf(novelfile)
refgtf = mdl.getGtf(reffile)
newGtf = dict()
comm1 = 0; comm2= 0; allcomm = 0
for chr in novelgtf.keys():
	trxs = novelgtf[chr]
	reftrxs = refgtf[chr]
	antisense = filter(lambda x: x.genetype() == 'antisense', trxs)
	otherTrxs = filter(lambda x: x.genetype() != 'antisense', trxs)
	newGtf[chr] = []
	for anti in antisense:
		ovrefs = filter(lambda x: x.overlaps3(anti) == True, reftrxs)
		exonStr = map(lambda x: (x.start(), x.end()), anti.exons())
		exonStr = mdl.comp(exonStr, lambda x: x[0], True)
		for ov in ovrefs:
			ovexonStr = map(lambda x: (x.start(), x.end()), ov.exons())
			ovexonStr = mdl.comp(ovexonStr, lambda x: x[0], True)
			commonStr = list(set(exonStr).intersection(ovexonStr))
			if len(commonStr) >= 1:
#				newGtf[chr].append(anti)
				comm1 += 1
				continue
			elif len(commonStr) == 2:
				comm2 += 1
		newGtf[chr].append(anti)
print comm1
print comm2
mdl.writeGtf(newGtf, outfile)
