import sys
import module as mdl

infile = sys.argv[1]
#refFile = sys.argv[2]
outfile = sys.argv[2]
ingtf = mdl.getGtf(infile)
#ref = mdl.getGtf(refFile)

newGtf = dict()

for chr in ingtf.keys():
	intrxs = filter(lambda x: 'Mo.known' in x.geneid(), ingtf[chr])
	newGtf[chr] = intrxs
mdl.writeGtf(newGtf, outfile)

#	reftrxs = ref[chr]
#	for trx in intrxs:
#		geneid = trx.geneid().split('lncRNA-')[-1]
#		start = trx.start(); end = trx.end()
#		refgene = filter(lambda x: x.geneid().split('.')[0] == geneid and x.start() == start and x.end() == end, reftrxs)
#		if len(refgene) != 1:
##		trx.setSymb(symb)
#			print geneid
