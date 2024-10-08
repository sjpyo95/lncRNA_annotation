import sys, os
import module as mdl

gtf1 = mdl.getGtf(sys.argv[1])
gtf2 = mdl.getGtf(sys.argv[2])
outfile = sys.argv[3]
geneids1 = map(lambda x: x.geneid(), sum(gtf1.values(), []))
geneids2 = list(set(map(lambda x: x.geneid(), sum(gtf2.values(), []))))

#include = filter(lambda x: x in geneids2, geneids1)
#exclude = filter(lambda x: x not in geneids1, geneids2)
#print exclude, len(exclude)
newgtf = dict()
for chrom in gtf1.keys():
	trxs = gtf1[chrom]
	newtrxs = filter(lambda x: x.geneid() in geneids2, trxs)
#	print map(lambda x: x.geneid(), newtrxs)
	newgtf[chrom] = newtrxs
mdl.writeGtf(newgtf, outfile)
