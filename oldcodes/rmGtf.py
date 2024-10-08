import sys
import module as mdl
gtf = sys.argv[1]
rm = sys.argv[2]
outfile = sys.argv[3]
gtfdic = mdl.getGtf(gtf)
rmdic = mdl.getGtf(rm)
newgtf = dict()
for chrom in rmdic.keys():
	trxs = gtfdic[chrom]
	rmgene = list(set(map(lambda x: x.geneid(), rmdic[chrom])))
	ftrxs = list(set(filter(lambda x: x.geneid() not in rmgene, trxs)))
	newgtf[chrom] = ftrxs
mdl.writeGtf(newgtf, outfile)


