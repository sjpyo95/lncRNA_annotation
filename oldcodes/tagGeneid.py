import sys
import module as mdl

infile = sys.argv[1]
outfile = sys.argv[2]
tag = sys.argv[3]
gtfdic = mdl.getGtf(infile)
newgtf = dict()
for chrom in gtfdic.keys():
	gtf = gtfdic[chrom]
	if not newgtf.has_key(chrom):
		newgtf[chrom] = []
	for i in xrange(len(gtf)):
		trx = gtf[i]
		geneid = trx.geneid()
		newid = tag+geneid
		trx.setGeneid(newid)
		newgtf[chrom].append(trx)

mdl.writeGtf(newgtf, outfile)
