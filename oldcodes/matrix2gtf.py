import sys
import module as mdl
matrixfile = sys.argv[1]
gtffile = sys.argv[2]
outfile = sys.argv[3]
def parseMatrix(infile):
	matrix = open(infile, 'r')
	lines = matrix.readlines(); matrix.close()
	ids  = []
	for i in range(1, len(lines)):
		line = lines[i].strip()
		geneid = line.split('\t')[0]
		ids.append(geneid)
	return ids

gtfdic = mdl.getGtf(gtffile)
newgtf = dict()
ids = parseMatrix(matrixfile)
geneids = map(lambda x: x.geneid().split('.')[0], sum(gtfdic.values(), []))
#geneid2 = map(lambda x: x.geneid(), sum(gtfdic.values(), []))
#print len(geneids)
#print len(list(set(geneids)))
#print len(list(set(geneid2)))
#test = filter(lambda x: x not in geneids, ids)
#print test
#exit()
for chrom in gtfdic.keys():
	newtrxs = []
	trxs = gtfdic[chrom]
	for trx in trxs:
		geneid = trx.geneid()
		sid = geneid
		if sid in ids:
			if trx not in newtrxs:
				newtrxs.append(trx)
	
	newtrxs = filter(lambda x: x.geneid().split('.')[0] in ids, trxs)
	newgtf[chrom] = newtrxs

mdl.writeGtf(newgtf, outfile)
