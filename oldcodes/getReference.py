import sys
import module as mdl
matrix = open(sys.argv[1], 'r')
gtf = sys.argv[2]
outfile = sys.argv[3]
gtfdic = mdl.getGtf(gtf)
lines = matrix.readlines(); matrix.close()
mgenes = []
for i in xrange(1, len(lines)):
	line = lines[i].strip().split('\t')
	geneid = line[0]
	mgenes.append(geneid)

#allgenes = [x.geneid().split('.')[0] for x in sum(gtfdic.values(), [])]

#test = filter(lambda x: x in allgenes, mgenes)
#print len(test)
#print test
#print len(test)
newGtf = dict()
for chrom in gtfdic.keys():
	trxs = gtfdic[chrom]
	newTrxs = list(set(filter(lambda x: x.geneid().split('.')[0] in mgenes and not x.geneid().split('.')[-1].find('PAR') > -1, trxs)))
	
	newGtf[chrom] = newTrxs
mdl.writeGtf(newGtf, outfile)
