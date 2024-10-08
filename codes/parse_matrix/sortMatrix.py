import sys, os
import module as mdl

matrixfile = sys.argv[1]
metafile = sys.argv[2]
outfile = sys.argv[3]
matrix, samples = mdl.getMatrix(matrixfile)
meta = open(metafile, 'r')
lines = meta.readlines(); meta.close()
metdic = dict()
samOrders = []
for i in xrange(len(lines)):
	line = lines[i].strip()
	if not line.startswith('sample'):
		col = line.split('\t')
		sample = col[1]; rawname = col[0]
		metdic[rawname] = sample
		samOrders.append(sample)

newMatrix = dict()
expressingGenes = 0; no_expressingGenes = 0
for idsymb in matrix.keys():
	exps = matrix[idsymb]
#	if len(filter(lambda x: float(x) >= 1, exps.values())) < 1: no_expressingGenes += 1; continue
#	expSum = sum(map(lambda x: float(x), exps.values()))
#	if expSum == 0: no_expressingGenes += 1; continue
	expressingGenes += 1
	if not newMatrix.has_key(idsymb):
		newMatrix[idsymb] = dict()
	for sample in exps.keys():
		exp = exps[sample]
		name = metdic[sample]
		newMatrix[idsymb][name] = exp
print 'less than 1 TPM  expression in all samples:',no_expressingGenes
mdl.writeMatrix(newMatrix, samOrders, outfile)
