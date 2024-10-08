import sys, os
import module as mdl

matrixfile = sys.argv[1]
metafile = sys.argv[2]
groupNum = int(sys.argv[3])
outfile = sys.argv[4]

def getMeta(infile, groupNum):
	mfile = open(infile, 'r')
	lines = mfile.readlines(); mfile.close()
	meta = dict(); groups = []
	for i in range(len(lines)):
		line = lines[i].strip()
		if not line.startswith('sample'):
			col = line.split('\t')
			sample = col[1]; group = col[groupNum+2]
			if not meta.has_key(group):
				meta[group] = []
			meta[group].append(sample)
			if not group in groups:
				groups.append(group)
	return meta, groups

meta,groups= getMeta(metafile, groupNum)
matrix, samples = mdl.getMatrix(matrixfile)

newM = dict()
for gene in matrix.keys():
	if not newM.has_key(gene):
		newM[gene] = {}
	for group in meta.keys():
		sams = meta[group]
		exps = map(lambda x: float(matrix[gene][x]), sams)
		avgExp = sum(exps)/len(exps)
		newM[gene][group] = avgExp

mdl.writeMatrix(newM, groups, outfile)
