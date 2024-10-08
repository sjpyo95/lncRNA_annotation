import sys
import module as mdl
matrixfile = sys.argv[1]
gtfFile = sys.argv[2]
outfile = open(sys.argv[3], 'w')
geneids = list(set([x.geneid().split('-')[-1].split('.')[0] for x in sum(mdl.getGene(gtfFile).values(), [])]))

def getMatrix(filename):
	filein = open(filename, 'r')
	lines = filein.readlines()
	samples = lines[0].strip().split('\t')[2:]
	geneDic = dict()
	for i in xrange(1, len(lines)):
		line = lines[i].strip()
		geneid = line.split('\t')[0];
		geneDic[geneid] = line
	return samples, geneDic

samples, matrix = getMatrix(matrixfile)
outfile.write('ID\tsymbol\t'+'\t'.join(samples)+'\n')
for g in matrix.keys():
	line = matrix[g]
	geneid = g.split('-')[-1].split('.')[0]
	if geneid in geneids:
		outfile.write(line+'\n')
outfile.close()
