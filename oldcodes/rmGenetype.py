import sys
import module as mdl

matrixfile = sys.argv[1]
gtffile = sys.argv[2]
#type = sys.argv[3]
outfile = sys.argv[3]
rmOrselect = sys.argv[4]	#remove or select
matrix = mdl.getMatrix(matrixfile)
genes =  map(lambda x: x.geneid().split('.')[0], sum(mdl.getGene(gtffile).values(),[]))
#genes_type = map(lambda x: x.geneid().split('.')[0], filter(lambda x: x.genetype() == type, genes))
newMatrix = dict()
for idsym in matrix.keys():
	id = idsym[0]; symb = idsym[1]
	express = matrix[idsym]
	if rmOrselect == 'rm' or rmOrselect == 'remove':
		if id not in genes:
			newMatrix[idsym] = express
	elif rmOrselect == 'select':
		if id in genes:
			newMatrix[idsym] = express
mdl.writeMatrix(newMatrix, outfile)


