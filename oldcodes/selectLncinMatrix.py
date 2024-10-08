import sys
import module as mdl
matrixfile = sys.argv[1]
gtffile = sys.argv[2]
genetype = sys.argv[3]	#protein_coding, antisense, lincRNA, ... all
outfile = sys.argv[4]
gtf = sum(mdl.getGene(gtffile).values(),[])
#genes = filter(lambda x: x.genetype() == genetype, gtf)
#print len(genes)
matrix = mdl.getMatrix(matrixfile)
samples = matrix.values()[0].keys()
print samples
newMatrix = dict()
for idsym in matrix.keys():
	id = idsym[0]; symb = idsym[1]
	expDic = matrix[idsym]
	if id in [x.geneid().split('.')[0] for x in gtf]:
		if not newMatrix.has_key((id, symb)):
			newMatrix[(id,symb)] = dict()
		newMatrix[(id,symb)] = expDic

mdl.writeMatrix(newMatrix, outfile)
