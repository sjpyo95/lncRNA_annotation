import sys, os
import module as mdl

infile = sys.argv[1]
outfile = sys.argv[2]

matrix = mdl.getMatrix(infile)
newMatrix = dict()
for idsym in matrix.keys():
	express = matrix[idsym]
	totalcount = sum(float(x) for x in express.values())
	if totalcount == 0: continue
	elif 'MIR' in idsym[-1]: continue
#	print totalcount
	if not newMatrix.has_key(idsym):
		newMatrix[idsym] = dict()
	newMatrix[idsym] = express
mdl.writeMatrix(newMatrix, outfile)
