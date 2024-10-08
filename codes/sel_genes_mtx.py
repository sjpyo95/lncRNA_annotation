import sys, os
import module as mdl

infile = sys.argv[1]
genes = sys.argv[2]
genelist = genes.split(',')
outfile = sys.argv[3]
mtx, samples = mdl.getinfoMatrix(infile)
newMtx = dict()
for key in mtx.keys():
	symb = key[1]
	val = mtx[key]
	if symb not in genelist: continue
	newMtx[key] = val

mdl.writeinfoMatrix(newMtx, samples, outfile)
