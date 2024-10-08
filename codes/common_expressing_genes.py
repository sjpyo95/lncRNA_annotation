import sys, os
import module as mdl
import numpy as np

matrix_file = sys.argv[1]
rank_list = sys.argv[2]
out_file1 = sys.argv[3]
out_file2 = sys.argv[4]
mtx, samples = mdl.getMatrix(matrix_file)

def getlist(infile):
	filein = open(infile, 'r')
	lines = filein.readlines(); filein.close()
	sample = lines[0].split('\t')[-1]
	rankdic = dict()
	rankdic[sample] = []
	for i in range(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symbol = line[1]
		exp = float(line[2])
		if exp < 1: continue
		rankdic[(id,symbol)] = exp
	return rankdic, sample

ranks, sample = getlist(rank_list)

#common
common = dict()
spec = dict()
for idsymb in mtx.keys():
	exps = mtx[idsymb].values()
	if len(filter(lambda x: float(x) >= 1, exps)) >= len(samples)/2.0:
		avgExp = np.mean(exps)
		common[idsymb] = mtx[idsymb]
	elif idsymb in ranks.keys() and len(filter(lambda x: mtx[idsymb][x] > ranks[idsymb], samples)) <= len(samples)/4:
		spec[idsymb] = ranks[idsymb]

outfile=  open(out_file1, 'w')
outfile.write('ID\tsymbol\t'+'\t'.join(samples)+'\n')
for idsymb in common.keys():
	exps = common[idsymb]
	outfile.write(idsymb[0]+'\t'+idsymb[1])
	for sam in samples:
		outfile.write('\t'+str(exps[sam]))
	outfile.write('\n')
outfile.close()



outfile = open(out_file2, 'w')
outfile.write('ID\tsymbol\t'+sample)
for idsymb in spec.keys():
	id = idsymb[0]; symb = idsymb[1]
	outfile.write(id+'\t'+symb+'\t'+str(spec[idsymb])+'\n')
outfile.close()
	
	




