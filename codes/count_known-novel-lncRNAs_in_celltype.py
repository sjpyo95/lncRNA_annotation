import sys, os
import module as mdl

def get_infoMatrix(filename):
	filein = open(filename)
	lines= filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[4:]
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id  = line[0]; symb = line[1]; type = line[2]; genetype = line[3]
		exps = line[4:]
		for x in xrange(len(exps)):
			sample = samples[x]; exp = float(exps[x])
			if not geneDic.has_key(sample):
				geneDic[sample] = dict()
			if not geneDic[sample].has_key(type):
				geneDic[sample][type] = []
			if exp == 0 : continue
			geneDic[sample][type].append(exp)

	return geneDic

info_mtx_file = sys.argv[1]
out_file = sys.argv[2]
outfile = open(out_file, 'w')
outfile.write('cell-type\ttype\tcount\n')

mtx = get_infoMatrix(info_mtx_file)
for sample in mtx.keys():
	typedic = mtx[sample]
	for type in typedic.keys():
		expN = len(typedic[type])
		outfile.write(sample+'\t'+type+'\t'+str(expN)+'\n')
