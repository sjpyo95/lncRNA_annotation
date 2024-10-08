import sys, os
import module as mdl

def get_infoMatrix(filename, cutoff):
	filein = open(filename)
	lines= filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[4:]
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id  = line[0]; symb = line[1]; type = line[2]; genetype = line[3]
		exps = line[4:]
#		if type == 'Known_lncRNA' or type == 'Novel_lncRNA': type = 'lncRNA'
		for x in range(len(samples)):
			sample = samples[x]; exp = float(exps[x])
			if exp < cutoff: continue

			if not geneDic.has_key(sample):
				geneDic[sample] = dict()
			if not geneDic[sample].has_key(type):
				geneDic[sample][type] = dict()

			if not geneDic[sample][type].has_key(genetype):
				geneDic[sample][type][genetype] = 0
			geneDic[sample][type][genetype] += 1

	return geneDic, samples

mtx_file = sys.argv[1]
out_file = sys.argv[2]
expcutoff = float(sys.argv[3])
mtx, samples = get_infoMatrix(mtx_file, expcutoff)


outfile = open(out_file, 'w')
outfile.write('sample\ttype\tgenetype\tcount\n')

for sample in samples:
	typedic = mtx[sample]

	for type in typedic.keys():
		genetypedic = typedic[type]
		for genetype in genetypedic.keys():
			count = genetypedic[genetype]
			outfile.write(sample+'\t'+type+'\t'+genetype+'\t'+str(count)+'\n')
outfile.close()
