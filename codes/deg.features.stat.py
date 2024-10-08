import sys, os
import module as mdl

def getdegMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[5:]
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; lgfc = float(line[4]);#padj = float(line[5])
		exp = line[5:]; type = line[2]; genetype = line[3]
		if 'PAR_Y' in id: continue
		elif lgfc < 2 : continue
#		elif padj > 0.05 : continue
		if not geneDic.has_key((id,symb,type,genetype,lgfc)):
			geneDic[(id,symb,type,genetype,lgfc)] = dict()
		for x in xrange(len(exp)):
			sample = samples[x]; express = float(exp[x])
			geneDic[(id,symb,type,genetype,lgfc)][sample] = express
	return geneDic, samples


inputdir = sys.argv[1]
`
degfiles = filter(lambda x: '.txt' in x, os.listdir(inputdir))

for i in range(len(degfiles)):
	degfile = inputdir + degfiles[i]
	groupname = degfiles[i].split('.degs')[0]

