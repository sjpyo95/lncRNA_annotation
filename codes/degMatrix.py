import sys,os
import module as mdl

celltype = sys.argv[1]
mtx_file = sys.argv[2]
deg_file = sys.argv[3]
out_file = sys.argv[4]

def getDEG(infile):
	degfile = open(infile, 'r')
	lines = degfile.readlines(); degfile.close()
	degdic = dict()
	for i in range(1,len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		cell = tmp[0]; geneid = tmp[1]; type = tmp[2]; expression = tmp[3]
		if not degdic.has_key(cell):
			degdic[cell] = []
		degdic[cell].append(geneid)
	return degdic

def getMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[4:]
	print len(samples)
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; exp = line[4:]; type = line[2]
		if 'PAR_Y' in id: continue
		if not geneDic.has_key((id,symb,type)):
			geneDic[(id,symb,type)] = dict()
		for x in xrange(len(exp)):
			sample = samples[x]; express = float(exp[x])
			geneDic[(id,symb,type)][sample] = express
	return geneDic, samples

mtx,samples = getMatrix(mtx_file)
degs = getDEG(deg_file)[celltype]

outfile = open(out_file, 'w')

outfile.write('ID\tsymbol\ttype\t'+'\t'.join(samples)+'\n')

for idsymb in mtx.keys():
	id = idsymb[0]
	if id in degs:
		exps = mtx[idsymb]
		outfile.write('\t'.join(idsymb))
		for sample in samples:
			outfile.write('\t'+str(exps[sample]))
		outfile.write('\n')
print samples
print celltype
print 'degs: ',len(degs)
outfile.close()
