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
		id = line[0]; symb = line[1]; exp = line[5:]; lgfc = float(line[4]); type = line[2]; genetype = line[3]
		if 'PAR_Y' in id: continue
		if abs(lgfc) < 1: continue
		if genetype == 'intervening': continue
		if type == 'PCG': continue
		if not geneDic.has_key((id,symb,type,genetype,'DEG')):
			geneDic[(id,symb,type,genetype,'DEG')] = dict()
		for x in xrange(len(exp)):
			sample = samples[x]; express = float(exp[x])
			geneDic[(id,symb,type,genetype, 'DEG')][sample] = express
	return geneDic, samples

def getgmtMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[5:]
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		pathway = line[0]; id = line[1]; symb = line[2]; type = line[3]; genetype = line[4]
		exps =  line[5:]
		if 'PAR_Y' in id: continue
		if not geneDic.has_key((id,symb,type,genetype,pathway)):
			geneDic[(id,symb,type,genetype,pathway)] = dict()
		for x in xrange(len(exps)):
			sample = samples[x]; expr = float(exps[x])
			geneDic[(id,symb,type,genetype,pathway)][sample] = expr
	return geneDic, samples

inputdir = sys.argv[1]
gmtmtx_file = sys.argv[2]
outputdir = sys.argv[3]

if not os.path.exists(outputdir): os.makedirs(outputdir)

gmtmtx, samples = getgmtMatrix(gmtmtx_file)
gmtname = gmtmtx_file.split('/')[-1].split('.')[0]
mtxfiles = filter(lambda x: '.txt' in x, os.listdir(inputdir))
skipedcells = []
for i in xrange(len(mtxfiles)):
	mtxfile = inputdir + mtxfiles[i]
	celltype = mtxfiles[i].split('.deg')[0]
	mtx, samples = getdegMatrix(mtxfile)
	if len(mtx) == 0: skipedcells.append(celltype); continue	
	print celltype, len(mtx.keys())

	mtx.update(gmtmtx)

	outfile = open(outputdir + celltype + '_deg_lncRNAs.'+gmtname + '.mtx.txt', 'w')
	outfile.write('ID\tsymbol\ttype\tgenetype\tpathway\t'+'\t'.join(samples))
	lines = ''
	for gene in mtx.keys():
		expdic = mtx[gene]
		line = '\t'.join(gene)
		for sample in samples:
			line += '\t' + str(expdic[sample])
		lines += '\n' + line
	outfile.write(lines)
	outfile.close()
print 'No DEG lncRNAs logFC >= 2:', skipedcells
