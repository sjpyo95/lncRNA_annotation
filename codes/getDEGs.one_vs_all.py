import sys, os
import module as mdl


def getDEG(infile):
	filein = open(infile, 'r')
	sample = infile.split('/')[-1].split('.deseq2')[0]
	lines = filein.readlines(); filein.close()
	degs = dict()
	for i in range(1, len(lines)):
		line = lines[i].strip()
		col = line.split('\t')
		id = col[0]; lgfc = float(col[2]); padj = float(col[-1])
		if abs(lgfc) >= 2 and padj <= 0.05:
			 degs[id] = lgfc
	return degs

def getMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[4:]
	print len(samples)
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; exp = line[4:]; type = line[2]; genetype = line[3]
		if 'PAR_Y' in id: continue
		if not geneDic.has_key((id,symb,type,genetype)):
			geneDic[(id,symb,type,genetype)] = dict()
		for x in xrange(len(exp)):
			sample = samples[x]; express = float(exp[x])
			geneDic[(id,symb,type,genetype)][sample] = express
	return geneDic, samples

degdir = sys.argv[1]
mtx_file = sys.argv[2]
outputdir = sys.argv[3]
expcut = float(sys.argv[4])
if not os.path.exists(outputdir): os.makedirs(outputdir)
mtx, samples = getMatrix(mtx_file)

degfiles = filter(lambda x: '.txt' in x and 'deseq2' in x, os.listdir(degdir))

degdic = dict()

for i in range(len(degfiles)):
	degfile = degdir + degfiles[i]
	celltype = degfile.split('/')[-1].split('.deseq2')[0]
	degids = getDEG(degfile)
	outfile = open(outputdir + celltype + '.degs.tpm.mtx.txt', 'w')
	outfile.write('ID\tsymbol\ttype\tgenetype\tlogFC\t'+'\t'.join(samples))
	lines = ''
	degN = 0
	for gene in mtx.keys():
		geneid = gene[0]; symbol = gene[1]; type = gene[2]; genetype = gene[3]
		expdic = mtx[gene]
		if geneid not in degids.keys(): continue
		if float(expdic[celltype]) < expcut: continue
		lgfc = str(degids[geneid])
		line = '\t'.join(gene)+'\t'+lgfc
		for sample in samples:
			line += '\t'+str(expdic[sample])
		lines += '\n'+line
		degN += 1
	print celltype, degN
	outfile.write(lines)
	outfile.close()
