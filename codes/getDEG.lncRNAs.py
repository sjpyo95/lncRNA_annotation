import sys, os
import module as mdl

def getDegmtxfile(infile):
	filein = open(infile, 'r')
	lines = filein.readlines(); filein.close()
	samples = lines[0].strip().split('\t')[5:]
	expdic = dict()
	for i in range(1, len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		id = tmp[0]; symb = tmp[1]; type = tmp[2]; genetype = tmp[3]
		logfc = tmp[4]
		exps = tmp[5:]
		if 'PAR_Y' in id: continue
		if 'lncRNA' not in type: continue
		if not expdic.has_key((id, symb, type,genetype,logfc)):
			expdic[(id, symb, type,genetype,logfc)] = dict()

		for j in range(len(exps)):
			sample = samples[j]; exp = float(exps[j])
			expdic[(id, symb, type,genetype,logfc)][sample] = exp

	return expdic,samples

degmtxdir = sys.argv[1]
outputdir = sys.argv[2]
filenames = filter(lambda x: '.txt' in x, os.listdir(degmtxdir))

for i in range(len(filenames)):
	filename = filenames[i]
	sample = filename.split('.degs')[0]
	infile = degmtxdir + filename
	degdic,samples = getDegmtxfile(infile)
	outfile = open(outputdir + sample + '.lncRNA_degs.mtx.txt', 'w')
	outfile.write('ID\tsymbol\ttype\tgenetype\tlogFC\t'+'\t'.join(samples))
	lines = ''
	degN = 0
	for gene in degdic.keys():
		samexps = degdic[gene]
		line = '\t'.join(gene)
		for sam in samples:
			line += '\t'+str(samexps[sam])
		lines += '\n'+line
		degN += 1
	print sample, degN
	outfile.write(lines)
	outfile.close()

