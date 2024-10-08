import sys, os
import module as mdl

def getCorrelation(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	tgenes = lines[0].split('\t')
	pcordic = dict(); ncordic = dict()
	for i in range(1,len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		cgene = tmp[0]; cor = tmp[1:]
		for t in range(len(tgenes)):
			tg = tgenes[t]; corScore = float(cor[t])
			if corScore >= 0.90: 
				if not pcordic.has_key(tg):
					pcordic[tg] = []
				pcordic[tg].append(cgene)
			elif corScore <= -0.90:
				if not ncordic.has_key(tg):
					ncordic[tg] = []
				ncordic[tg].append(cgene)
	return pcordic, ncordic

mtx_file = sys.argv[1]
cor_file = sys.argv[2]
target_genes = sys.argv[3]
out_file = sys.argv[4]


mtx, samples = mdl.getinfoMatrix(mtx_file)

pcordic. ncordic = getCorrelation(cor_file)

if target_genes == 'all': targets = cordic.keys()
else: targets = target_genes.split(',')
print 'targets genes No.:', len(targets)

outfile = open(out_file, 'w')
outfile.write('ID\tsymbol\ttype\tgenetype\t'+'\t'.join(samples))
lines = ''

checkgenes = []
for target in targets:
	corgenes = cordic[target]
	for gene in mtx.keys():
		symbol = gene[1]
		if symbol not in corgenes: continue
		if symbol in checkgenes: continue
		checkgenes.append(symbol)

		expdic = mtx[gene]
		line = '\t'.join(gene)
		for sample in samples:
			line += '\t'+str(expdic[sample])
		lines += '\n'+line
	outfile.write(lines)
outfile.close()


