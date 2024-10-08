import sys, os
import module as mdl
import numpy as np
degfile = sys.argv[1]
matrixfile = sys.argv[2]
contrast = sys.argv[3]	# example: BM_hi.vs.BM_low
samples = contrast.split('.vs.')
print samples
outfile = open(sys.argv[4], 'w')
pcgs = sys.argv[5]
lncs = sys.argv[6]
topN = sys.argv[7]	#int or all
program = sys.argv[8]

pcgGtf = sum(mdl.getGene(pcgs).values(), [])
lncGtf = sum(mdl.getGene(lncs).values(), [])
pcgid = map(lambda x: x.geneid().split('.')[0], pcgGtf)
#lncid = map(lambda x: x.geneid().split('.')[0], lncGtf)
lncid = []
for gene in lncGtf:
	if not 'ENSM' in gene.geneid():
		id = gene.geneid()
	else:
		id = gene.geneid().split('.')[0]
	lncid.append(id)


def getDEG(filename):
	infile = open(filename)
	lines = infile.readlines();infile.close()
	degdic = dict()
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('ID'):
			line = line.split('\t')
			if program == 'deseq':
				id = line[0]; baseMean = line[1]; log2fc = line[2]; lfcse = line[3]; stat = line[4]
				pval = line[5]; padj = line[6]
			elif program == 'edgeR':
				id = line[0]; log2fc = line[1]
				pval = line[3]; padj = line[-1]
			degdic[id] = [log2fc, padj, pval]
	return degdic
def calculateVarianceMean(samDic, group):
	var = []
	for g in group:
		if 'non' in g:
			samples = filter(lambda x: g in x and 'non' in x, samDic.keys())
		else:
			samples = filter(lambda x: g in x and not 'non' in x, samDic.keys())
#		samples = filter(lambda x: g in x, samDic.keys())
		express = map(lambda x: float(samDic[x]), samples)
#		sam_var = np.var(express)
		var.append((max(express)+1)/(min(express)+1))
	maxVar = max(var)
	return maxVar
def getAbsoluteExpression(id, matrix, samples):
	expresses = dict()
	for idsym in matrix.keys():
		if id == idsym[0]:
			express = matrix[idsym]
#			print express.keys()
			for sam in express.keys():
#				if 'non' in sam: continue
				if samples[0] in sam and not 'non' in sam:
					expresses[sam] = express[sam]
				elif samples[-1] in sam:
					expresses[sam] = express[sam]
	return expresses

matrix = mdl.getMatrix(matrixfile)
degs = getDEG(degfile)
outfile.write('ID\tsymbol\ttype\tlog2FoldChange\tpadj\t')

degFilt_up= []
degFilt_down = []
vars = []

for id in degs.keys():
	log2fc = degs[id][0]; pval = degs[id][2]; padj = degs[id][1]
	if not id in map(lambda x: x[0], matrix.keys()): continue
	symb = filter(lambda x: id == x[0], matrix.keys())[0][1]
	if symb == 'NA':
		symb = id
	if id in pcgid: type = 'PCG'
	elif id in lncid: type = 'lncRNA'
	else: type = 'other'
	if float(log2fc) >= 2 and float(padj) <= 0.05:
		realExpress = getAbsoluteExpression(id, matrix, samples)
		meanVar = calculateVarianceMean(realExpress, samples)
		vars.append(meanVar)
		degFilt_up.append([id,symb, type, log2fc, padj, realExpress, meanVar])
	elif float(log2fc) <= -2 and float(padj) <= 0.05:
		realExpress = getAbsoluteExpression(id, matrix, samples)
		meanVar = calculateVarianceMean(realExpress, samples)
		vars.append(meanVar)
		degFilt_down.append([id,symb, type, log2fc, padj, realExpress, meanVar])

varcutoff = np.quantile(vars, 0.75)
print varcutoff
degFilt_up = filter(lambda x: x[-1] <= varcutoff, degFilt_up)
degFilt_down = filter(lambda x: x[-1] <= varcutoff, degFilt_down)
sams = realExpress.keys()
print sams
if topN == 'all':
#	degFilt_up = mdl.comp(degFilt_up, lambda x: float(max(map(lambda y: float(y), x[-2].values()))), True)
#	degFilt_down = mdl.comp(degFilt_down, lambda x: float(max(map(lambda y: float(y), x[-2].values()))), False)
	degFilt_up = mdl.comp(degFilt_up, lambda x: float(x[-1]), False)
	degFilt_down = mdl.comp(degFilt_down, lambda x: float(x[-1]), True)
else:
	degFilt_up = mdl.comp(degFilt_up, lambda x: float(x[3]), True)[:int(topN)]
	degFilt_down = mdl.comp(degFilt_down, lambda x: float(x[3]), False)[:int(topN)]
outfile.write('\t'.join(sams)+'\n')
for deg in degFilt_up:
	outfile.write('\t'.join(deg[:5]))
	rexp = deg[5]
	for sam in sams:
		outfile.write('\t'+rexp[sam])
#	outfile2.write(str(deg[-1])+'\n')
	outfile.write('\n')

for deg in degFilt_down:
	outfile.write('\t'.join(deg[:5]))
	rexp = deg[5]
	for sam in sams:
		outfile.write('\t'+rexp[sam])
#	outfile2.write(str(deg[-1])+'\n')
	outfile.write('\n')
outfile.close()
#outfile2.close()
