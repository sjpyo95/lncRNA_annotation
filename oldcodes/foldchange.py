import sys
import module as mdl
import math
import numpy as np
from scipy import stats
infile = sys.argv[1]
metafile = sys.argv[2]
target = sys.argv[3]
background = sys.argv[4]
matrix = mdl.getMatrix(infile)
cutoff = sys.argv[5] # float or none
outputdir = sys.argv[6]
#def reverse_numeric(x,y):
#	return y - x
def comp(sampleList, stand):
	return sorted(sampleList, key=stand)#, cmp = reverse_numeric)
def log(x,n):
	return math.log(x,n)
def parseMeta(filein):
	infile = open(metafile, 'r')
	lines = infile.readlines(); infile.close()
	for i in xrange(1, len(lines)):
		line = lines[i].strip().split('\t')
		sample = line[0]; group = line[1]
def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
	from numpy import array, empty
	pvalues = array(pvalues)
	n = int(pvalues.shape[0])
	new_pvalues = empty(n)
	if correction_type == "Bonferroni":
		new_pvalues = n * pvalues
	elif correction_type == "Bonferroni-Holm":
		values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
		values.sort()
		for rank, vals in enumerate(values):
			pvalue, i = vals
			new_pvalues[i] = (n-rank) * pvalue
	elif correction_type == "Benjamini-Hochberg":
		values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
		values.sort()
		values.reverse()
		new_values = []
		for i, vals in enumerate(values):
			rank = n - i
			pvalue, index = vals
			new_values.append((float(n)/rank) * pvalue)

#		print new_values
		for i in xrange(0, int(n)-1):
			if new_values[i] < new_values[i+1]:
				new_values[i+1] = new_values[i]
		for i, vals in enumerate(values):
			pvalue, index = vals
			new_pvalues[index] = new_values[i]
	return new_pvalues

onlyGene = []
foldchange = []
pvals = []
for idsymb in matrix.keys():
	geneid = idsymb[0]; symb = idsymb[1]
	print geneid
	sampledic = matrix[idsymb]
	t_samples = filter(lambda x: target in x, sampledic.keys())
	b_samples = filter(lambda x: background in x, sampledic.keys())
#	print 'target: ' + ' '.join(t_samples)
#	print 'background: ' + ' '.join(b_samples)
	backExpress = map(lambda x: float(sampledic[x])+0.1, b_samples)
	targetExp = map(lambda x: float(sampledic[x])+0.1, t_samples)
#	if sum(backExpress)+sum(targetExp)/4 < 28.28: continue
#	otherExpress = map(lambda x: float(sampledic[x]), filter(lambda x: x != target, sampledic.keys()))
	print backExpress
	print targetExp
	ctrl_av = float(sum(backExpress))/len(backExpress)
	target_av = float(sum(targetExp))/len(targetExp)
	print ctrl_av
	print target_av
	fc = target_av/ctrl_av
	logfc = math.log(fc, 2)
	print fc
	print logfc
	if backExpress[0] == backExpress[1] and targetExp[0] == targetExp[1]:
		print 'no variance:'
		print targetExp[0], targetExp[1], backExpress[0], backExpress[1]
		continue
	if abs(logfc) <= 0:
		ttestR = [1,1]
#		print 'logFC = 0:'
		print targetExp[0], targetExp[1], backExpress[0], backExpress[1]
	else:
		ttestR = stats.ttest_ind(targetExp, backExpress, equal_var=False)
	foldchange.append([geneid, symb, logfc, ttestR[0], ttestR[-1]])
exit()
print len(foldchange)
padjs = correct_pvalues_for_multiple_testing(map(lambda x: x[-1], foldchange), correction_type = "Benjamini-Hochberg")
result = []
outfile = open(outputdir + target + '.' + background + '.fc.ttest.v12.txt', 'w')
outfile.write('ID\tsymbol\tlogFC\tTtest\tpval\tpadj\n')
for i in xrange(len(foldchange)):
	gene = foldchange[i]
	padj = padjs[i]
	if padj >= 1: padj = 1
	result.append(gene)
	outfile.write(gene[0]+'\t'+gene[1]+'\t'+str(gene[2])+'\t'+str(gene[3])+'\t'+str(gene[4])+'\t'+str(padj)+'\n')
#print len(filter(lambda x: x[-1] == 1, result))
	
	
	
#	if sum(otherExpress) == 0 and float(sampledic[target]) != 0:
#		
#		onlyGene.append([geneid, float(sampledic[target])])
#	else:
#		ctrl_av = float(sum(otherExpress))/len(otherExpress)
#		fchange = float(sampledic[target])/ctrl_av
#		logfc = math.log(fchange, 2)
#		foldchange.append([geneid, fchange, logfc])
#print onlyGene[0:10]

#sort_fc = comp(foldchange, lambda x: (float(x[1])))
#sort_logfc = comp(foldchange, lambda x: (float(x[2])))
#outfile = open(outputdir + target + '.' + background + '.fold_change.txt', 'w')
#if cutoff == 'none':
#	cutlogfc = filter(lambda x: x, foldchange)
#else:
#	cutlogfc = filter(lambda x: x[1] > float(cutoff), foldchange)
#cutlogfc = comp(cutlogfc, lambda x: (x[1]))
#outfile.write('ID\tfoldchange\tlog2fc\n')
#for i in xrange(len(cutlogfc)):
#	geneid = cutlogfc[i][0]; fc = cutlogfc[i][1]; logfc = cutlogfc[i][2]
#	outfile.write(geneid[0]+'\t'+str(fc)+'\t'+str(logfc)+'\n')
outfile.close()

