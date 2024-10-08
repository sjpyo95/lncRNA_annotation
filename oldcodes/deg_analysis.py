import sys
import module as mdl
def comp(samplelist, stand, re):
	return sorted(samplelist, key=stand, reverse = re)

infile = sys.argv[1]
#matrix = open(sys.argv[2], 'r')
fc_cut = float(sys.argv[2])
pval_cut = float(sys.argv[3])
#outfile = open(sys.argv[4], 'w')
outputdir = sys.argv[4]
pcgs = sys.argv[5]
lnc1 = sys.argv[6]
lnc2 = sys.argv[7]
lines = open(infile, 'r').readlines()
deglist = []
pcgGtf = sum(mdl.getGene(pcgs).values(), [])
lncGtf1 = sum(mdl.getGene(lnc1).values(), [])
lncGtf2 = sum(mdl.getGene(lnc2).values(), [])
pcgid = map(lambda x: x.geneid(), pcgGtf)
lncid1 = map(lambda x: x.geneid().split('.')[0], lncGtf1)
lncid2 = map(lambda x: x.geneid(), lncGtf2)

for i in xrange(len(lines)):
	line = lines[i].strip()
	if i != 0:
		fc = float(line.split('\t')[4]); 
		logfc = float(line.split('\t')[5])
#		lfc = float(line.split('\t')[3])
		pval = line.split('\t')[6]#; padj = line.split('\t')[6]
		if pval == 'NA' : continue
		if abs(logfc) < fc_cut: continue
		elif float(pval) > pval_cut: continue
		genetype = 'Other'
		geneid = line.split('\t')[0]
		symbol = 'NA'; subtype = 'NA'
		if geneid in pcgid:
			genetype = 'PCG'
			gene = filter(lambda x: x.geneid() == geneid, pcgGtf)[0]
			symb = gene.symb(); subtype = gene.genetype()
		elif geneid.split('-')[-1] in lncid1:
			genetype = 'lncRNA'
			gene = filter(lambda x: x.geneid().split('.')[0] == geneid.split('-')[-1], lncGtf1)[0]
			symb = gene.symb(); subtype = gene.genetype()
		elif geneid in lncid2:
			genetype = 'lncRNA'
			gene = filter(lambda x: x.geneid() == geneid, lncGtf2)[0]
			symb = gene.symb(); subtype = gene.genetype()

#		for gene in genes:
#			print gene.genetype()
#			if geneid == gene.geneid().split('.')[0]:
#				genetype = gene.genetype()
#				break
		deglist.append([geneid, symb, fc, logfc, pval, genetype, subtype])
deglist = comp(deglist, lambda x: (x[2], x[4]), False)
outfile1 = open(outputdir + infile.split('/')[-1].split('.txt')[0] + '.pcg.txt','w')
outfile2 = open(outputdir + infile.split('/')[-1].split('.txt')[0] + '.lncRNA.txt','w')
outfile1.write('#|Foldchange| >= ' + str(fc_cut) + ' & P value <= ' + str(pval_cut) + '\n')
outfile1.write('id\tsymbol\tfoldchange\tlog2FoldChange\tpval\tgenetype\tsubtype\n')
outfile2.write('#|Foldchange| >= ' + str(fc_cut) + ' & P value <= ' + str(pval_cut) + '\n')
outfile2.write('id\tsymbol\tfoldchange\tlog2FoldChange\tpval\tgenetype\tsybtype\n')

for i in xrange(len(deglist)):
	deg = map(str, deglist[i])
	if deg[-2] == 'PCG':
		outfile1.write('\t'.join(deg)+'\n')
	elif deg[-2] == 'lncRNA':
		outfile2.write('\t'.join(deg)+'\n')
outfile1.close()
outfile2.close()
