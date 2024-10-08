import sys

degfile = sys.argv[1]
outPre = sys.argv[2]
lncGtf = sys.argv[3]

def getGtf(infile):
	inf = open(infile, 'r')
	lines = inf.readlines(); inf.close()
	genes = dict()
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('#'):
			tline = line.split('\t')
			feature = tline[2]
			if feature != 'gene': continue
			tmp = tline[-1]
			geneid = tmp.split('gene_id')[1].split('"')[1]
			genetype = tmp.split('gene_type')[1].split('"')[1]
			if not geneid.startswith('ENSG'): 
				genes[geneid] = genetype
			else:
				geneN = tmp.split('gene_name')[1].split('"')[1]
				genes[geneid] = [geneN, genetype]
#			print geneid
#			print geneN
	return genes

knownlnc = getGtf(lncGtf)
#print knownlnc
def getDeg(degfile):
	filein = open(degfile, 'r')
	lines = filein.readlines(); filein.close()
	degdic = dict()
	for i in xrange(2, len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('#'):
			tline = line.split('\t')
			id = tline[0]; fc = tline[1]; logfc = tline[2]; pval = tline[3]; genetype = tline[-1]
			geneName = 'NA'; subtype = 'NA'
			if id in knownlnc.keys():
				if id.startswith('ENSG'):
					geneName = knownlnc[id][0]; subtype = knownlnc[id][1]
				else:
					subtype = knownlnc[id]
			if geneName.startswith('MIR'): continue
			degdic[id] = [fc, logfc, pval, genetype, geneName, subtype]
	return degdic

def comp(samplelist, stand, re):#re = True or False
	return sorted(samplelist, key=stand, reverse = re)


degs = getDeg(degfile)
#print degs
lncup = []; lncdown = []
pcgup = []; pcgdown = []
for id in degs.keys():
	val = degs[id]
	logfc = val[1]; pval = val[2]; genetype = val[3]; geneName = val[4];subtype = val[-1]
	if genetype == 'lncRNA':
		if float(logfc) < 0:
			lncup.append([id,logfc,pval,subtype, geneName])
		elif float(logfc)> 0:
			lncdown.append([id,logfc,pval,subtype, geneName])
	elif genetype == 'PCG':
		if float(logfc) < 0:
			pcgup.append([id,logfc,pval,subtype, geneName])
		elif float(logfc) > 0:
			pcgdown.append([id,logfc,pval,subtype, geneName])
lncup = comp(lncup, lambda x: (float(x[1])), False)
lncdown = comp(lncdown, lambda x: (float(x[1])),True)
lncout = open(outPre + '.lncRNA_degs.txt', 'w')
lncout.write('id\tlogFC\tpval\tgenetype\tgene_name\tregulation\n')
for val in lncup:
	lncout.write('\t'.join(val)+ '\tup\n')
for val in lncdown:
	lncout.write('\t'.join(val)+'\tdown\n')
lncout.close()


