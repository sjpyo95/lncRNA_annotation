import sys, os
import module as mdl

matrixfile = sys.argv[1]
fcountfile = sys.argv[2]
metafile = sys.argv[3]
outputdir = sys.argv[4] + '/'
outfilePrefix = sys.argv[5]
gtffile = sys.argv[6]
gtf = mdl.getGene(gtffile)

def getMeta(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	samOrders = []
	for i in range(1, len(lines)):
		line = lines[i].strip().split('\t')
		sample = line[0]
		samOrders.append(sample)
	return samOrders

def getIdSymb(gtf):
	idsymb = dict()
	for chr in gtf.keys():
		genes = gtf[chr]
		for gene in genes:
			id = gene.geneid(); symb = gene.symb()
			idsymb[id] = symb
	return idsymb

def getGeneLen(fcountfile):
	fcount = open(fcountfile, 'r')
	lines = fcount.readlines(); fcount.close()
	countDic = dict(); totalcount = 0
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('#') and not line.startswith('Geneid'):
			tmp = line.split('\t')
#			print tmp
			geneid = tmp[0]; 
			if not geneid in matrix.keys(): continue
			length = int(tmp[5])
			if not countDic.has_key(geneid):
				countDic[geneid] = {}
			countDic[geneid] = length
	return countDic


def getRcount(matrix, sample, geneLens):
	totalcount = 0
	for id in matrix.keys():
		geneLen = geneLens[id]
		count = int(matrix[id][sample])
		totalcount += count
		if not countDic.has_key(geneid):
			countDic[geneid] = {}
		countDic[geneid][sample] = [count, length]
	return countDic, totalcount

def getMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[2:]
	totalcounts = dict()
#	print samples
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symbol = line[1]; exp = line[2:]
#		print exp
		if not geneDic.has_key(id):
			geneDic[id] = dict()
		for x in xrange(len(exp)):
			sample = samples[x]; express = int(exp[x])
			geneDic[id][sample] = express
			if not totalcounts.has_key(sample):
				totalcounts[sample] = 0
			totalcounts[sample] += express
	return geneDic, totalcounts

def getFPKM(matrix, geneLens, totalcounts):
	samfpkm = dict(); totalfpkm = dict()
	for geneid in matrix.keys():
		samples = matrix[geneid].keys()
		gl = geneLens[geneid]
		for sample in samples:
			tmp = totalcounts[sample]
			mr = matrix[geneid][sample]
			fpkm = (mr*1000000.0)/(gl*tmp)
			if not totalfpkm.has_key(sample):
				totalfpkm[sample] = 0
			totalfpkm[sample] += fpkm
			if not samfpkm.has_key(geneid):
				samfpkm[geneid] = {}
			samfpkm[geneid][sample] = fpkm
	return samfpkm, totalfpkm

def getTPM(samfpkm, geneLens, totalfpkm):
	samtpm = dict()
	for geneid in samfpkm.keys():
		samples = samfpkm[geneid].keys()
		gl = geneLens[geneid]
		for sample in samples:
			tf = totalfpkm[sample]
			fpkm = samfpkm[geneid][sample]
			tpm = (fpkm*1000000.0)/tf
			if not samtpm.has_key(geneid):
				samtpm[geneid] = {}
			samtpm[geneid][sample] = tpm
	return samtpm

def mergeSamples(alldic, samdic):
	for geneid in samdic.keys():
		sam = samdic[geneid]
		if not alldic.has_key(geneid):
			alldic[geneid] = {}
		alldic[geneid].update(sam)		
	return alldic

def writeMatrix(expdic, idsymb, samOrders):
	lines = ''
	for geneid in expdic.keys():
		symbol = idsymb[geneid]
		if symbol == '':
			symbol = geneid
		line = geneid + '\t' + symbol
		dic = expdic[geneid]
		for sample in samOrders:
			exp = expdic[geneid][sample]
			line += '\t' + str(exp)
		lines += '\n'+line
	return lines
idsymb = getIdSymb(gtf)
matrix, totalcounts = getMatrix(matrixfile)
samOrders = getMeta(metafile)
 
print '\t'.join(samOrders)
#fcounts = filter(lambda x: '.txt' in x and not '.summary' in x and not '.log' in x, os.listdir(fcountdir))
geneLens = getGeneLen(fcountfile)

fpkms, totalfpkm  = getFPKM(matrix, geneLens, totalcounts)
tpms = getTPM(fpkms, geneLens, totalfpkm)

#samples = allcount_dic[allcount_dic.keys()[0]].keys()
firstline = 'ID\tsymbol\t'+'\t'.join(samOrders)

fpkmfile = open(outputdir + outfilePrefix + '.fpkm.txt', 'w')
fpkmlines = writeMatrix(fpkms, idsymb, samOrders)
fpkmfile.write(firstline+fpkmlines)
fpkmfile.close()

tpmfile = open(outputdir + outfilePrefix + '.tpm.txt', 'w')
tpmlines = writeMatrix(tpms, idsymb, samOrders)
tpmfile.write(firstline+tpmlines)
tpmfile.close()


		

