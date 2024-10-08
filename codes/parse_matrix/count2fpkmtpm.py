import sys, os
import module as mdl

matrixfile = sys.argv[1]
fcountfile = sys.argv[2]
outputdir = sys.argv[3] + '/'
outfilePrefix = sys.argv[4]
outfilePrefix = matrixfile.split('/')[-1].split('.RC')[0]
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
#			if not geneid in matrix.keys(): continue
			length = int(tmp[-2])
			if not countDic.has_key(geneid):
				countDic[geneid] = {}
			countDic[geneid] = length
	return countDic

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
		if not geneDic.has_key((id,symbol)):
			geneDic[(id,symbol)] = dict()
		for x in xrange(len(exp)):
			sample = samples[x]; express = int(exp[x])
			geneDic[(id,symbol)][sample] = express
			if not totalcounts.has_key(sample):
				totalcounts[sample] = 0
			totalcounts[sample] += express
	return geneDic, totalcounts, samples

def getFPKM(matrix, geneLens, totalcounts):
	samfpkm = dict(); totalfpkm = dict()
	for idsym in matrix.keys():
		id = idsym[0]; symbol = idsym[1]
		samples = matrix[idsym].keys()
		gl = geneLens[id]
		for sample in samples:
			tmp = totalcounts[sample]
			mr = matrix[idsym][sample]
			fpkm = (mr*1000000.0)/(gl*tmp)
			if not totalfpkm.has_key(sample):
				totalfpkm[sample] = 0
			totalfpkm[sample] += fpkm
			if not samfpkm.has_key(idsym):
				samfpkm[idsym] = {}
			samfpkm[idsym][sample] = fpkm
	return samfpkm, totalfpkm

def getTPM(samfpkm, geneLens, totalfpkm):
	samtpm = dict()
	for idsym in samfpkm.keys():
		samples = samfpkm[idsym].keys()
		gl = geneLens[idsym[0]]
		for sample in samples:
			tf = totalfpkm[sample]
			fpkm = samfpkm[idsym][sample]
			tpm = (fpkm*1000000.0)/tf
			if not samtpm.has_key(idsym):
				samtpm[idsym] = {}
			samtpm[idsym][sample] = tpm
	return samtpm

matrix, totalcounts, samOrders = getMatrix(matrixfile)

print '\t'.join(samOrders)
#fcounts = filter(lambda x: '.txt' in x and not '.summary' in x and not '.log' in x, os.listdir(fcountdir))
geneLens = getGeneLen(fcountfile)
fpkms, totalfpkm  = getFPKM(matrix, geneLens, totalcounts)
tpms = getTPM(fpkms, geneLens, totalfpkm)

#samples = allcount_dic[allcount_dic.keys()[0]].keys()
firstline = 'ID\tsymbol\t'+'\t'.join(samOrders)

fpkmfile = outputdir + outfilePrefix + '.fpkm.txt'
fpkmlines = mdl.writeMatrix(fpkms, samOrders, fpkmfile)

tpmfile = outputdir + outfilePrefix + '.TPM.txt'
tpmlines = mdl.writeMatrix(tpms, samOrders, tpmfile)
