import sys, os
import module as mdl
'''
	FPKM = (MR * 10^(6))/(GL*TMP)
	TPM = (FPKM * 10^(6))/sum(FPKM)
	
	MR: Number of Mapped Read to a gene
	GL: A gene's Gene Length
	TMP: Number of Total Mapped Read
'''

fcountdir = sys.argv[1] + '/'
outputdir = sys.argv[2] + '/'
if not os.path.exists(outputdir): os.makedirs(outputdir)
outfilePrefix = sys.argv[3]
gtffile = sys.argv[4]
gtf = mdl.getGene(gtffile)

def getIdSymb(gtf):
	idsymb = dict()
	for chr in gtf.keys():
		genes = gtf[chr]
		for gene in genes:
			id = gene.geneid(); symb = gene.symb()
			idsymb[id] = symb
	return idsymb

def getFcount(fcountfile, sample):
	fcount = open(fcountfile, 'r')
	lines = fcount.readlines(); fcount.close()
	countDic = dict(); totalcount = 0
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('#') and not line.startswith('Geneid'):
			tmp = line.split('\t')
#			print tmp
			geneid = tmp[0]; length = int(tmp[5]); count = int(tmp[-1])
			totalcount += count
			if not countDic.has_key(geneid):
				countDic[geneid] = {}
			countDic[geneid][sample] = [count, length]
	return countDic, totalcount

def getFPKM(samcount, tmp):
	samfpkm = dict(); totalfpkm = 0
	for geneid in samcount.keys():
		sample = samcount[geneid].keys()[0]
		mr = samcount[geneid][sample][0]; gl = samcount[geneid][sample][1]
		fpkm = (mr*1000000.0)/(gl*tmp)
		totalfpkm += fpkm
		if not samfpkm.has_key(geneid):
			samfpkm[geneid] = {}
		samfpkm[geneid][sample] = [fpkm, gl]
	return samfpkm, totalfpkm

def getTPM(samfpkm, totalfpkm):
	samtpm = dict(); totaltpm = 0
	for geneid in samfpkm.keys():
		sample = samfpkm[geneid].keys()[0]
		fpkm = samfpkm[geneid][sample][0]; gl = samfpkm[geneid][sample][1]
		tpm = (fpkm*1000000.0)/totalfpkm
		totaltpm += tpm
		if not samtpm.has_key(geneid):
			samtpm[geneid] = {}
		samtpm[geneid][sample] = [tpm, gl]
	return samtpm, totaltpm

def mergeSamples(alldic, samdic):
	for geneid in samdic.keys():
		sam = samdic[geneid]
		if not alldic.has_key(geneid):
			alldic[geneid] = {}
		alldic[geneid].update(sam)		
	return alldic

def writeMatrix(alldic, samples):
	lines = ''
	for geneid in alldic.keys():
		symbol = idsymb[geneid]
		if symbol == '':
			symbol = geneid
		line = geneid + '\t' + symbol
		countdic = alldic[geneid]
		for sample in samples:
			count = countdic[sample][0]
			line += '\t' + str(count)
		lines += '\n'+line
	return lines

fcounts = filter(lambda x: '.txt' in x and not '.summary' in x and not '.log' in x, os.listdir(fcountdir))
allcount_dic = dict(); allfpkm_dic = dict(); alltpm_dic = dict()
for i in xrange(len(fcounts)):
	sample = fcounts[i].split('.count')[0]
	fcountfile = fcountdir + fcounts[i]
	count_dic, totalcount = getFcount(fcountfile, sample)
	fpkm_dic, totalfpkm  = getFPKM(count_dic, totalcount)
	tpm_dic, totaltpm = getTPM(fpkm_dic, totalfpkm)

	print sample, totalcount
	print sample, totalfpkm
	print sample, totaltpm,'\n'
#####Merge Count dictionary#####
	allcount_dic = mergeSamples(allcount_dic, count_dic)
	allfpkm_dic = mergeSamples(allfpkm_dic, fpkm_dic)
	alltpm_dic = mergeSamples(alltpm_dic, tpm_dic)


#	for geneid in sam_dic.keys():
#		samCount = sam_dic[geneid]
#		if not count_dic.has_key(geneid):
#			count_dic[geneid] = {}
#		count_dic[geneid].update(samCount)
#print allcount_dic['ENSMUSG00000103006-4933417C20Rik']
#print allfpkm_dic['ENSMUSG00000103006-4933417C20Rik']
#print alltpm_dic['ENSMUSG00000103006-4933417C20Rik']
idsymb = getIdSymb(gtf)
samples = allcount_dic[allcount_dic.keys()[0]].keys()
firstline = 'ID\tsymbol\t'+'\t'.join(samples)



countfile = open(outputdir + outfilePrefix + '.rCount.txt', 'w')
countlines = writeMatrix(allcount_dic, samples)
countfile.write(firstline+countlines)
countfile.close()

#fpkmfile = open(outputdir + outfilePrefix + '.fpkm.txt', 'w')
#fpkmlines = writeMatrix(allfpkm_dic, samples)
#fpkmfile.write(firstline+fpkmlines)
#fpkmfile.close()

#tpmfile = open(outputdir + outfilePrefix + '.tpm.txt', 'w')
#tpmlines = writeMatrix(alltpm_dic, samples)
#tpmfile.write(firstline+tpmlines)
#tpmfile.close()

