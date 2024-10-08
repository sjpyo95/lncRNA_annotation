import sys, os
import module as mdl

inputdir = sys.argv[1]
type = sys.argv[2]	#count, fpkm, tpm, all
type = type.lower()
outputdir = sys.argv[3]
tag = sys.argv[4]
gtf = sys.argv[5]
anno = sum(mdl.getGene(gtf).values(), [])
ids = dict()
for gene in anno:
	if not 'ENSM' in gene.geneid():
		id = gene.geneid()
	else:
		id = gene.geneid().split('.')[0]
	symb = gene.symb()
	ids[id] = symb
if not os.path.exists(outputdir): os.makedirs(outputdir)
samples = filter(lambda x: '.txt' in x and not 'summary' in x, os.listdir(inputdir))
def getfpkm(count, totalcount, length):
	return float((10**9)*count)/(totalcount*length)

def getTotalfpkm(countDic):
	totalcount = sum([x[-1] for x in countDic.values()])
	totalfpkm = 0
	for geneid in countDic.keys():
		sample = countDic[geneid][0]
		length = countDic[geneid][1]
		count = countDic[geneid][2]
		totalfpkm += getfpkm(count, totalcount, length)
	return totalfpkm

def gettpm(count, totalfpkm, length):
	fpkm = getfpkm(count, totalcount, length)
	return float(fpkm*10**6)/(totalfpkm)

def parseFeatureCounts(infile):
	text = open(infile, 'r')
	lines = text.readlines(); text.close()
	countDic = dict()
#	name = infile.split('/')[-1].split('.')[0].split('_')[1]
	name = infile.split('/')[-1].split('.count')[0]
#	print name
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('#') and not line.startswith('Geneid'):
			tline = line.split('\t')
			geneid = tline[0]; chrom = tline[1]; start = tline[2]; end = tline[3] 
			sense = tline[4]; length = int(tline[5])
			count = int(tline[-1])
			countDic[geneid]=[name,length,count]
	return countDic

def makeMatrix(dic, outname):
	outfile = open(outname, 'w')
	names = sorted(dic.values()[0].keys())
#	print names
	outfile.write('ID\tsymbol\t'+'\t'.join(names)+'\n')
	for geneid in dic.keys():
		symb = 'NA'
#		if geneid.split('-')[-1].split('.')[0] in ids:
		if not'ENSM' in geneid:
			id = geneid
		else: id = geneid.split('.')[0]
		if id in ids:
#			symb = ids[geneid.split('-')[-1].split('.')[0]]
			symb = ids[id]
			if len(symb) == 0:
				symb = id
#		outfile.write(geneid.split('.')[0])
		if 'Mir' in symb: continue
#		if symb == 'NA': continue
		outfile.write(id+'\t'+symb)
		expressdic = dic[geneid]
#		print expressdic
		for name in names:
			express = expressdic[name]
			outfile.write('\t'+str(express))
		outfile.write('\n')
	outfile.close()

countMatrix = dict(); fpkmMatrix = dict(); tpmMatrix = dict()
for i in xrange(len(samples)):
	samplefile = inputdir + '/' + samples[i]
	countDic = parseFeatureCounts(samplefile)
	totalcount = sum([x[-1] for x in countDic.values()])
	totalfpkm = getTotalfpkm(countDic)
	print samples[i], str(totalcount)
	for geneid in countDic.keys():
		name = countDic[geneid][0]

		length = countDic[geneid][1]
		count = countDic[geneid][-1]
		if type == 'count' or type == 'all':
			if not countMatrix.has_key(geneid):
				countMatrix[geneid] = dict()
		 	countMatrix[geneid].update({name:count})
		if type == 'fpkm' or type == 'all':
			if not fpkmMatrix.has_key(geneid):
				fpkmMatrix[geneid] = dict()
			fpkmMatrix[geneid].update({name:getfpkm(count,totalcount,length)})
		if type == 'tpm' or type == 'all':
			if not tpmMatrix.has_key(geneid):
				tpmMatrix[geneid] = dict()
			tpmMatrix[geneid].update({name:gettpm(count,totalfpkm,length)})
if type == 'count' or type == 'all':
	outname = outputdir+tag+'.rCount.txt'
	makeMatrix(countMatrix, outname)
if type == 'fpkm' or type == 'all':
	outname = outputdir + tag + '.fpkm.txt'
	makeMatrix(fpkmMatrix, outname)
if type == 'tpm' or type == 'all':
	outname = outputdir+ tag + '.tpm.txt'
	makeMatrix(tpmMatrix, outname)
