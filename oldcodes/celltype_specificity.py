import sys, os
import module as mdl
infile = sys.argv[1]
cutoff = float(sys.argv[2])
rangeN = (sys.argv[3])	# 1~5
rangeN = sorted(range(int(rangeN.split('~')[0]),int(rangeN.split('~')[1])+1))
outfile = open(sys.argv[4], 'w')

pcg = mdl.getGtf(sys.argv[5])
lnc = mdl.getGtf(sys.argv[6])
lnc2 = sys.argv[7]
#num = int(sys.argv[4])
pcgids = list(set([x.geneid().split('.')[0] for x in sum(pcg.values(), [])]))
lncids = list(set([x.geneid().split('.')[0] for x in sum(lnc.values(), [])]))
lnc2ids = [x.geneid() for x in sum(mdl.getGtf(lnc2).values(), [])]

def float2string(n):
	return str(n)

def parseMatrix(infile):
	mat = open(infile, 'r')
	lines = mat.readlines(); mat.close()
	samples = lines[0].strip().split('\t')[2:]
	result = dict()
	for i in range(len(samples)):
		sample = samples[i]
		for j in xrange(1,len(lines)):
			line = lines[j].strip()
			tline = line.split('\t')
			geneid = tline[0]; expression = tline[2:]
			if not result.has_key(geneid):
				result[geneid] = {}
			express = float(expression[i])
			result[geneid].update({sample:express})
#			print result
	return result
#lnctype = getGenetype(lnc)
expressDic = parseMatrix(infile)
gcount = 0
#outfile.write('expressN\tPCGs\tGENCODE\tnovel\n')
#outfile.write('n_type\tPCGs\tlncRNAs\n')
countDic = dict()
for geneid in expressDic.keys():
	sampleDic = expressDic[geneid]
	satSamples = filter(lambda x: sampleDic[x] >= cutoff, sampleDic.keys())
	satCount = len(satSamples)
	if not countDic.has_key('PCGs'):
		countDic['PCGs'] = dict.fromkeys(rangeN, 0)
	if not countDic.has_key('lncRNAs'):
		countDic['lncRNAs'] = dict.fromkeys(rangeN, 0)
#	if not countDic.has_key('GENCODE'):
#		countDic['GENCODE'] = dict.fromkeys(rangeN, 0)
	if not countDic.has_key('novel'):
		countDic['novel'] = dict.fromkeys(rangeN, 0)
	if satCount == 0 : continue
	if satCount < rangeN[-1]:
	 	if geneid in pcgids:
			countDic['PCGs'][satCount] += 1
		elif geneid in lncids:
			countDic['lncRNAs'][satCount] += 1
#			countDic['GENCODE'][satCount] += 1
		elif geneid in lnc2ids:
			countDic['novel'][satCount] += 1
	elif satCount >= rangeN[-1]:
	 	if geneid in pcgids:
	 		countDic['PCGs'][rangeN[-1]] += 1
	 	elif geneid in lncids:
	 		countDic['lncRNAs'][rangeN[-1]] += 1
#			countDic['GENCODE'][rangeN[-1]] += 1
		elif geneid in lnc2ids:
			countDic['novel'][rangeN[-1]] += 1

#print countDic
totalpcg = sum(countDic['PCGs'].values())
totallnc = sum(countDic['lncRNAs'].values())
totalnovel = sum(countDic['novel'].values())
totalgenes = totalpcg+totallnc+totalnovel
print totalpcg
print totallnc
print totalnovel
#totallnc = sum(countDic['GENCODE'].values()) + sum(countDic['novel'].values())
for type in countDic.keys():
	geneN = countDic[type]
	for n in geneN.keys():
		count = geneN[n]
		if n == rangeN[-1]:
			if type == 'PCGs':
				outfile.write(type+ '\t>' + str(n)+'_type\t' + str(float(count)/totalpcg)+ '\n')
			elif type == 'novel':
				outfile.write(type+'\t>' + str(n)+'_type\t' + str(float(count)/totallnc)+ '\n')
			else:
				outfile.write(type+ '\t>' + str(n)+'_type\t' + str(float(count)/totalnovel)+ '\n')
		else:
			if type == 'PCGs':
				outfile.write(type+ '\t' + str(n)+'_type\t' + str(float(count)/totalpcg)+ '\n')
			elif type == 'novel':
				outfile.write(type+'\t' + str(n)+'_type\t' + str(float(count)/totallnc)+ '\n')
			elif type == 'lncRNAs':
				outfile.write(type+ '\t' + str(n)+'_type\t' + str(float(count)/totalnovel)+ '\n')
outfile.close()
