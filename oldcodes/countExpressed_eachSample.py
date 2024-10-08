import sys, os
import module as mdl
matrix = sys.argv[1]
cutoff = float(sys.argv[2])
outTable = open(sys.argv[3], 'w')
pcg = sys.argv[4]
lnc = sys.argv[5]
lnc2 = sys.argv[6]
#num = int(sys.argv[4])
pcgids = list(set([x.geneid().split('.')[0] for x in sum(mdl.getGtf(pcg).values(), [])]))
lncids = list(set([x.geneid().split('.')[0] for x in sum(mdl.getGtf(lnc).values(), [])]))
lnc2ids = [x.geneid() for x in sum(mdl.getGtf(lnc2).values(), [])]
def parseMatrix(infile):
	mat = open(infile, 'r')
	lines = mat.readlines(); mat.close()
	expressDic = dict()
	for i in xrange(1, len(lines)):
		line = lines[i].strip()
		tline = line.split('\t')
		geneid = tline[0]
		express = tline[2:]
		expressDic[geneid] = express
	samples = lines[0].strip().split('\t')[2:]
	result = dict()
	for i in range(len(samples)):
		sample = samples[i]
		for geneid in expressDic.keys():
			express = float(expressDic[geneid][i])
			if not result.has_key(sample):
				result[sample] = []
			result[sample].append([geneid, express])
	return result

expressDic = parseMatrix(matrix)
csegCount = dict()
testN = 0
totalpcg = 0; totallnc = 0; totallnc2 = 0

totalGenes = sum(expressDic.values(), [])
overGenes = []
for i in xrange(len(totalGenes)):
	geneid = totalGenes[i][0]; express = totalGenes[i][1]
	if geneid in overGenes: continue
	if geneid in lncids:
		testN += 1
	if express >= cutoff:
		overGenes.append(geneid)
		if geneid in pcgids:
			totalpcg += 1
		elif geneid in lncids:
		 	totallnc += 1
#		else:
#		 	totallnc2 += 1
#print testN
#exit()
for sample in expressDic.keys():
	expression = expressDic[sample]
	if not csegCount.has_key(sample):
#		csegCount[sample] = {'PCGs':0,'GENCODE':0, 'novel':0}
		csegCount[sample] = {'PCGs':0,'lncRNAs':0, 'novel':0}
	for gene in expression:
		geneid = gene[0]; express = gene[1]
		if express > cutoff:
			if geneid in lncids:
				csegCount[sample]['lncRNAs'] += 1
#				csegCount[sample]['GENCODE'] += 1
			elif geneid in pcgids:
				csegCount[sample]['PCGs'] += 1
			else:
				csegCount[sample]['novel'] += 1
#outTable.write('samples\tPCGs\tGENCODE\tnovel\n')
outTable.write('samples\tPCGs\tlncRNAs\tlncRNAs/totalPCGs\tnovel/totalPCGs\tlncRNAs/PCGs10K\tnovel/PCGs10K\n')
for sample in csegCount.keys():
	pcgN = csegCount[sample]['PCGs']
	lncN = csegCount[sample]['lncRNAs']
	lnc_tpcg = float(lncN)/totalpcg
	lnc_10k = float(10000*lncN)/pcgN
#	lncN = csegCount[sample]['GENCODE']
	lnc2N = csegCount[sample]['novel']
	lnc2_tpcg = float(lnc2N)/totalpcg
	lnc2_10k = float(1000*lncN)/pcgN
#	outTable.write(sample+'\t'+ str(pcgN)+'\t'+str(lncN)+'\t'+str(lnc2N)+'\n')
#	pcgN = csegCount[sample]['PCGs']
#	lncN = csegCount[sample]['lncRNAs']
	outTable.write(sample+'\t'+str(pcgN)+'\t'+str(lncN) + '\t' + str(lnc_tpcg) + '\t' + str(lnc2_tpcg) + '\t' + str(lnc_10k) + '\t' + str(lnc2_10k) + '\n')
#outTable.write('Total\t'+ str(totalpcg)+'\t'+str(totallnc))
#outTable.write('Total\t'+ str(totalpcg)+'\t'+str(totallnc) + '\t' + str(totallnc2))
outTable.close()
