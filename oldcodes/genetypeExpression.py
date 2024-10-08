import sys, os
import module as mdl

matrixfile = sys.argv[1]
pcgfile = sys.argv[2]
knownlncfile = sys.argv[3]
novelfile = sys.argv[4]
outnamePrefix = sys.argv[5]
pcgs = [x.geneid() for x in sum(mdl.getGene(pcgfile).values(), [])]
knownlnc = [x.geneid() for x in sum(mdl.getGene(knownlncfile).values(), [])]
novel = [x.geneid() for x in sum(mdl.getGene(novelfile).values(), [])]

matrix = mdl.getMatrix(matrixfile)
outfile = open(outnamePrefix+'.expression.txt', 'w')
outfile.write('samples\ttypes\texpression\n')
samples = matrix.values()[0].keys()

pcgCountdic = {'1':0,'2':0,'3':0,'4':0,'>=5':0}
knownCountdic = {'1':0,'2':0,'3':0,'4':0,'>=5':0}
novelCountdic = {'1':0,'2':0,'3':0,'4':0,'>=5':0}
for idsymb in matrix.keys():
	geneid = idsymb[0]; symbol = idsymb[1]
	samplesEx = matrix[idsymb]
	for sample in samples:
		express = samplesEx[sample]
		if geneid in pcgs:
			outfile.write(sample+'\tpcg\t'+express+'\n')
		elif geneid in knownlnc:
			outfile.write(sample+'\tknown_lncRNA\t'+express+'\n')
		elif geneid in novel:
			outfile.write(sample+'\tnovel_lncRNA\t'+express+'\n')
	allexpress = samplesEx.values()
	if geneid in pcgs:
		pcgcount = sum(map(lambda x : x >= 0.5, allexpress))
		if pcgcount < 5:
			pcgCountdic[str(pcgcount)] += 1
		else: pcgCountdic['>=5'] += 1
	elif geneid in knownlnc:
		knowncount = sum(map(lambda x: x >= 0.5, allexpress))
		if knowncount < 5:
			knownCountdic[str(knowncount)] += 1
		else: knownCountdic['>=5'] += 1
	elif geneid in novel:
		novelcount = sum(map(lambda x: x >= 0.5, allexpress))
		if novelcount < 5:
			novelCountdic[str(novelcount)] += 1
		else: novelCountdic['>=5'] += 1


outfile2 = open(outnamePrefix+'.count.txt', 'w')
outfile2.write('type\t1\t2\t3\t4\t>=5\n')
outfile2.write('PCGs\t'+str(pcgCountdic['1'])+'\t'+str(pcgCountdic['2'])+'\t'+str(pcgCountdic['3'])+'\t'+str(pcgCountdic['4'])+'\t'+str(pcgCountdic['>=5'])+'\n')
outfile2.write('known_lncRNAs\t'+str(knownCountdic['1'])+'\t'+str(knownCountdic['2'])+'\t'+str(knownCountdic['3'])+'\t'+str(knownCountdic['4'])+'\t'+str(knownCountdic['>=5'])+'\n')
outfile2.write('novel_lncRNAs\t'+str(novelCountdic['1'])+'\t'+str(novelCountdic['2'])+'\t'+str(novelCountdic['3'])+'\t'+str(novelCountdic['4'])+'\t'+str(novelCountdic['>=5']))


outfile.close()
outfile2.close()



