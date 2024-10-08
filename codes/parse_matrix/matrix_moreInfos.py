import sys, os
import module as mdl

matrixfile = sys.argv[1]
allfile = mdl.getGtf(sys.argv[2])
#pcgfile = mdl.getGtf2(sys.argv[3])
#lncfile = mdl.getGtf2(sys.argv[4])
novelfile = mdl.getGtf2(sys.argv[5])
outfilename = sys.argv[2]

#pcg = map(lambda x: x.geneid(), sum(mdl.getGtf2(pcgfile).values(), []))
#knownlnc = map(lambda x: x.geneid(), sum(mdl.getGtf2(lncfile).values(), []))
#novellnc = map(lambda x: x.geneid(), sum(mdl.getGtf2(novelfile).values(), []))

matrix, samples = mdl.getMatrix(matrixfile)

def geneinfo(gtf):
	genedic = dict()
	for chr in gtf.keys():
		genes = gtf[chr]
		for gene in genes:
			geneid = gene.geneid()
			if gene.elseinfo().has_key('pre_names'):
				prename = gene.elseinfo()['pre_names']
			else:
				prename = 'NA'
			trxs = gene.transcripts()
#			exons_nums = map(lambda x: x.exonNum(), trxs)
#			meanEx_num = float(sum(exons_nums))/len(exons_nums)
			genetype = gene.genetype()
			genedic[geneid] = genetype
	return genedic

pcg = geneinfo(pcgfile)
knownlnc = geneinfo(lncfile)
novellnc = geneinfo(novelfile)
newMatrix = dict()
for gene in matrix.keys():
	id = gene[0]; symbol = gene[1]
	exps = matrix[gene]
	if id in pcg.keys():
		infos = pcg[id]
		newMatrix[(id,symbol,'PCG', infos)] = exps
	elif id in knownlnc.keys():
		infos = knownlnc[id]
		newMatrix[(id,symbol,'Known_lncRNA', infos)] = exps
	elif id in novellnc.keys():
		infos = novellnc[id]
		newMatrix[(id,symbol,'Novel_lncRNA', infos)] = exps
	else:
		infos = 
expCount = dict()
for gene in newMatrix.keys():
	type = gene[2]
	if not expCount.has_key(type):
		expCount[type] = {}
	for sample in samples:
		exp = newMatrix[gene][sample]
		if float(exp) >= 1:
			if not expCount[type].has_key(sample):
				expCount[type][sample] = 0
			expCount[type][sample] += 1

outfile = open(outfilename, 'w')
outfile.write('ID\tsymbol\ttype\tgenetype\t'+'\t'.join(samples))
lines = ''
for gene in newMatrix.keys():
	id = gene[0]; symbol = gene[1]; type = gene[2]; genetype = gene[3]
	line = id+'\t'+symbol+'\t'+type+'\t'+genetype
	express = newMatrix[gene]
	for sample in samples:
		line += '\t'+ str(express[sample])
	lines += '\n'+line
outfile.write(lines)
outfile.close()

#total = dict()
#for type in expCount.keys():
#	for sample in expCount[type].keys():
#		count = expCount[type][sample]
#		if not total.has_key(sample):
#			total[sample] = 0.0
#		total[sample] += count

#outfile2 = open(outfilename[:-4]+'.stat.txt', 'w')
#outfile2.write('type\tsample\tcount\tpercentage\n')
#for type in expCount.keys():
#	for sample in expCount[type].keys():
#		count = expCount[type][sample]
#		totalcount = total[sample]
#		percentage = count/totalcount*100
#		outfile2.write(type +'\t'+sample+'\t'+str(count)+'\t'+str(percentage)+'\n')
#outfile2.close()
