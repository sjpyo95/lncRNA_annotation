import sys
import module as md

gtffile = sys.argv[1]
outfile = open(sys.argv[2],'w')
#matrix = sys.argv[3]
gtf = md.getGene(gtffile)
#mat = md.getMatrix(matrix)
genes = sum(gtf.values(), [])
print len(genes)
outfile.write('sense_exonic\tsense_intronic\tantisense\tbidirectional\tintergenic\tnovel_antisense\tnovel_bidirectional\tnovel_intergenic\n')
typecount = {'sense_exonic':0,'sense_intronic':0,'antisense':0,'bidirectional':0,'intergenic':0, 'novel_antisense':0, 'novel_bidirectional':0, 'novel_intergenic':0}
for i in xrange(len(genes)):
	genetype = genes[i].genetype()
	geneid = genes[i].geneid().split('.')[0]; symbol = genes[i].symb()
	attri = genes[i].attri()
#	if (geneid,symbol) in mat.keys():
#		if sum(map(lambda x: float(x), mat[(geneid,symbol)].values())) == 0: 
#			print geneid
#			continue
	if not 'novel' in attri:
		if 'sense_exonic' in genetype:
			typecount['sense_exonic'] += 1
		elif 'sense_intronic' in genetype:
			typecount['sense_intronic'] += 1
		elif 'antisense' in genetype:
			typecount['antisense'] += 1
		elif 'bidirectional' in genetype:
			typecount['bidirectional'] += 1
		elif 'lincRNA' in genetype:
			typecount['intergenic'] += 1
		else:
			print genetype
	elif 'novel' in attri:
		if 'antisense' in genetype:
			typecount['novel_antisense'] += 1
		elif 'bidirectional' in genetype:
			typecount['novel_bidirectional'] += 1
		elif 'lincRNA' in genetype:
			typecount['novel_intergenic'] += 1
		else:
			print genetype
outfile.write(str(typecount['sense_exonic'])+'\t'+str(typecount['sense_intronic'])+'\t'+str(typecount['antisense'])+'\t'+str(typecount['bidirectional'])+ '\t' + str(typecount['intergenic'])+'\t' +str(typecount['novel_antisense'])+'\t'+str(typecount['novel_bidirectional'])+'\t'+str(typecount['novel_intergenic']))
outfile.close()

