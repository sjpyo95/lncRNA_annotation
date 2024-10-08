import sys

deg1 = sys.argv[1]
deg2 = sys.argv[2]
outputdir = sys.argv[3]

def getDegGenes(degfile):
	deg = open(degfile,'r')
	lines= deg.readlines(); deg.close()
	degdic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; gtype = line[-3]
		name = line[-2]
		reg = line[-1]
		logfc = line[1]
		degdic[id] = [gtype,reg, name, logfc]
	return degdic
def countType(degs, ids):
	inter = {'antisense': 0, 'lincRNA':0,'bidirectional':0, 'sense_intronic':0, 'etc':0}
	for id in ids:
		gtype = degs[id][0]
		reg1 = degs[id][1]
#		reg2 = degs[id][1]
#		if reg1 != reg2: print id
		if 'antisense' in gtype:
			inter['antisense'] += 1
		elif 'bidirectional_promoter_lncRNA' in gtype:
			inter['bidirectional'] += 1
		elif 'lincRNA' in gtype:
			inter['lincRNA'] += 1
		elif 'sense_intronic' in gtype:
			inter['sense_intronic'] += 1
		else:
			inter['etc'] += 1
	return inter
deg1Name = deg1.split('/')[-1]
deg2Name = deg2.split('/')[-1]
deg1genes = getDegGenes(deg1)
deg2genes = getDegGenes(deg2)
ids1 = deg1genes.keys()
ids2 = deg2genes.keys()
typecount = {}
interids = list(set(deg1genes.keys()).intersection(set(deg2genes.keys())))
unids = list(set(deg1genes.keys()).union(set(deg2genes.keys())))
interN = len(interids); uniN = len(unids)

typecount[deg1Name] = countType(deg1genes, ids1)
typecount[deg2Name] = countType(deg2genes, ids2)
typecount[deg1Name+ ' n ' + deg2Name] = countType(deg1genes, interids)

print 'sample\tantisense\tlincRNA\tbidirectional\tsense_intronic\tETC.'
for sample in typecount.keys():
	counts = typecount[sample]
	print sample+ '\t' +str(counts['antisense'])+'\t'+str(counts['lincRNA'])+'\t'+str(counts['bidirectional'])+'\t'+str(counts['sense_intronic'])+'\t'+str(counts['etc'])
outfile1 = open(outputdir + 'intersection.genelist.txt', 'w')
outfile2 = open(outputdir + 'union.genelist.txt', 'w')

outfile1.write('ID\ttype\tName\tlogFC\n')
outfile2.write('ID\ttype\tName\tlogFC\n')
for id in interids:
	if deg1genes[id][1] != 'up':
		continue
	elif deg1genes[id][1] != deg2genes[id][1]:
		print id, deg1genes[id], deg2genes[id]
		continue
	outfile1.write(id + '\t' +  deg1genes[id][0]+ '\t' + deg1genes[id][2]+ '\t' + deg1genes[id][-1]+'\n')
outfile1.close()
for id in unids:
	if id in ids1:
		if deg1genes[id][1] != 'up':	
			continue
		outfile2.write(id+'\t' +  deg1genes[id][0]+ '\t' + deg1genes[id][2]+ '\t' + deg1genes[id][-1] + '\n')
	elif id in ids2:
		if deg2genes[id][1] != 'up':
			continue
		outfile2.write(id+'\t' +  deg2genes[id][0]+ '\t' + deg2genes[id][2]+ '\t' + deg2genes[id][-1] + '\n')
outfile2.close()
