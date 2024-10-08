import sys
import module as mdl

#matrixfile = sys.argv[1]
gofile = sys.argv[1]
#terms = sys.argv[2]	#cytokine,chemokine,inflammatory,...
#terms = terms.split(',')
outputdir = sys.argv[2]
def getGOterms(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	go = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip()
		if len(line) == 0 : continue
		col = line.split('\t')
		source = col[0].split(':')[-1]; termName = col[1]; neglog10padj = col[4]
		genes = col[-1]
		if not go.has_key(source):
			go[source] = []
		go[source].append([termName, neglog10padj, genes])
#		print termName
	return go

goterms = getGOterms(gofile)
#print len(set(goterms))
#newMatrix = dict()
#for idsymb in matrix.keys():
#	id = idsymb[0]; symb = idsymb[1]
#	if symb in goterms:
#		newMatrix[idsymb] = matrix[idsymb]

#print len(newMatrix.keys())
#noGenes = filter(lambda x: x not in map(lambda x:x[1],newMatrix.keys()), goterms)
#print set(noGenes)
#mdl.writeMatrix(newMatrix, outfile)

for source in goterms.keys():
	go = goterms[source]
	sortGo = mdl.comp(go, lambda x: float(x[1]), True)[:20]
	outfile = open(outputdir + source+'.'+gofile, 'w')
#	print sortGo
	for line in sortGo:
		outfile.write(source+'\t'+'\t'.join(line)+'\n')
outfile.close()
