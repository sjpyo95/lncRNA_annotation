import sys
import module as mdl

infile = sys.argv[1]
annofile = sys.argv[2]
term = sys.argv[3]
def parseGOterm(infile):
	filein = open(infile, 'r')
	lines = filein.readlines(); filein.close()
	godic = dict()
	for i in xrange(len(lines)):
		line = lines[i].strip()
		goterm = line.split('\t')[1]
		genes = line.split('\t')[-1]
		godic[goterm] = genes
	return godic
anno = sum(mdl.getGene(annofile).values(), [])
godic = parseGOterm(infile)
for goterm in godic.keys():
	genes = godic[goterm]
	if term in goterm:
		ids = map(lambda x: x.geneid().split('.')[0], filter(lambda x: x.symb().upper() in genes, anno))
		print goterm 
		print ids
		print genes

		
