import sys
import module as mdl

degfile = sys.argv[1]
gofile = sys.argv[2]
term = sysar.gv[3]
outfile = sys.argv[4]

def getGOterms(filename, term):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	gene = []
	for i in xrange(len(1,len(lines)):
		line=  lines[i].strip()
		if len(line) == 0: continue
		colm = line.split('\t'
		sym = colm[1]; cterm = colm[7]
		if cterm == term:
			genes.append(sym)
	return genes
def getfiltDEG(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	degs = dict()
	for i in xrange(1,len(lines)):
		tline = lines[i].strip()
		line = tline.split('\t')
		id = line[0]; symb = line[1]; logfc = line[3]; padj = line[4]
		degs[id] = [symb, logfc, padj]
goterms = getGOterms(gofile)
for id

