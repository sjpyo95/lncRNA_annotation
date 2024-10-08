import sys, os
import module as mdl

def getGMT(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	gmt = dict()
	for i in range(len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		term = tmp[0]; geneset = tmp[2:]
		gmt[term] = geneset
	return gmt

gmt_file = sys.argv[1]
mtx_file = sys.argv[2]
out_file = sys.argv[3]

gmt = getGMT(gmt_file)
mtx, samples = mdl.getinfoMatrix(mtx_file)
print len(gmt.keys()), 'PATHWAYS'

outfile = open(out_file, 'w')
outfile.write('Pathway\tID\tsymbol\ttype\tgenetype\t'+'\t'.join(samples))
genescheck = []
for term in gmt.keys():
	gene_symbols = gmt[term]
	sel_keys = filter(lambda x: x[1] in gene_symbols, mtx.keys())
	lines = ''
	pathway = '_'.join(term.split('_')[1:])
	for key in sel_keys:
		if key[1] in genescheck: continue
		genescheck.append(key[1])
		expdic = mtx[key]
		line = pathway+'\t'+'\t'.join(key)
		for sample in samples:
			line += '\t'+str(expdic[sample])
		lines += '\n'+line
	outfile.write(lines)
outfile.close()

