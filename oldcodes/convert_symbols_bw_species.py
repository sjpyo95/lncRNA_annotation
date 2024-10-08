import sys, os

gmtfile = sys.argv[1]
bioMartfile = sys.argv[2]
outfile = open(sys.argv[3], 'w')
def getBiomart(filename):
	filein = open(filename, 'r')
	lines = filein.readlines(); filein.close()
	symbols = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip()
		genes = line.split(',')
		human = genes[0]; mouse = genes[1]
		symbols[human] = mouse
	return symbols

def getGmt(filename):
	filein = open(filename, 'r')
	lines = filein.readlines(); filein.close()
	gmt = dict()
	for i in xrange(len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		phenotype = tmp[0]; web = tmp[1]
		genes = tmp[2:]
		gmt[(phenotype,web)] = genes
	return gmt

gmt = getGmt(gmtfile)
mart = getBiomart(bioMartfile)

for p in gmt.keys():
	genes = gmt[p]
	phenotype = p[0]; web = p[1]
	outfile.write(phenotype+'\t'+web)
	for gene in genes:
		if gene in mart.keys():
			conv_gene = mart[gene]
			outfile.write('\t'+conv_gene)
	outfile.write('\n')
outfile.close()
