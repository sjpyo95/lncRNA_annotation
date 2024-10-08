import sys
import module as mdl

tablefile = sys.argv[1]
gtffile = sys.argv[2]
outfile = open(sys.argv[3], 'w')
table = mdl.getMatrix(tablefile)
gtf = sum(mdl.getGene(gtffile).values(), [])
typecount = dict()
for gene in gtf:
	geneid = gene.geneid().split('.')[0]
	symb = gene.symb()
	if symb == '': symb = 'NA'
	genetype = gene.genetype()
	if not (geneid,symb) in table.keys():
		continue
	express = table[(geneid,symb)].values()
	if len(filter(lambda x: float(x) >= 0.5, express)) > 0:
		if not typecount.has_key(genetype):
			typecount[genetype] = 0
		typecount[genetype] += 1

for type in typecount.keys():
	count = typecount[type]
	outfile.write(type+'\t'+str(count)+'\n')
outfile.close()
