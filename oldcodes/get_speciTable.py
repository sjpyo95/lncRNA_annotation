import sys
import module as mdl

tablefile = sys.argv[1]
gtffile = sys.argv[2]
outfile = open(sys.argv[3],'w')
gtf = sum(mdl.getGene(gtffile).values(), [])
ids = []
for gene in gtf:
	geneid = gene.geneid()
	if not 'ENSMUS' in geneid:
		ids.append(geneid)
	else:
		ids.append(geneid.split('.')[0])
table = open(tablefile, 'r')
lines = table.readlines(); table.close()
for i in xrange(len(lines)):
	line = lines[i].strip()
	if line.startswith('ID'):
		outfile.write(line+'\n')
	geneid = line.split('\t')[0]
	if geneid in ids:
		outfile.write(line+'\n')
outfile.close()
