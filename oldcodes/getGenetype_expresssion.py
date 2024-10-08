import sys, os
import module as mdl
exptable = open(sys.argv[1], 'r')
gtffile = sys.argv[2]
genetype = sys.argv[3]
outfile = open(sys.argv[4], 'w')
gtf = mdl.getGene(gtf)
sp_genes = filter(lambda x: genetype in x.genetype(), sum(gtf.values(),[]))

lines=  exptable.readlines(); exptable.close()
for i in xrange(len(lines)):
	line = lines[i].strip()
	if not line.startswith('ID'):
		id = line.split('\t')[0]
		if id in map(lambda x: x.geneid().split('.')[0], sp_genes):
			outfile.write(line+'\n')
	else:
		outfile.write(line+'\n')

