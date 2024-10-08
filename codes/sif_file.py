import sys, os
import module as mdl

def parsefile(filename):
	infile= open(filename, 'r')
	lines = infile.readlines(); infile.close()
	cordic = dict()
	for i in range(1, len(lines)):
		line = lines[i].strip().split('\t')
		lncRNA = line[0]; cor = line[1]
		pcgs = line[2:]
		if not cordic.has_key(lncRNA):
			cordic[lncRNA] = dict()
		if not cordic[lncRNA].has_key('pos'):
			cordic[lncRNA]['pos'] = []
		if not cordic[lncRNA].has_key('neg'):
			cordic[lncRNA]['neg'] = []
		if cor == 'positive':
			cordic[lncRNA]['pos'].extend(pcgs)
		elif cor == 'negative':
			cordic[lncRNA]['neg'].extend(pcgs)
	return cordic


inputfile = sys.argv[1]
outSIFfile = open(sys.argv[2], 'w')

genes = parsefile(inputfile)

for lnc in genes.keys():
	pos_pcgs = genes[lnc]['pos']
	neg_pcgs = genes[lnc]['neg']
	if len(pos_pcgs) > 0:
		outSIFfile.write(lnc+'\tpositive\t'+'\t'.join(pos_pcgs)+'\n')
	elif len(neg_pcgs) > 0 :
		outSIFfile.write(lnc+'\tnegative\t'+'\t'.join(pos_pcgs)+'\n')
	
outSIFfile.close()
