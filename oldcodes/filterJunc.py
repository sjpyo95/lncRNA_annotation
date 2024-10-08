import sys

infile = sys.argv[1]
outPrefix = sys.argv[2]
junc = open(infile, 'r')
lines = junc.readlines(); junc.close()
newdic = dict()
for i in xrange(len(lines)):
	line = lines[i].strip()
	tline = line.split('\t')
	chr = tline[0]; start = tline[1]; end = tline[2]; hit = tline[3]; sense = tline[4]
	if int(hit) < 2: continue
	if not newdic.has_key(chr): newdic[chr] = []
	newdic[chr].append(chr+'\t'+start+'\t'+end+'\t'+hit+'\t'+sense+'\n')
for chr in newdic.keys():
	juncs = newdic[chr]
	outfile = open(outPrefix + chr + '.bed', 'w')
	for j in juncs:
		outfile.write(j)
	outfile.close()
	

