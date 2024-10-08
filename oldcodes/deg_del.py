import sys
import module as mdl
degfile = sys.argv[1]
annofile = sys.argv[2]
chr = sys.argv[3]
outfile = open(sys.argv[4], 'w')
anno = mdl.getGtf(annofile)
chrAnno = list(set(map(lambda x: x.geneid().split('.')[0], anno[chr])))
def parseDEG(degfile):
	infile = open(degfile, 'r')
	lines = infile.readlines(); infile.close()
	degdic = dict()
	col = lines[0].strip()
	for i in xrange(1, len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		id = tmp[0]
		degdic[id] = line
	return col, degdic

col, deg = parseDEG(degfile)
outfile.write(col+'\n')
chrGenes = []
for id in deg.keys():
	if id in chrAnno: 
		chrGenes.append(id)
		continue
	outfile.write(deg[id]+'\n')
outfile.close()
print len(chrGenes)
