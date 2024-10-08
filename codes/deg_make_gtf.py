import sys, os
import module as mdl

def getdegMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[6:]
	genes = []
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]
		genes.append(id)
	return genes

inputdir = sys.argv[1]
files = filter(lambda x: '.txt' in x, os.listdir(inputdir))

gtffile = sys.argv[2]
gtf = mdl.getGtf2(gtffile)
outputdir = sys.argv[3]

if not os.path.exists(outputdir) : os.makedirs(outputdir)
for i in range(len(files)):
	file = inputdir + files[i]
	genetype = files[i].split('_degs')[0]
	degenes = getdegMatrix(file)
	degGtf = dict()
	for chr in gtf.keys():
		genes = gtf[chr]
		degGtfgenes = filter(lambda x: x.geneid() in degenes, genes)
		degGtf[chr] = degGtfgenes
	outfile = outputdir + genetype + '.degs.gtf'
	mdl.writeGtf3(degGtf, outfile)

