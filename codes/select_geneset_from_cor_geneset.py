import sys, os
import module as mdl

def getGeneset(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	geneset = dict()
	for i in xrange(len(lines)):
		line = lines[i].strip().split('\t')
		lnc = line[0]; pcgs = line[1:]
		if len(pcgs) == 0: continue
		geneset[lnc] = pcgs

	return geneset
inputdir = sys.argv[1]
#in_file = sys.argv[1]

cutnum = int(sys.argv[2])

outputdir = sys.argv[3]#+'/'+in_file.split('/')[-1].split('.')[0]+'/'
if not os.path.exists(outputdir): os.makedirs(outputdir)

samples = filter(lambda x: '.txt' in x, os.listdir(inputdir))

for i in range(len(samples)):
	in_file = inputdir + samples[i]
	geneset = getGeneset(in_file)
	soutputdir = outputdir + '/'+samples[i].split('.')[0]+'/'
	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	j = 0
	for lnc in sorted(geneset, key=lambda lnc: len(geneset[lnc]), reverse=True):
		if j >= cutnum:break
		genes = geneset[lnc]
		outfile = open(soutputdir + lnc + '.cor_gene_list.txt', 'w')
		outfile.write('\n'.join(genes))
		outfile.close()
		j +=1
