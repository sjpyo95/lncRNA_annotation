import sys, os
import module as mdl

def gethighcor(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	lncRNAs = lines[0].strip().split('\t')
	cordic = dict()
	for i in range(1,len(lines)):
		line = lines[i].strip().split('\t')
		pcgid = line[0]; cor_scores = line[1:]
		for l in range(len(lncRNAs)):
			lnc = lncRNAs[l]; cor_score = float(cor_scores[l])
			if not cordic.has_key(lnc):
				cordic[lnc] = []
			if abs(cor_score) >= 0.9:
				cordic[lnc].append(pcgid)
	return cordic

cordir = sys.argv[1]
mtxfile = sys.argv[2]
outputdir = sys.argv[3]

if not os.path.exists(outputdir): os.makedirs(outputdir)

mtx, samples = mdl.getinfoMatrix(mtxfile)

infiles = filter(lambda x: '.txt' in x, os.listdir(cordir))

for i in range(len(infiles)):
	corfile = cordir + infiles[i]
	celltype = infiles[i].split('.')[0]
	cordic = gethighcor(corfile)
	outfile = open(outputdir + celltype + '.high_cor_genelist.txt', 'w')
	lncs = cordic.keys()
	for lnc in lncs:
		lncsymb = filter(lambda x: x[0] == lnc, mtx.keys())[0][1]
		pcgs = cordic[lnc]
		outfile.write(lncsymb)
		for pcg in pcgs:
			pcgsymb = filter(lambda x: x[0] == pcg, mtx.keys())[0][1]
			outfile.write('\t'+pcgsymb)
		outfile.write('\n')
	outfile.close()


