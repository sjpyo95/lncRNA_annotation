import sys, os
import module as mdl

def countElements(l):
	countdic = dict()
	elements = list(set(l))
	for i in elements:
		count = l.count(i)
		countdic[i] = count
	return countdic

cvfile = sys.argv[1]
infofile = sys.argv[2]
topN = int(sys.argv[3])
outfile = open(sys.argv[4], 'w')

infomtx = mdl.getinfoMatrix(infofile)

infile = open(cvfile, 'r')
lines = infile.readlines(); infile.close()
cols = lines[0]
samples = cols.split('\t')[1:-1]
cvs = list()

for i in range(1,len(lines)):
	line = lines[i].strip()
	tmp = line.split('\t')
	id = tmp[0]; exps = [float(x) for x in tmp[1:-1]]; cv = float(tmp[-1])
	hiexpIndex = exps.index(max(exps))
	hiexpSam = samples[hiexpIndex]
	cvs.append([id, cv, hiexpSam, exps[hiexpIndex]])
mdl.comp(cvs, lambda x: x[1], True)

topCVs = cvs[:topN]
sams = [x[2] for x in topCVs]
countSams = countElements(sams)

outfile.write('Cell-types\tcount\n')
for sam in countSams.keys():
	count = countSams[sam]
	outfile.write(sam + '\t' + str(count) + '\n')
outfile.close()
