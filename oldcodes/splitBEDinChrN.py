import sys, os

filename = sys.argv[1]
outputdir = sys.argv[2]
if not os.path.exists(outputdir): os.makedirs(outputdir)
bedfile = open(filename, 'r')
lines = bedfile.readlines(); bedfile.close()
chrDic = dict()
for i in xrange(len(lines)):
	line = lines[i].strip()
	tmp =  line.split('\t')
	chr = tmp[0]
	if not chrDic.has_key(chr):
		chrDic[chr] = []
	chrDic[chr].append(line)

for chr in chrDic.keys():
	newfile = open(outputdir + filename.split('/')[-1][:-4]+'.'+chr+'.bed', 'w')
	lines = chrDic[chr]
	newfile.write('\n'.join(lines))
	newfile.close()

