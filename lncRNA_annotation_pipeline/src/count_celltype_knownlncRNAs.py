import sys, os
import module as mdl

def readCol(infile):
	filein = open(infile, 'r')
	lines = filein.readlines(); filein.close()
	groupdic = dict()
	for i in range(1,len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		sample = tmp[0].split('_rep')[0]; symbol = tmp[1]
		group1 = tmp[3]
		if not groupdic.has_key(group1):
			groupdic[group1] = []
		if not sample in groupdic[group1]:
			groupdic[group1].append(sample)
	return groupdic

gtf = mdl.getGtf2(sys.argv[1])
coldata = sys.argv[2]
outfile = open(sys.argv[3], 'w')
groupdic = readCol(coldata)
print groupdic

countdic = dict()
for chr in gtf.keys():
	genes = gtf[chr]
	for gene in genes:
		celltypes = gene.elseinfo()['celltypes'].split(',')
		refnames = gene.elseinfo()['REF'].split(',')
		for group in groupdic.keys():
			samples = groupdic[group]
			for sample in samples:
					if sample in celltypes:
						if not countdic.has_key(group):
							countdic[group] = 0
						countdic[group] += 1
						break
outfile.write('sample\tcount\n')
for group in countdic.keys():
	count = countdic[group]
	outfile.write(group + '\t' + str(count) + '\n')
outfile.close()
