import sys, os
import module as mdl

def getMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[4:]
	print len(samples)
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; exp = line[4:]; type = line[2]

		for x in xrange(len(exp)):

			sample = samples[x]; express = float(exp[x])
			if express >= 10:
				if not geneDic.has_key(sample):
					geneDic[sample] = []
				geneDic[sample].append([id,express,type])
	return geneDic

mtx_file = sys.argv[1]
out_file = sys.argv[2]

mtx = getMatrix(mtx_file)

outfile = open(out_file, 'w')
outfile.write('samples\ttype\texpression\n')

for sample in mtx:
	explist = mtx[sample]
	for expinfo in explist:
		id = expinfo[0]; express = expinfo[1]; type = expinfo[2]
		outfile.write(sample+'\t'+type+'\t'+str(express)+'\n')

outfile.close()



