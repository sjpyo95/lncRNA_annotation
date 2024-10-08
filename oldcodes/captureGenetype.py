import sys
import module as mdl
infile = sys.argv[1]
#matrix = open(sys.argv[2], 'r')
type = sys.argv[2]
outprefix = sys.argv[3]
#outfile = sys.argv[3]
gtf = mdl.getGtf(infile)
capgene = dict()
exgene = dict()
for chr in gtf.keys():
	trxs = gtf[chr]
#	trxs_type = filter(lambda x: type not in x.genetype() , trxs)
	for trx in trxs:
		if type not in trx.genetype(): 
			if not exgene.has_key(chr):
				exgene[chr] = []
			exgene[chr].append(trx)
		else:
			if not capgene.has_key(chr):
				capgene[chr] = []
			capgene[chr].append(trx)
#lines = matrix.readlines(); matrix.close()
#col = lines[0].strip()
#for i in xrange(1,len(lines)):
#	line = lines[i].strip()
#	id = line.split('\t')[0]
#	if id 

#outfile = outputdir + infile.split('/')[-1][:-4] + '_' + type + '.gtf'
mdl.writeGtf(capgene, outprefix+'.'+type+'.gtf')
mdl.writeGtf(exgene, outprefix+'.no_'+type+'.gtf')
#mdl.writeGtf(exgene, outfile)
