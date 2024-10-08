import sys, os
import module as mdl

inputdir = sys.argv[1]
#outfile = open(sys.argv[2], 'w')
def countGeneTrx(infile):
#	genes = dict()
	refs = dict()
	gt = mdl.getGtf(infile)
	for chr in gt.keys():
		trxs = gt[chr]
		for trx in trxs:
			geneid = trx.geneid()
			elseinfo = trx.elseinfo()
			ref = elseinfo['REF']
#			if not genes.has_key(geneid):
#				genes[geneid] = []
#			genes[geneid].append(trx)
			if not refs.has_key(ref):
				refs[ref] = []
			refs[ref].append(trx)
#	return genes,refs
	return refs

#outfile.write('sample\tgene\ttranscript\n')

gtffiles = filter(lambda x: '.gtf' in x, os.listdir(inputdir))
for i in range(len(gtffiles)):
	gtf = gtffiles[i]
	sample = gtf.split('.')[0]
	gtffile = inputdir + gtf
	refs = countGeneTrx(gtffile)
	geneN = len(refs.keys())
	trxN = len(sum(refs.values(), []))
#	outfile = open(outputdir + sample + '.known-lncRNA.count.txt', 'w')
#	outfile.write(sample+'\t'+str(geneN)+'\t'+str(trxN)+'\n')
#	print len(genes.keys())
	print len(refs.keys())
#outfile.close()
