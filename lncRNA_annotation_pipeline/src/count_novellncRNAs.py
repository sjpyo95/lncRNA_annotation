import sys, os
import module as mdl

def countGeneTrx(infile):
#	genes = dict()
	genedic = dict()
	gt = mdl.getGtf2(infile)
	for chr in gt.keys():
		genes = gt[chr]
		for gene in genes:
			trxs = gene.transcripts()
#			if not genes.has_key(geneid):
#				genes[geneid] = []
#			genes[geneid].append(trx)
			if not genedic.has_key(gene):
				genedic[gene] = []
			genedic[gene] += trxs
	return genedic

inputdir = sys.argv[1]
outfile = open(sys.argv[2], 'w')

outfile.write('sample\tgene\ttranscript\n')

gtffiles = filter(lambda x: '.gtf' in x, os.listdir(inputdir))
for i in range(len(gtffiles)):
	gtf = gtffiles[i]
	sample = gtf.split('.')[0]
	gtffile = inputdir + gtf
	genes = countGeneTrx(gtffile)
	geneN = len(genes.keys())
	trxN = len(sum(genes.values(), []))
#	outfile = open(outputdir + sample + '.known-lncRNA.count.txt', 'w')
	outfile.write(sample+'\t'+str(geneN)+'\t'+str(trxN)+'\n')
#	print len(genes.keys())
#	print len(refs.keys())
outfile.close()
