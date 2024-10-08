import sys, os
import module as mdl

gtffile = sys.argv[1]
celltype = sys.argv[2]
outfile = sys.argv[3]

gtf = mdl.getGtf2(gtffile)

for chr in gtf.keys():
	genes = gtf[chr]
	for gene in genes:
		gene.setElseinfo({'celltypes':celltype})
		trxs = gene.transcripts()
		for trx in trxs:
			trx.setElseinfo({

