import sys
import module as mdl
gtf = sys.argv[1]
junction = sys.argv[2]
tssAnno = '/export/home/sjpyo/packages/cafe/src.1.0.1/data/anno/cage/human/hg19/hg19.cage_peak_coord_robust.bed'
cpsAnno = '/export/home/sjpyo/packages/cafe/src.1.0.1/data/anno/polya/human/hg19/hg19_all.15.sumCM.bed'

gtfdic = mdl.getGtf(gtf)
jundic = mdl.readingBed(junction)

for chr in gtfdic.keys():
	trxs = gtfdic[chr]
	juncs = jundic[chr]
	for trx in trxs:
		start = trx.start(); end = trx.end()
		for junc in juncs:
			junc
