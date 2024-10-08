import module as mdl
import sys, os

gtffile = sys.argv[1]
newfile = sys.argv[2]
reftypes = ['lncRNA','bidirectional_promoter_lncRNA','antisense','3prime_overlapping_ncRNA','lincRNA', 'macro_lncRNA', 'non_coding', 'processed_transcript', 'sense_intronic', 'sense_overlapping']
gtf = mdl.getGtf2(gtffile)
newgtf = dict()
for chr in gtf.keys():
	if not newgtf.has_key(chr):
		newgtf[chr] = []
	genes = gtf[chr]
	for gene in genes:
		genetype = gene.genetype()
		start = gene.start(); end = gene.end()
		genelen = abs(start-end+1)

		if genetype in reftypes:
			if genelen < 200:
				print gene.geneid()
				print genetype
			newgtf[chr].append(gene)

mdl.writeGtf3(newgtf, newfile)
