import sys, os
import module as mdl

gtf1_file = sys.argv[1]
gtf2_file = sys.argv[2]

gtf1 = mdl.getGtf2(gtf1_file)
gtf2 = mdl.getGtf2(gtf2_file)

difgenes = []
for chr in gtf1.keys():
	genes1 = [x.geneid() for x in gtf1[chr]]
	genes2 = [x.geneid() for x in gtf2[chr]]
	difgenes += (filter(lambda x: x not in genes2, genes1))

print difgenes
