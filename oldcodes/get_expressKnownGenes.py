import sys, os
import module as mdl

anno = sys.argv[1]
reference = sys.argv[2]
outfile = sys.argv[3]

annodic = mdl.getGtf(anno)
refdic = mdl.getGtf(reference)
newDic = dict()
for chr in annodic.keys():
	trxs = annodic[chr]
	refs = refdic[chr]
	geneids = [x.geneid() for x in trxs]
	refids = [x.geneid() for x in refs]
	RefGenes = filter(lambda x: x.geneid() in geneids, refs)
	otherGenes = filter(lambda x: x.geneid() not in refids, trxs)
	expressedGenes = RefGenes + otherGenes
	newDic[chr] = expressedGenes
mdl.writeGtf(newDic, outfile)
