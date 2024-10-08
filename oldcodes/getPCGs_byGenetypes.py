import sys
import module as mdl
import re
infile = sys.argv[1]
gtf = mdl.getGtf(infile)
outfile = sys.argv[2]
Genetypes = ['protein_coding', 'nonsense_mediated_decay', 'non_stop_decay', 'polymorphic_pseudogene']
immune = [re.compile('IG_(\S+)_gene'), re.compile('TR_(\S+)_gene'), re.compile('IG_(\S+)_pseudogene')]

def immuneTypeSelect(genetype, immune):
	for im in immune:
		if im.match(genetype):
			return True
	return False
newGtf = dict()
for chrom in gtf.keys():
	trxs = gtf[chrom]
	pcgs = filter(lambda x: x.genetype() in Genetypes, trxs)
	others = filter(lambda x: x.genetype() not in Genetypes, trxs)
	impcgs = filter(lambda x: immuneTypeSelect(x.genetype(), immune), others)
	allpcgs = pcgs+impcgs
	newGtf[chrom] = allpcgs
mdl.writeGtf(newGtf, outfile)
