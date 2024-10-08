import sys
import module as mdl
infile = sys.argv[1]
pcgs = sys.argv[2]
lncs = sys.argv[3]
lncs2 = sys.argv[4]
outfile = open(sys.argv[5], 'w')
cutoff = float(sys.argv[6]) #0.5
pcgid = list(set([x.geneid().split('.')[0] for x in sum(mdl.getGtf(pcgs).values(), [])]))
lncid = list(set([x.geneid().split('.')[0] for x in sum(mdl.getGtf(lncs).values(), [])]))
lncid2 = list(set([x.geneid().split('.')[0] for x in sum(mdl.getGtf(lncs2).values(), [])]))

def NA(l,index):
	for i in index:
		l[i] = 'NA'
	return l
		
def parseMatrix(infile):
	mat = open(infile, 'r')
	lines = mat.readlines(); mat.close()
	newlines = [lines[0].strip()]
	samples = lines[0].strip().split('\t')[2:]
	expdic = dict()
	for i in xrange(1, len(lines)):
		line = lines[i].strip()
		tline = line.split('\t')
		geneid = tline[0]
		express = tline[2:]
		nexpress = map(float, express)
		expdic[geneid] = nexpress
	samdic = dict()
	for s in xrange(len(samples)):
		sample = samples[s]
		if not samdic.has_key(sample):
#			samdic[sample] = {'PCGs':[], 'lncRNAs':[]}
			samdic[sample] = {'PCGs':[], 'known_lncRNAs':[], 'novel_lncRNAs':[]}
		for geneid in expdic.keys():
			express = expdic[geneid][s]
			if express >= cutoff:
				if geneid in pcgid:
					samdic[sample]['PCGs'].append(express)
				elif geneid in lncid:
#					samdic[sample]['lncRNAs'].append(express)
					samdic[sample]['known_lncRNAs'].append(express)
				elif geneid in lncid2:
					samdic[sample]['novel_lncRNAs'].append(express)
	print len(samdic.keys())
	return samdic

samdic = parseMatrix(infile)
#outfile = open(outprefix + '.'+str(cutoff) + '.txt', 'w')
outfile.write('samples\ttypes\ttpm\n')
for sample in samdic.keys():
	pcgE = samdic[sample]['PCGs']
#	lncE = samdic[sample]['lncRNAs']
	lncE = samdic[sample]['known_lncRNAs']
	lncE2 = samdic[sample]['novel_lncRNAs']

	for p in pcgE:
		outfile.write(sample+'\tPCG\t'+ str(p)+'\n')
	for l in lncE:
#		outfile.write(sample+'\tlncRNA\t'+ str(l) +'\n')	
		outfile.write(sample+'\tknown_lncRNA\t'+ str(l) +'\n')
	for l2 in lncE2:
		outfile.write(sample+'\tnovel_lncRNA\t'+ str(l2) +'\n')
#outfile.write('\n'.join(newlines))
outfile.close()
