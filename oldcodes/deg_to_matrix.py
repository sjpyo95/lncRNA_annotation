import sys
import module as mdl
#gtf = mdl.getGtf(sys.argv[1])
deg = sys.argv[1]
matrix = open(sys.argv[2], 'r')
outfile = open(sys.argv[3], 'w')

def getDeg(degfile):
	degf = open(degfile, 'r')
	lines = degf.readlines(); degf.close()
	ids = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		ids[line[0]] = line[-2]
	return ids
genes = getDeg(deg)
print len(genes)
#genes = list(set([x.geneid() for x in sum(gtf.values(), [])]))
lines = matrix.readlines(); matrix.close()
samples = lines[0].strip().split('\t')
outfile.write(samples[0] + '\tsymbol\t' + '\t'.join(samples[1:])+'\n')
for i in xrange(1, len(lines)):
	line = lines[i].strip().split('\t')
	geneid = line[0]
	if geneid in genes.keys():
		outfile.write(geneid+'\t'+genes[geneid]+'\t'+'\t'.join(line[1:])+'\n')
outfile.close()

