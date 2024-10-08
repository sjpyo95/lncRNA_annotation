import sys, os

deg1 = sys.argv[1]
deg2 = sys.argv[2]
outfile = open(sys.argv[3], 'w')
def getDeseq(filename):
	infile = open(filename)
	lines = infile.readlines();infile.close()
	degdic = dict()
	col = lines[0]
	for i in xrange(1,len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('gene'):
			tline = line.split('\t')
#			id = tline[0]; log2fc = tline[3]; padj = tline[4]
			id = tline[0]; baseMean = tline[1]; log2fc = tline[2]; lfcse = tline[3]; stat = tline[4]
			pval = tline[5]; padj = tline[6]
			if float(padj) < 0.01 and abs(float(log2fc)) >= 4:
				degdic[id] = [log2fc, padj]
	return degdic
def getEdgeR(filename):
	infile = open(filename)
	lines = infile.readlines();infile.close()
	degdic = dict()
	col = lines[0]
	for i in xrange(1,len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('ID'):
			tline = line.split('\t')
			id = tline[0]; log2fc = tline[1]; padj = tline[-1]
#			id = line[0]; baseMean = line[1]; log2fc = line[2]; lfcse = line[3]; stat = line[4]
#			pval = line[5]; padj = line[6]
			if float(padj) < 0.01 and abs(float(log2fc)) >= 4:
				degdic[id] = [log2fc, padj]
	return degdic

def common_member(a,b):
	return list(set(a).intersection(b))

deg1 = getDeseq(deg1)
deg2 = getEdgeR(deg2)
print len(deg1)
print len(deg2)
deg1ids = deg1.keys()
deg2ids = deg2.keys()
outfile.write('gene\tdeseq_log2fc\tdeseq_padj\tedger_log2fc\tedger_padj\n')

comm_ids = common_member(deg1ids, deg2ids)
for id in comm_ids:
	if id == 'ID': continue
	deseq = deg1[id]
	edger = deg2[id]
	outfile.write(id+'\t'+'\t'.join(deseq)+'\t'+'\t'.join(edger)+'\n')
outfile.close()
