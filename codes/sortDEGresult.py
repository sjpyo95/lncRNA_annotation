import sys, os
import module as mdl

def parseResult(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	col = lines[0]
	degdic = dict()
	for i in range(1,len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		if 'NA' in tmp : continue
		id = tmp[0]; lgfc = float(tmp[2]); padj = float(tmp[-1]); baseMean = float(tmp[1])
		if lgfc < cfc: continue
		elif padj > cpadj : continue

		degdic[id] = [baseMean,lgfc, padj]
	return degdic

def rewriteResult(res, outfilename):
	outFile = open(outfilename, 'w')
	outFile.write('ids\tsymbols\ttype\tbaseMean\tlog2FoldChange\tpadj')
	for id in res.keys():
		vals = res[id]
		if len(vals) < 4: continue
		baseMean = vals[0];logfc = vals[1]; padj = vals[2]; symbol = vals[3]; type = vals[4]
		if 'ENSG' in id: id = id.split('.')[0]
#		if type == 'PCG' : continue
		outFile.write('\n' + id + '\t' + symbol + '\t' + type + '\t' + str(baseMean) + '\t' + str(logfc) + '\t' + str(padj))
	outFile.close()

resdir = sys.argv[1]
resfiles = [resdir + '/' + x for x in filter(lambda x: 'raw.txt' in x, os.listdir(resdir))]
mtxfile = sys.argv[2]
outputdir = sys.argv[3]
cfc = float(sys.argv[4])
cpadj = float(sys.argv[5])

mtx, samples = mdl.getinfoMatrix(mtxfile)

for resfile in resfiles:
	res = parseResult(resfile)
	for key in mtx.keys():
		id = key[0]; symb = key[1]; type = key[2]; genetype = [3]
		exps = mtx[key].values()

		if id not in res.keys(): continue
		elif genetype == 'intervening': 
			mtx.pop[id]
		else:
			res[id] += [symb, type]
	outfile = outputdir + resfile.split('/')[-1].split('.raw')[0] + '.result.txt'
	rewriteResult(res, outfile)
