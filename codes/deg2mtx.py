import sys, os
import module as mdl

def parseRes(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	degdic = dict()
	for i in range(1,len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		id = tmp[0]; symbol = tmp[1]; type = tmp[2]; baseMean = tmp[3]
		fc = float(tmp[4]); padj = float(tmp[5])
		degdic[symbol] = [fc, padj]
		
	return degdic

def selmtx(res, mtx):
	newmtx = dict()
	for info in mtx.keys():
		id = info[0]; symbol = info[1]; type = info[2]; genetype = info[3]
		expression = mtx[info]
		if not symbol in res.keys(): continue
		fc = res[symbol][0]; padj = res[symbol][1]
		newinfo = info + (fc, padj)
		newmtx[newinfo] = expression
	return newmtx

def writedegMatrix(matdic, samples, filename):
	outfile = open(filename, 'w')
	outfile.write('ID\tsymbol\ttype\tgenetype\tlogFC\t'+'\t'.join(samples))
	lines = ''
	for gene in matdic.keys():
		geneid = gene[0]; symbol = gene[1]; type = gene[2]; genetype = gene[3]; fc = gene[4]
		expdic = matdic[gene]
		line = geneid+'\t'+symbol+'\t'+type+'\t'+genetype+'\t'+str(fc)
		for sample in samples:
			line += '\t'+str(expdic[sample])
		lines += '\n'+line
	outfile.write(lines)
	outfile.close()

inputdir = sys.argv[1]
mtxfile = sys.argv[2]
outputdir = sys.argv[3]
infiles = [inputdir + '/' + x for x in filter(lambda x: '.result.txt' in x, os.listdir(inputdir))]
mtx, samples = mdl.getinfoMatrix(mtxfile)

for infile in infiles:
	res = parseRes(infile)
	newmtx = selmtx(res, mtx)
	sample = infile.split('/')[-1].split('.result.')[0]
	outfile = outputdir + '/'+sample+'.DEG_info.TPM.txt'
	writedegMatrix(newmtx, samples, outfile)
