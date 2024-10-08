import sys, os
import module as mdl

def get_type_mtx(filename,gtype,ctype):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[5:]
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; logfc = line[4]; exp = line[5:]; type = line[2]; genetype = line[3]
#		if gtype not in type: continue
		if 'PAR_Y' in id: continue
		if not geneDic.has_key((id,symb,type,genetype,logfc, ctype)):
			geneDic[(id,symb,type,genetype,logfc,ctype)] = dict()
		for x in xrange(len(exp)):
			sample = samples[x]; express = float(exp[x])
			geneDic[(id,symb,type,genetype,logfc, ctype)][sample] = express
	return geneDic,samples

def sumDuplicated_mtx(mtx):
	dup = []
	onedic = dict(); dupdic = dict()
	for gene in mtx.keys():
		geneid = gene[0]; logfc = gene[4]; celltype = gene[5]
		expdic = mtx[gene]
		if geneid in dup:
			if not dupdic.has_key(geneid):
				dupdic[geneid] =  [[],[]]
			dupdic[geneid][0].append(logfc)
			dupdic[geneid][1].append(celltype)

		dup.append(geneid)
	
	newdic = dict()
	for gene in mtx.keys():
		geneid = gene[0]
		expdic = mtx[gene]
		symb = gene[1]; type = gene[2]; genetype = gene[3]
		if geneid in dupdic.keys():
			logfc = ','.join(dupdic[geneid][0])
			celltype = ','.join(dupdic[geneid][1])
		else:
			logfc = gene[4]
			celltype = gene[5]
		newdic[(geneid,symb,type,genetype,logfc,celltype)] = expdic
	return newdic
		

	
def writeDegmtx(matdic, samples, filename):
	outfile = open(filename, 'w')
	outfile.write('ID\tsymbol\ttype\tgenetype\tlogFC\tdeg\t'+'\t'.join(samples))
	lines = ''
	for gene in matdic.keys():
		geneid = gene[0]; symbol = gene[1]; type = gene[2]; genetype = gene[3]
		expdic = matdic[gene]
		line = '\t'.join(gene)
		for sample in samples:
			line += '\t'+str(expdic[sample])
		lines += '\n'+line
	outfile.write(lines)
	outfile.close()


inputdir = sys.argv[1]
#degdir = sys.argv[2]
type = sys.argv[2]
outfile = sys.argv[3]

files = filter(lambda x: 'result.txt' in x, os.listdir(inputdir))

degmtx = dict()
for i in range(len(files)):
	file = files[i]
	celltype = file.split('.')[0]
	infile = inputdir + file
	mtx,samples = get_type_mtx(infile,type,celltype)
	print celltype, len(mtx.keys())
	degmtx.update(mtx)
srtMtx = sumDuplicated_mtx(degmtx)
writeDegmtx(srtMtx, samples, outfile)


