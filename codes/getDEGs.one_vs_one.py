import sys, os
import module as mdl

def getDEG(infile):
	filein = open(infile, 'r')
	sample = infile.split('/')[-1].split('.deseq2')[0]
	lines = filein.readlines(); filein.close()
	degs = []
	for i in range(1, len(lines)):
		line = lines[i].strip()
		col = line.split('\t')
		id = col[0]; lgfc = float(col[2]); padj = float(col[-1])
		if lgfc >= 1 and padj <= 0.05:
			 degs.append(id)
	return degs

def getMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[4:]
	print len(samples)
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; exp = line[4:]; type = line[2]; genetype = line[3]
		if 'PAR_Y' in id: continue
		if not geneDic.has_key((id,symb,type,genetype)):
			geneDic[(id,symb,type,genetype)] = dict()
		for x in xrange(len(exp)):
			sample = samples[x]; express = float(exp[x])
			geneDic[(id,symb,type,genetype)][sample] = express
	return geneDic, samples

#def getgroup(filename, groupname):
#	filein = open(filename, 'r')
#	lines = filein.readlines(); filein.close()
#	groupnames = lines[0].strip().split('\t')[3:]
#	groupindex = groupnames.index(groupname)
#	groupdic = dict()
#
#	for i in range(1, len(lines)):
#		line = lines[i].strip()
#		tmp = line.split('\t')
#		symbol = tmp[1]; groups = tmp[3:]
#		group = groups[groupindex]
#		if not groupdic.has_key(group):
#			groupdic[group] = []
#		groupdic[group].append(symbol)
#	return groupdic
		
degdir = sys.argv[1]
mtx_file = sys.argv[2]
#groupname = sys.argv[3]
#group_file = sys.argv[4]
outputdir = sys.argv[3]
if not os.path.exists(outputdir): os.makedirs(outputdir)

mtx, samples = getMatrix(mtx_file)

#groupdic = getgroup(group_file, groupname)

celltypes = filter(lambda x: '.txt' not in x, os.listdir(degdir))

for i in xrange(len(celltypes)):
	celltype = celltypes[i]
#	sams = groupdic[celltypes]

	sdegdir = degdir + celltype + '/'
	degresults = filter(lambda x: '.txt' in x, os.listdir(sdegdir))

	alldegs = []; countdegs = dict()
	for result in degresults:
		degfile = sdegdir + result
		degs = getDEG(degfile)
		for deg in degs:
			if not countdegs.has_key(deg):
				countdegs[deg] = 0
			countdegs[deg] += 1
			if deg not in alldegs:
				alldegs.append(deg)
	realdegs = []
	for deg in countdegs.keys():
		count = countdegs[deg]
		if count >= len(celltypes)*7.0/10:
			realdegs.append(deg)
	print celltype, len(realdegs)
	outfile = open(outputdir + celltype + '.degs.tpm.mtx.txt', 'w')
	outfile.write('ID\tsymbol\ttype\tgenetype\t'+'\t'.join(samples))
	lines = ''
	for gene in mtx.keys():
		expdic = mtx[gene]
		geneid = gene[0]; symbol = gene[1]; type = gene[2]; genetype = gene[3]
		if geneid not in realdegs: continue
		line = '\t'.join(gene)
		for sample in samples:
			line += '\t'+str(expdic[sample])
		lines += '\n'+line
	outfile.write(lines)
	outfile.close()
