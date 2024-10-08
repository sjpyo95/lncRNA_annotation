import sys, os
import module as mdl
import math
import numpy as np

def getGroupinfo(filename, groupname):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	col = lines[0].strip().split('\t')
	groupindex = col.index(groupname)
	groupdic = dict()
	groups = []
	for i in range(1, len(lines)):
		line = lines[i]
		tmp = line.strip().split('\t')
		sample = tmp[0]
		app_group = tmp[groupindex]
		if app_group == 'None': continue
		if not groupdic.has_key(app_group):
			groupdic[app_group] = []
		groupdic[app_group].append(sample)
#		if not app_group in groups:
#			groups.append(app_group)
	return groupdic#, groups


def getdegMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[5:]
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; lgfc = float(line[4]);#padj = float(line[5])
		exp = line[5:]; type = line[2]; genetype = line[3]
		if 'PAR_Y' in id: continue
		elif lgfc < 2 : continue
#		elif padj > 0.05 : continue
		if not geneDic.has_key((id,symb,type,genetype,lgfc)):
			geneDic[(id,symb,type,genetype,lgfc)] = dict()
		for x in xrange(len(exp)):
			sample = samples[x]; express = float(exp[x])
			geneDic[(id,symb,type,genetype,lgfc)][sample] = express
	return geneDic, samples

def getTops(deglist, stand, topN):
	
	if stand == 'expression':
		deglist = sorted(sorted(deglist, key = lambda x: float(x[1]), reverse = True), key = lambda x: x[-1], reverse = False)
	elif stand == 'logFC':
		deglist = sorted(sorted(deglist, key = lambda x: float(x[2]), reverse = True), key = lambda x: x[-1], reverse = False)
	#deglist = deglist[:topN*5]
	#deglist = mdl.comp(deglist, lambda x: x[-1], False)
	return deglist[:topN]

def writeTopmtx(toplist, mtx, soutputdir, type, sams):
	outfile = open(soutputdir + type+'_degs.top'+str(top)+'.txt', 'w')
	outfile.write('ID\tsymbol\ttype\tgenetype\tlogFC\tCV\tcelltype\t'+'\t'.join(sams))
	lines = ''
	check = []
	for i in toplist:
		gene = i[0]; meanexp = i[1]; cv = i[-2]; group = i[-1]
		geneid = gene[0]
		mtxgene = filter(lambda x: x[0] == geneid, mtx.keys())[0]
		expdic = mtx[mtxgene]
		if geneid in check: continue
		line = '\t'.join([str(x) for x in gene])+'\t'+str(cv)+'\t'+group
		for sample in sams:
			line += '\t'+str(expdic[sample])
		lines += '\n'+line
		check.append(geneid)
	outfile.write(lines)
	outfile.close()

inputdir = sys.argv[1]
groupname = sys.argv[2]
group_file = sys.argv[3]
outputdir = sys.argv[4]
standard = sys.argv[5]	#expression or logFC
top = int(sys.argv[6])
allmtxfile = sys.argv[7]
allmtx, samples = mdl.getinfoMatrix(allmtxfile)

if not os.path.exists(outputdir): os.makedirs(outputdir)

infiles = filter(lambda x: '.txt' in x, os.listdir(inputdir))
groups = getGroupinfo(group_file, groupname)
allnoveldegs = []; allknowndegs = []
for i in range(len(infiles)):
	infile = inputdir + infiles[i]
	group = infiles[i].split('.')[0]
	print group
	mtx, samples = getdegMatrix(infile)
	groupsams = groups[group]
	allgroupsams = sum(groups.values(), [])
	pcgGenes = []; lncknownGenes = []; lncnovelGenes = []
	for info in mtx.keys():
		type = info[2]; lgfc = info[-1]
		expdic = mtx[info]
		exps = [expdic[x] for x in groupsams]
		meanexp = sum(exps)/len(exps)
		if meanexp <= 10 : continue
		cv_val = np.std(exps, ddof=1)/np.mean(exps)*100

		if type == 'PCG':
			pcgGenes.append([info, meanexp, lgfc, cv_val, group])

		elif type == 'Known_lncRNA':
			lncknownGenes.append([info,meanexp,lgfc, cv_val, group])

		elif type == 'Novel_lncRNA':
			lncnovelGenes.append([info,meanexp, lgfc, cv_val, group])
	
	pcgTop = getTops(pcgGenes, standard, top)
	lncknownTop = getTops(lncknownGenes, standard, top)
	lncnovelTop = getTops(lncnovelGenes, standard, top)
	allknowndegs.extend(lncknownTop)
	allnoveldegs.extend(lncnovelTop)
#	soutputdir = outputdir + group + '/'
#	if not os.path.exists(soutputdir) : os.makedirs(soutputdir)
#	writeTopmtx(pcgTop, mtx, soutputdir, 'PCG', samples)
#	writeTopmtx(lncknownTop, mtx, soutputdir, 'known_lncRNA', samples)
#	writeTopmtx(lncnovelTop, mtx, soutputdir, 'novel_lncRNA', samples)
writeTopmtx(allknowndegs,allmtx, outputdir, 'known_lncRNA', samples)
writeTopmtx(allnoveldegs, allmtx, outputdir, 'novel_lncRNA', samples)
