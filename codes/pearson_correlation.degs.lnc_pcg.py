import sys, os
import module as mdl
from scipy import stats

def getMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[4:]
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; lgfc = float(line[4]);#padj = float(line[5])
		exp = [float(x) for x in line[4:]]; type = line[2]; genetype = line[3]
		geneDic[(id,symb,type,genetype)] = exp
	return geneDic


def getdegs(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[7:]
	genes = []
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; lgfc = float(line[4]);#padj = float(line[5])
		exp = [float(x) for x in line[7:]]; type = line[2]; genetype = line[3]
		genes.append((id,symb,type,genetype))
	return genes

def pearson_corrcoef(exps1, exps2):
	return stats.pearsonr(exps1,exps2)

inputdir = sys.argv[1]
allmtxfile = sys.argv[2]
outputdir = sys.argv[3]
if not os.path.exists(outputdir): os.makedirs(outputdir)
infomtx,samples = mdl.getinfoMatrix(allmtxfile)
allmtx = getMatrix(allmtxfile)
pcgmtx = {key: allmtx[key] for key in allmtx.keys() if key[2] == 'PCG'}

mtxfiles = filter(lambda x: '.txt' in x and not 'PCG' in x, os.listdir(inputdir))

outfile = outputdir + 'lncRNA_PCG.high_corr.mtx.txt'

corGenes = []
for i in range(len(mtxfiles)):
	mtxfile = inputdir + mtxfiles[i]
	lncs = getdegs(mtxfile)
#	lnc_pcg = {}
	corGenes.extend(lncs)
	for deg in lncs:
		symb = deg[1]
		degexps = allmtx[deg]
		for pcg in pcgmtx.keys():
			pcgexps = pcgmtx[pcg]
			pcc = pearson_corrcoef(degexps, pcgexps)
			if abs(pcc[0]) > 0.90 and pcc[1] < 0.05:
				if pcg not in corGenes:
					corGenes.append(pcg)

corGenes_mtx = {key: infomtx[key] for key in corGenes}
mdl.writeinfoMatrix(corGenes_mtx, samples, outfile)

#				if not lnc_pcg.has_key('neg'):
#					lnc_pcg['neg'] = dict()
#				if not lnc_pcg.has_key('pos'):
#					lnc_pcg['pos'] = dict()
#
#				if pcc[0] < 0:
#					print symb, pcc
#					if not lnc_pcg['neg'].has_key(symb):
#						lnc_pcg['neg'][symb] = []
#					lnc_pcg['neg'][symb].append([pcg[1],pcc[0]])
#
#				elif pcc[0] > 0:
#					if not lnc_pcg['pos'].has_key(symb):
#						lnc_pcg['pos'][symb] = []
#					lnc_pcg['pos'][symb].append([pcg[1],pcc[0]])
#		if symb in lnc_pcg.keys():
#			corpcgs = lnc_pcg[symb]
#			lnc_pcg[symb] = sorted(corpcgs, key = lambda x: abs(x[1]), reverse = True)
#	for lnc in lnc_pcg['pos'].keys():
#		pcgs = lnc_pcg['pos'][lnc]
#		for p in pcgs:
#			symb = p[0]; score = p[1]
#			wr.writerow([lnc,symb,'p',str(score)])
#	
#	for lnc in lnc_pcg['neg'].keys():
#		pcgs = lnc_pcg['neg'][lnc]
#		for p in pcgs:
#			symb = p[0]; score = p[1]
#			wr.writerow([lnc,symb,'n',str(score)])
