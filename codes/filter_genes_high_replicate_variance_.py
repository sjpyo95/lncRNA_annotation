import sys, os
import module as mdl
import numpy as np
import math

def getColdata(filename, groupname):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	cols = lines[0].strip().split('\t')
	colindex = cols.index(groupname)
	coldic = dict()
	for i in range(1,len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		sam = tmp[1]
		group = tmp[colindex]
		if not coldic.has_key(group):
			coldic[group] = []
		coldic[group].append(sam)
	return coldic

cv = lambda x: np.std(x, ddof=1) / np.mean(x) * 100

mtx_file = sys.argv[1]
col_file = sys.argv[2]
out_file = sys.argv[3]

mtx, samples = mdl.getinfoMatrix(mtx_file)
coldic = getColdata(col_file, 'group1')

newMtx = dict()
highVarN = 0
print 'Total genes:', len(mtx.keys())
for gene in mtx.keys():
	expdic = mtx[gene]
	expvar = []
	for group in coldic.keys():
		gsams = coldic[group]
		exps = [expdic[x] for x in gsams]
		var = cv(exps)
		if math.isnan(var): continue
		expvar.append(var)
	avgV = np.mean(expvar)
	if math.isnan(avgV): continue
	if avgV >= 200: highVarN += 1; continue
	newMtx[gene] = expdic

print 'Removed genes:',highVarN

mdl.writeinfoMatrix(newMtx, samples, out_file)
