import sys, os
import module as mdl
import math

mtxinfofile = sys.argv[1]
groupfile = sys.argv[2]
group = sys.argv[3] #group1, group2 or group3
tpmcutoff = float(sys.argv[4])
logfc_cutoff = float(sys.argv[5])
outputdir = sys.argv[6]
if not os.path.exists(outputdir): os.makedirs(outputdir)
	
mtx, samples = mdl.getinfoMatrix(mtxinfofile)

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

def getLogfc(groupExpdic, maingroup):
	mainExp = groupExpdic[maingroup]
	groupNames = groupExpdic.keys()
	groupNames.remove(maingroup)
	otherExps = [groupExpdic[x] for x in groupNames]
	loglist = []
	for o in otherExps:
		logval = math.log(mainExp/o, 2)
		loglist.append(logval)
	logmean = sum(loglist)/len(loglist)
	return logmean

groupinfo = getGroupinfo(groupfile, group)

grouplogdic = dict()
for info in mtx.keys():
	expsdic = mtx[info]
	groupexpdic = dict()
	for g in groupinfo.keys():
		sams = groupinfo[g]
		groupExps = [expsdic[x] for x in  sams]
		meanExp = sum(groupExps)/len(groupExps)+0.5
		groupexpdic[g] = meanExp

	for mg in groupexpdic.keys():
		logmean = getLogfc(groupexpdic, mg)
		mainExpmean = groupexpdic[mg]
		if mainExpmean < tpmcutoff: continue
		if logmean < logfc_cutoff: continue
		if not grouplogdic.has_key(mg):
			grouplogdic[mg] = dict()
		l = list(info)
		l.append(str(logmean))

		grouplogdic[mg][tuple(l)] = expsdic
	
for mg in grouplogdic.keys():
	degmtx = grouplogdic[mg]
	outfilename = outputdir + mg + '.degs.txt'
	outfile = open(outfilename, 'w')
	outfile.write('ID\tsymbol\ttype\tgenetype\tlogFC\t'+'\t'.join(samples))
	lines = ''
	degN = 0
	for gene in degmtx.keys():
		line = '\t'.join(gene)
		degexpdic = degmtx[gene]
		for sample in samples:
			line += '\t' + str(degexpdic[sample])
		lines += '\n'+line
		degN += 1
	print mg + '\t' + str(degN)
	outfile.write(lines)
	outfile.close()

