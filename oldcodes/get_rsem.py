import sys, os
import module as mdl
import numpy as np
#def getIsoform_pct(filename):
#	infile = open(filename, 'r')
#	lines = infile.readlines(); infile.close()
#	tpmdict = dict()
#	for i in xrange(len(lines)):
#		line = lines[i]
#		if not line.startswith('transcript_id'):
#			tmp = line.split('\t')
#			trxid = tmp[0]; geneid = tmp[1]; tpm = float(tmp[-3])
#			if not tpmdict.has_key(geneid):
#				tpmdict[geneid] = {

inputdir = sys.argv[1]
targets = sys.argv[2]
gtf = mdl.getGtf(sys.argv[3])
outfile = open(sys.argv[3],'w')
Tgeneids = targets.split(',')
samples = filter(lambda x: not 'log' in x, os.listdir(inputdir))
samdict = dict()
for i in xrange(len(samples)):
	sample = samples[i]
	sinputdir = inputdir + sample + '/'
	result = open(sinputdir + filter(lambda x: '.isoforms.results' in x, os.listdir(sinputdir))[0], 'r')
	lines = result.readlines(); result.close()
	for i in xrange(len(lines)):
		line = lines[i]
		if not line.startswith('transcript_id'):
			tmp = line.split('\t')
			trxid = tmp[0]; geneid = tmp[1]; tpm = float(tmp[-3])
			if not samdict.has_key(geneid):
				samdict[geneid] = {}
			if not samdict[geneid].has_key(trxid):
				samdict[geneid][trxid] = {}
			samdict[geneid][trxid][sample] = tpm

outfile.write('transcript_id\tgene_id\tTPM_mean\tIso_pct\n')
genedict = dict()
for geneid in samdict.keys():
	trx = samdict[geneid]
	for trxid in trx.keys():
		tpms = trx[trxid].values()
		meanTpm = np.mean(tpms)
		if not genedict.has_key(geneid):
			genedict[geneid] = {}
		genedict[geneid][trxid] = meanTpm
	totlTpm = sum(genedict[geneid].values())
	for trxid in genedict[geneid].keys():
		if meanT == 0:
			isopct = 0
		else:
			isopct = meanT/totalTpm*100
		if not genedict.has_key(geneid):
			genedict[geneid] = {}
		genedict[geneid][trxid] = isopct
		outfile.write(trxid+'\t'+geneid+'\t'+str(meanT)+'\t'+str(isopct)+'\n')
maintrxs = []
for geneid in genedict.keys():
	trxs = genedict[geneid]
	maintrxid = filter(lambda x: trxs[x] == max(trxs.values()), trxs.keys())[0]
	maintrxs.append(maintrxs)

newgtf = dict()
for chr in gtf.keys():
	trxs = gtf[chr]
	mains = filter(lambda x: x.trxid() in maintrxs, trxs)
	for main in mains:
		main.attri('main transcript')
	

#for geneid in Tgeneids:
#	trx = genedict[geneid]
#	totalTpm = sum(trx.values())
#	for trxid in trx.keys():
#		meanT = trx[trxid]
#		if meanT == 0:
#			isopct = 0
#		else:
#			isopct = meanT/totalTpm*100
#		outfile.write(trxid+'\t'+geneid+'\t'+str(meanT)+'\t'+str(isopct)+'\n')
	

outfile.close()
