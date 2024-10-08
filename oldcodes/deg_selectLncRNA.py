import sys, os

def comp(samplelist, stand, re):
	return sorted(samplelist, key=stand, reverse = re)

def getGtf(infile):
	inf = open(infile, 'r')
	lines = inf.readlines(); inf.close()
	genes = dict()
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('#'):
			tline = line.split('\t')
			feature = tline[2]
			if feature != 'gene': continue
			tmp = tline[-1]
			geneid = tmp.split('gene_id')[1].split('"')[1].split('.')[0]
			geneN = tmp.split('gene_name')[1].split('"')[1]
			geneN = ''
			genes[geneid] = geneN
#			print geneid
#			print geneN
	return genes

def selectDEG(infile, fcut, pcut, ref, updown):
	rfile = open(infile, 'r')
	lines = rfile.readlines(); rfile.close()
	deglist = []
	for i in xrange(1,len(lines)):
		line = lines[i].strip()
		if i != 0:
			#fc = float(line.split('\t')[1]); 
#			logfc = float(line.split('\t')[2])
			logfc = float(line.split('\t')[5])
#			lfc = float(line.split('\t')[3])
			pval = line.split('\t')[-2]#; padj = line.split('\t')[6]
			if pval == 'NA' : continue
			geneid = line.split('\t')[0]
#			if geneid in pcgid: genetype = 'PCG'
#			elif geneid in lncid: genetype = 'lncRNA'
#		for gene in genes:
#			print gene.genetype()
#			if geneid == gene.geneid().split('.')[0]:
#				genetype = gene.genetype()
#				break
			if not geneid in ref.keys(): continue
			geneN = ref[geneid]
#			print geneN
			if abs(logfc) < fcut: continue
			elif float(pval) > pcut: continue
			if updown == 'up':
				if logfc >= 0 : continue
			elif updown == 'down':
				if logfc < 0: continue
			print geneid+'\t'+geneN
			deglist.append([infile.split('/')[-1].split('.')[0], geneid, geneN, logfc, pval])
	deglist = comp(deglist, lambda x: (x[3], x[4]), True)
	return deglist


inputdir = sys.argv[1]
gtf = sys.argv[2]
fcut = float(sys.argv[3])
pcut = float(sys.argv[4])
top = sys.argv[5]	#num or all
if top != 'all': ttop = top
updown = sys.argv[6] #up or down (regulation)
ref = getGtf(gtf)
samples = filter(lambda x: '.txt' in x and 'deseq' in x, os.listdir(inputdir))
alldeg = []
for i in xrange(len(samples)):
	sample = samples[i]
	print sample
	infile = inputdir + sample
	deglist = selectDEG(infile, fcut, pcut, ref, updown)
	if top == 'all':
		ttop = len(deglist)
		print top
	if updown == 'up':
		degtop = deglist[-int(ttop):]
	elif updown == 'down':
		degtop = deglist[:int(ttop)]
	alldeg.extend(degtop)

outfile = open(inputdir+'allDEGs.lnc.top'+updown + str(top)+'.txt', 'w')
outfile.write('sample\tgeneid\tsymbol\tlogFC\tpval\n')
for deg in alldeg:
	deg = map(str, deg)
	outfile.write('\t'.join(deg) + '\n')
outfile.close()



