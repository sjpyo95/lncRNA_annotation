import sys, os
import module as mdl

def getdegMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[6:]
	genelist = []
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; lgfc = float(line[4]); cv = float(line[5])
		exp = line[6:]; type = line[2]; genetype = line[3]
		if 'PAR_Y' in id: continue
		genelist.append(id)
#		elif lgfc < 2 : continue
#		elif padj > 0.05 : continue
	return genelist

def getFeatures(gtf, genel, sample):
	gtfgenes = sum(gtf.values(),[])
	featuredic = dict()
	for gene in gtfgenes:
		geneid = gene.geneid()
		if geneid in genel:
			chr = gene.chr(); symbol = gene.symb()
			start = gene.start(); end = gene.end(); sense = gene.sense()
			genetype = gene.genetype()
			lenth = gene.len()

			trxs = gene.transcripts()
			trxN = len(trxs)
			exonNums = [x.exonNum() for x in trxs]
			exonMean = float(sum(exonNums))/len(exonNums)
			endsups = '|'.join(list(set([x.endsup() for x in trxs])))
			featuredic[geneid] = [symbol, genetype, chr, start, end, sense, lenth, trxN, exonMean, endsups, sample]
	return featuredic

def writeout(features, outfile):
	for id in features.keys():
		f = features[id]
		outfile.write(f[10] + '\t' + id+'\t'+str(f[0])+'\t'+ str(f[1])+ '\t' + str(f[2]) + ':'+str(f[3])+'-'+str(f[4])+ '('+f[5]+')'+'\t'+str(f[6])+'\t'+str(f[7])+'\t'+str(f[8])+'\t'+str(f[9]) + '\n')

inputdir = sys.argv[1]
gtf_file = sys.argv[2]
outputdir = sys.argv[3]
if not os.path.exists(outputdir): os.makedirs(outputdir)
novelout = open(outputdir + 'DEG.novel_lncRNAs.features.txt','w')
knownout = open(outputdir + 'DEG.known_lncRNAs.features.txt', 'w')
novelout.write('celltype\tID\tsymbol\tgenetype\tposition\tlength\ttrxN\texonN mean\tendsup\n')
knownout.write('celltype\tID\tsymbol\tgenetype\tposition\tlength\ttrxN\texonN mean\tendsup\n')


gtf = mdl.getGtf2(gtf_file)

samples = filter(lambda x: '.txt' not in x, os.listdir(inputdir))

for i in range(len(samples)):
	sample = samples[i]
	sinputdir = inputdir + sample + '/'
	novel_file = sinputdir + filter(lambda x: 'novel_lncRNA' in x and '.txt' in x, os.listdir(sinputdir))[0]
	known_file = sinputdir + filter(lambda x: 'known_lncRNA' in x and '.txt' in x, os.listdir(sinputdir))[0]
	novelgenes = getdegMatrix(novel_file)
	knowngenes = getdegMatrix(known_file)
	novelfeatures = getFeatures(gtf, novelgenes,sample)
	knownfeatures = getFeatures(gtf, knowngenes,sample)
#	novelout = outputdir + sample + '.DEG.novel_lncRNAs.features.txt'
#	knownout = outputdir + sample + '.DEG.known_lncRNAs.features.txt'
	writeout(novelfeatures, novelout)
	writeout(knownfeatures, knownout)
novelout.close()
knownout.close()
