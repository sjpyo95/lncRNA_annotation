import sys, os
import module as mdl

#def calStat(res):
#	typenum = {'PCG': 0, 'Known_lncRNA':0,'Novel_lncRNA':0,'others':0}
#	for type in res.keys():
#		num = len(res[type])
#		if type == 'PCG': typenum['PCG'] +=1
#		elif type == 'Known_lncRNA': typenum['Known_lncRNA'] += 1
#		elif type == 'Novel_lncRNA': typenum['Novel_lncRNA'] += 1
#		else:
#			typenum['others'] += 1
#	return typenum
		
def parseRes(filename, cutfc, cutp):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	degdic = dict()
	for i in range(1,len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		id = tmp[0]; symbol = tmp[1]; type = tmp[2]; baseMean = tmp[3]
		fc = float(tmp[4]); padj = float(tmp[5])
		if fc < cutfc : continue
		elif padj > cutp : continue
		if not degdic.has_key(type):
			degdic[type] = []
		degdic[type].append(symbol)

	return degdic
#def getinfoDegMatrix(filename):
#	filein = open(filename)
#	lines = filein.readlines(); filein.close()
#	header = lines[0].strip()
#	samples = header.split('\t')[6:]
#	geneDic = dict()
#	for i in xrange(1,len(lines)):
#		line = lines[i].strip().split('\t')
#		id = line[0]; symb = line[1]; exp = line[6:]; type = line[2]; genetype = line[3]; logFC = line[4]
#		if 'PAR_Y' in id: continue
#		if not geneDic.has_key((id,symb,type,genetype)):
#			geneDic[(id,symb,type,genetype)] = dict()
#		for x in xrange(len(exp)):
#			sample = samples[x]; express = float(exp[x])
#			geneDic[(id,symb,type,genetype)][sample] = express
#	return geneDic, samples


inputdir = sys.argv[1]
outfilename = sys.argv[2]
cutfc = float(sys.argv[3])
cutp = float(sys.argv[4])
filenames = filter(lambda x: '.result.txt' in x, os.listdir(inputdir))
statdic = dict()
for i in range(len(filenames)):
	filename = filenames[i]
	infile = inputdir + filename
	degdic = parseRes(infile, cutfc, cutp)
#	typenum = calStat(degdic)
	sampleName = filename.split('.result.')[0]
	statdic[sampleName] = degdic
#	statdic[sampleName] = typenum

outfile = open(outfilename, 'w')
outfile.write('Sample\tPCG\tKnown lncRNA\tNovel lncRNA\t')
#print statdic
for sample in statdic.keys():
	if sample not in statdic.keys(): continue
	typedic = statdic[sample]
	types = typedic.keys()
	pcg = 0; klnc = 0; nlnc = 0
	if 'PCG' in types:
		pcg = len(typedic['PCG'])
	if 'Known_lncRNA' in types:
		klnc = len(typedic['Known_lncRNA'])
	if 'Novel_lncRNA' in types:
		nlnc = len(typedic['Novel_lncRNA'])

	outfile.write('\n'+sample+'\t'+str(pcg)+'\t'+str(klnc)+'\t'+str(nlnc))
outfile.close()

