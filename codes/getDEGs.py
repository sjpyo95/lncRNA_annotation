import sys, os
import module as mdl

def getDeseq2(filename,cutfc, cutpadj):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	col = lines[0]

	degdic = dict()
	for i in range(1, len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		id = tmp[0]; lgfc = float(tmp[2]); padj = float(tmp[-1])
		if abs(lgfc) < cutfc: continue
		elif padj > cutpadj: continue
		degdic[id] = [lgfc, padj]
	return degdic

def parseGroup(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	col = lines[0].strip()
	groupNo = 3
	groupdic = {}
	for i in range(1,len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		sample = tmp[1]; group = tmp[groupNo]
		if not groupdic.has_key(group):
			groupdic[group] = []
		groupdic[group].append(sample)
	return groupdic
		

degdir = sys.argv[1]

groupfile = sys.argv[2]
groups = parseGroup(groupfile)

cutfc = float(sys.argv[3])
cutpadj = float(sys.argv[4])

info_mtxfile = sys.argv[5]
mtx, samples = mdl.getinfoMatrix(info_mtxfile)

outputdir = sys.argv[6]
if not os.path.exists(outputdir): os.makedirs(outputdir)
files = filter(lambda x: 'deseq2.txt' in x, os.listdir(degdir))

for i in range(len(files)):
	file = degdir + files[i]
	groupname = files[i].split('.deseq2.')[0]
	degs = getDeseq2(file, cutfc, cutpadj)
	print groupname, len(degs)
	outfile = open(outputdir + groupname + '.degs.mtx.txt', 'w')
	outfile.write('ID\tsymbol\ttype\tgenetype\tlog2FC\tpadj\t'+'\t'.join(samples))
	lines = ''
	for info in mtx.keys():
		id = info[0]; symb = info[1]; type = info[2]; genetype = info[3]
		if not id in degs.keys(): continue
		lgfc = degs[id][0]; padj = degs[id][1]
		expdic = mtx[info]
		line = '\t'.join(info)+'\t'+str(lgfc)+'\t'+str(padj)
		for sample in samples:
			line += '\t'+str(expdic[sample])
		lines += '\n'+line
	outfile.write(lines)
	outfile.close()

		
