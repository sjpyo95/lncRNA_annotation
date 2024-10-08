import sys
import module as mdl
import random
matrixfile1 = sys.argv[1]
matrixfile2 = sys.argv[2]
typefile = sys.argv[3]
bmifile = sys.argv[4]
type = sys.argv[5]
outputdir = sys.argv[6]

def getType(filename, type):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	typedic = dict()
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if line.startswith('GTEX'):
			col = line.split('\t')
			samid = col[0]; bodysite = col[-2]; partid = col[-3]
			if bodysite == type:
				if not typedic.has_key(partid):
					typedic[partid] = []
				typedic[partid].append(samid)
#			if not typedic.has_key(bodysite):
#				typedic[bodysite] = dict()
#			if not typedic[bodysite].has_key(partid):
#				typedic[bodysite][partid] = []
#			typedic[bodysite][partid].append(samid)
	return typedic

def getbmi(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	bmidic = dict()
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if line.startswith('GTEX'):
			col = line.split('\t')
			partid = col[0]; bmi = col[1]
			bmidic[partid] = bmi
	return bmidic

def getobesity(typed, bmidic):
	obesitydic = {'normal':[], 'obesity':[]}
	for partid in bmidic.keys():
		bmi = float(bmidic[partid])
		if partid in typed.keys():
			samids = typed[partid]
#			if bmi < 18.5:
#				obesitydic['under']+= samids
	#		elif bmi >= 18.5 and bmi < 25:
			if bmi < 22.5:
				obesitydic['normal'].append((samids, bmi))
#			elif bmi >= 25 and bmi < 30:
#				obesitydic['over']+=samids
#			else:
			elif bmi > 32.5:
				obesitydic['obesity'].append((samids, bmi))
	return obesitydic
def getspecificMatrix(filename, samples):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	cols = lines[0].split('\t')
	samIndexs = map(lambda x: cols.index(x), samples)
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		id = tmp[0]; symb = tmp[1]
		if not geneDic.has_key((id,symb)):
			geneDic[(id,symb)] = dict()
		for x in samIndexs:
			sample = cols[x]; express = tmp[x]
			geneDic[(id,symb)][sample] = express
	return geneDic
		
typedic = getType(typefile, type)
bmidic = getbmi(bmifile)
obesity = getobesity(typedic, bmidic)

outfile1 = outputdir +  type + '.reads.txt'
outfile2 = outputdir + type + '.tpm.txt'
metafile = open(outputdir + 'meta_data.'+type+'.txt', 'w')
metafile.write('sample\tgroup\n')

allsams = []
n = 100000000
for fat in obesity.keys():
	samids_bmi = obesity[fat]
	if n > len(samids_bmi):
		n = len(samids_bmi)
	else:
		n = n
	print fat
	if fat == 'normal':
		or_samids = mdl.comp(samids_bmi, lambda x: (x[1]), False)
	else:
		or_samids = mdl.comp(samids_bmi, lambda x: (x[1]), True)
	sampling = map(lambda x: x[0][0], or_samids[:n])
	bmis = map(lambda x: x[1], or_samids)
	print fat+'\t'+str(len(or_samids))

#	sampling = random.sample(samids, k = 100)
	allsams += sampling
	for samid in sampling:
		metafile.write(samid+'\t'+fat+'\n')
metafile.close()
print 'meta data created'
newmatrix1 = getspecificMatrix(matrixfile1, allsams)
newmatrix2 = getspecificMatrix(matrixfile2, allsams)

mdl.writeMatrix(newmatrix1,outfile1)
mdl.writeMatrix(newmatrix2, outfile2)
