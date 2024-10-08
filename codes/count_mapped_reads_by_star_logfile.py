import sys, os

def parse_Logfile(filename):
	filein = open(filename, 'r')
	lines = filein.readlines(); filein.close()
	mapdic = dict()
	for i in range(len(lines)):
		line = lines[i].strip()
		tmp = line.split('|')
		if len(tmp) < 2: continue
		name = tmp[0].strip(); value = tmp[1].strip()
		mapdic[name] = value

	return mapdic


inputdir = sys.argv[1]
info_wanted = sys.argv[2]
out_file = sys.argv[3]

samples = filter(lambda x: '.txt' not in x, os.listdir(inputdir))

info_wanted = ' '.join(info_wanted.split('_'))
outfile = open(out_file, 'w')
outfile.write('celltype\tvalue\n')

cellinfo = dict()

for i in range(len(samples)):
	sample = samples[i]
	infile = inputdir + sample + '/2.mapping/rpds_star2/' + filter(lambda x: 'Log.final.out' in x, os.listdir(inputdir + sample + '/2.mapping/rpds_star2/'))[0]
	celltype = sample.split('_rep')[0]
	mapinfo = parse_Logfile(infile)
	info_val = mapinfo[info_wanted]
	if not cellinfo.has_key(celltype):
		cellinfo[celltype] = []
	cellinfo[celltype].append(float(info_val))

for celltype in cellinfo.keys():
	vals = cellinfo[celltype]
	val = sum(vals)/len(vals)
	outfile.write(celltype+'\t'+str(val)+'\n')
outfile.close()
