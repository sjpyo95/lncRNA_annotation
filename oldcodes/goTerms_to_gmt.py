import sys, os

inputdir = sys.argv[1]
outfile = open(sys.argv[2], 'w')

def getGOterms(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
#	pathway = infile.split('/')[-1].split('.txt')[0]
	go = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip()
		if len(line) == 0 : continue
		col = line.split('\t')
		symbol = col[1]; pathway = '_'.join(col[5].split(' '))
		if not go.has_key(pathway):
			go[pathway] = []
		go[pathway].append(symbol)
#		print termName
	return go

terms = filter(lambda x: '.txt' in x, os.listdir(inputdir))

for i in xrange(len(terms)):
	termfile = terms[i]
	term = getGOterms(termfile)
	for path in term.keys():
		genes = term[path]
		outfile.write(path+'\t.\t'+'\t'.join(genes)+'\n')
outfile.close()

