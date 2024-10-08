import sys
import module as mdl

gofiles = sys.argv[1]
gofiles = gofiles.split(',')
def getGOterms(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	go = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip()
		if len(line) == 0 : continue
		col = line.split('\t')
		source = col[0].split(':')[-1]; termName = col[1]; neglog10padj = col[4]
		genes = col[-1]
		if float(neglog10padj) < 2: continue
		if not go.has_key(source):
			go[source] = []
		go[source].append(termName)
#		print termName
	return go
goRdic = {'BP':[], 'MF':[], 'KEGG':[]}
for i in xrange(len(gofiles)):
	gofile = gofiles[i]
	go = getGOterms(gofile)
	for source in go.keys():
		goRdic[source].append(go[source])
for source in goRdic.keys():
	terms = goRdic[source]
	comTerms = list(set.intersection(*map(set, terms)))
	print source
	print comTerms
	print len(comTerms)
	print '\n'


