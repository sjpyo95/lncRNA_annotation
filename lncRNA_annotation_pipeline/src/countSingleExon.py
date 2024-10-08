import sys, os
import module as mdl
inputdir = sys.argv[1]

samples = filter(lambda x: not '.log' in x, os.listdir(inputdir))

for i in xrange(len(samples)):
	sample = samples[i]
	print sample
	sinputdir = inputdir +'/' + sample + '/stringtieMerge/filt/'
	gtffile = sinputdir + filter(lambda x: '.gtf' in x and 'strMerge' in x, os.listdir(sinputdir))[0]
	gtf = mdl.getGtf(gtffile)
	singleN = 0; multiN = 0
	for chr in gtf.keys():
		trxs = gtf[chr]
		for trx in trxs:
			if len(trx.exons()) == 1:
				singleN += 1
			else:
				multiN += 1
	print 'single exon transcripts:',singleN
	print 'multi exon transcripts:', multiN, '\n'
