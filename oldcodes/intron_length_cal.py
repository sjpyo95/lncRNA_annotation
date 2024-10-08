import sys
import module as mdl
import numpy as np

infile = sys.argv[1]
cutoff = float(sys.argv[2])
gtf = mdl.getGtf(infile)
intronLen = []
for chr in gtf.keys():
	trxs = gtf[chr]
	for i in xrange(len(trxs)):
		trx = trxs[i]
		introns = trx.introns()	
		intronLen += map(lambda x: abs(float(x.start())-float(x.end())), introns)
maxlen = np.percentile(intronLen, 100-cutoff)
minlen = np.percentile(intronLen, cutoff)

print 'max ' + str(cutoff) + ' : ' + str(maxlen)
print 'min ' + str(cutoff) + ' : ' + str(minlen)
