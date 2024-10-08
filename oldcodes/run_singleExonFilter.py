import sys, os
import config as cf
exjFilter = '/export/home/sjpyo/new/codes/cafe/backup/single_junctionSup.py'
endFilter = '/export/home/sjpyo/new/codes/cafe/backup/single_endfilter.py'

inputdir = sys.argv[1]
bedfile = sys.argv[2]
type = sys.argv[3]
outputdir = sys.argv[4]
jobN = int(sys.argv[5])
hm_chrs = cf.hm_chrs

for i in xrange(len(hm_chrs)):
	chr = hm_chrs[i]
	gtf = filter(lambda x: chr+'_' in x or chr+'.' in x, os.listdir(inputdir))[0]
	print gtf
	if type == 'exj':
		command = '/usr/bin/python ' + exjFilter + ' ' + inputdir + gtf + ' ' + bedfile + ' ' + outputdir + ' ' + chr
	elif type == 'end':
		command = '/usr/bin/python ' + endFilter + ' ' + inputdir + gtf + ' ' + outputdir + ' ' + chr
	
	queue = cf.queue('biglab', i)
	cf.qsub_time('sTrxFilter', jobN, 100)
	cf.qsub_execute('sTrxFilter', queue, command, outputdir)
