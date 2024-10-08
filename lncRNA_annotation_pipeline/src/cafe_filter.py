import sys, os
import config as cf

assembly = 'hg19'
chrs = cf.hm_chrs
cafe_filter = '/usr/bin/python /home/sjpyo/new/codes/cafe/exj/filter.py'

inputdir = sys.argv[1]
outputdir = sys.argv[2]

samples = filter(lambda x: not '.txt' in x and not 'nohup' in x, os.listdir(inputdir))

for i in range(len(samples)):
	sample = samples[i]
	gtfs= filter(lambda x: '.gtf' in x, os.listdir(inputdir+sample+'/'))
	if len(gtfs) < 24: print 'GTF files chromosome missing'; exit()
	soutputdir = outputdir + '4.5.filter/' + sample + '/'
	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	for chr in chrs:
		if len(filter(lambda x: 'transcripts_'+chr+'.' in x, gtfs)) == 0: continue
		gtffile = inputdir + sample + '/' + filter(lambda x: 'transcripts_'+chr+'.' in x, gtfs)[0]
		options = assembly + ' ' + gtffile + ' ' + chr + ' ' + soutputdir
		command = cafe_filter + ' ' + options
		queue = cf.queue('biglab', i)
		cf.qsub_time('sjpyo', 60, 20)
		cf.qsub_execute('sjpyo', queue, command, soutputdir + sample + '.' + chr + '.')
