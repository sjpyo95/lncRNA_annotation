import sys, os
import config as cf
cafe_cps = '/usr/bin/python /home/sjpyo/new/codes/cafe/exj/cps.py'
chrs = cf.hm_chrs
assembly = 'hg19'
inputdir = sys.argv[1]
cps_type = sys.argv[2]
outputdir = sys.argv[3]

samples = filter(lambda x: not '.txt' in x and not 'nohup' in x, os.listdir(inputdir))

for i in range(len(samples)):
	sample = samples[i]
#	print sample
	gtfs = filter(lambda x: '.gtf' in x, os.listdir(inputdir+sample+'/'))
	if len(gtfs) < 24: print 'GTF files chromosome missing'; exit()
	soutputdir = outputdir + '4.4.cps_'+cps_type+'/'+sample+'/'
	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	for chr in chrs:
		if len(filter(lambda x: 'transcripts_'+chr+'.' in x, gtfs)) == 0: continue
		gtffile = inputdir + sample + '/' + filter(lambda x: 'transcripts_'+chr+'.' in x, gtfs)[0]
		options = assembly + ' ' + gtffile + ' ' + chr + ' ' + soutputdir + ' ' + cps_type
		command = cafe_cps + ' ' + options
		queue = cf.queue('biglab', i)
		cf.qsub_time('sjpyo', 60, 20)
		cf.qsub_execute('sjpyo', queue, command, soutputdir + sample + '.' + chr + '.')
