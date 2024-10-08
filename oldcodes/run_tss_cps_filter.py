import sys, os
import config as cf
import commands, time
tssPath = '/usr/bin/python /export/home/sjpyo/new/codes/cafe/exj/tss.v2.py'
cpsPath = '/usr/bin/python /export/home/sjpyo/new/codes/cafe/exj/cps.py'
filterPath = '/usr/bin/python /export/home/sjpyo/new/codes/cafe/exj/filter.py'

inputdir = sys.argv[1]
if not inputdir[-1] == '/': inputdir += '/'	
outputdir = sys.argv[2]
assembly = sys.argv[3]
if not os.path.exists(outputdir): os.makedirs(outputdir)
job = sys.argv[4]
type = sys.argv[5]
server = sys.argv[6]
jobN = int(sys.argv[7])
#samples = filter(lambda x: not '.txt' in x and not '.log' in x and not '.out' in x, os.listdir(inputdir))

#for i in xrange(len(samples)):
#	sample = samples[i]
#	sinputdir = inputdir + sample + '/'
#	print sinputdir
#	hm_chrs = cf.hm_chrs
#	ms_chrs = cf.ms_chrs
#	soutputdir = outputdir + sample + '/' + job + '_'+ type + '/'
#	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
if assembly == 'hg19':
	chroms = cf.hm_chrs
elif assembly == 'mm10':
	chroms = cf.ms_chrs
if job == 'tss':
	print 'CAFE TSS update start'
	for j in xrange(len(chroms)):
		chr = chroms[j]
		queue = cf.queue(server, j)
		soutputdir = outputdir + '/tss_'+type+'/'
		if not os.path.exists(soutputdir): os.makedirs(soutputdir)
#		print chr
#		print filter(lambda x: chr+'.' in x, os.listdir(sinputdir+'/exj_update/'))
		gtfFile = inputdir+'/exj_update/' + filter(lambda x: chr+'.' in x, os.listdir(inputdir+'/exj_update/'))[0]
		command = tssPath + ' '+ assembly+ ' '+ gtfFile + ' ' + chr + ' ' + soutputdir + ' ' + type
		cf.qsub_time('tss', jobN, 90)
		cf.qsub_execute('tss_' + chr, queue, command, soutputdir + chr + '.log')
	while commands.getoutput('qstat').count('tss') > 0: time.sleep(40)
if job == 'tss' or job == 'cps':
	print 'CAFE CPS update start'
	for j in xrange(len(chroms)):
		chr = chroms[j]
		queue = cf.queue(server, j)
		soutputdir = outputdir + '/cps_'+type+'/'
		if not os.path.exists(soutputdir): os.makedirs(soutputdir)
		gtfFile = inputdir+'/tss_'+type+'/'+filter(lambda x: chr+ '.' in x and 'tss.gtf' in x, os.listdir(inputdir+'/tss_'+type+'/'))[0]
		command = cpsPath + ' '+ assembly + ' ' + gtfFile + ' ' + chr + ' ' + soutputdir + ' ' + type
		cf.qsub_time('cps', jobN, 90)
		cf.qsub_execute('cps_' + chr, queue, command, soutputdir + chr + '.log')
	while commands.getoutput('qstat').count('cps') > 0: time.sleep(40)
if job == 'tss' or job == 'cps' or job == 'filter':
	print 'CAFE Filter start'
	for j in xrange(len(chroms)):
		chr = chroms[j]
		queue = cf.queue(server, j)
		soutputdir = outputdir + '/filter_'+type+'/'
		if not os.path.exists(soutputdir): os.makedirs(soutputdir)
		gtfFile = inputdir+'/cps_'+type+'/'+filter(lambda x: chr+'.tss.cps.gtf' in x, os.listdir(inputdir+'/cps_'+type+'/'))[0]
		command = filterPath + ' '+assembly + ' ' + gtfFile + ' ' + chr + ' ' + soutputdir
		cf.qsub_time('filter', jobN, 90)
		cf.qsub_execute('filter_' + chr, queue, command, soutputdir + chr + '.log')

#	commandline.close()

		
