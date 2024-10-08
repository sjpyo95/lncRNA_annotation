#!/usr/bin/python
import sys, os
from subprocess import Popen, PIPE, STDOUT
#import commands, time
import config as cf	#parameters for RNA-seq analysis
#stringtiePath = '/share/apps/programs/stringtie/stringtie-1.3.4a.Linux_x86_64/stringtie'
stringtie = '/share/apps/programs/python/2.7.14/bin/python /home/sjpyo/new/codes/qsub_stringtie.py'
#refgenPath = '/export/home/sjpyo/data/anno/human/gencode/hg19/gencode.v19.annotation.gtf'

inputdir = sys.argv[1]+'/'
rpds = sys.argv[2]	#rpds or star
strand = sys.argv[3]	#rf: fr-firststrand or fr: fr-secondstrand
refPath = sys.argv[4]
server = sys.argv[5]
jobN = int(sys.argv[6])
outputdir = sys.argv[7]
tissue = 'stringtie'

samples = os.listdir(inputdir)
samples = filter(lambda x: not '.t' in x and not '.out' in x and not '.OU' in x and not '.log' in x, samples)
for i in xrange(len(samples)):
	sample = samples[i]
	if rpds == 'rpds':
		sinputdir = inputdir+sample+'/2.mapping/rpds_star/'
	else:
		if os.path.exists(inputdir+sample+'/2.mapping/cutadpt_star/'): sinputdir = inputdir+sample+'/2.mapping/cutadpt_star/'
		else: sinputdir = inputdir+sample+'/2.mapping/star/'
#	sinputdir = inputdir + sample + '/'
	if os.path.exists(sinputdir):
		if rpds == 'rpds':
			bamfile = filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(sinputdir))[0]
		else:
#		if len(filter(lambda x:'.bam' in x and not '.bai' in x, os.listdir(sinputdir))) == 0: continue
			bamfile = filter(lambda x:'.bam' in x and not '.bai' in x, os.listdir(sinputdir))[0]
		if rpds == 'rpds':
			soutputdir = inputdir + sample + '/3.assembly/rpds_stringtie/'
		else:
			soutputdir = inputdir + sample + '/3.stringtie/'
#		soutputdir = outputdir + sample+ '/3.stringtie/'
		if not os.path.exists(soutputdir): os.makedirs(soutputdir)
		queue = cf.queue(server, i)
		qsub_stringtie = stringtie + ' ' + sinputdir + bamfile + ' ' + strand + ' ' + refPath + ' ' + soutputdir
#		stringtie=stringtiePath+' '+sinputdir+bamfile+' -p '+ config.runThread + ' -v -f '+config.minIsoformFraction+' -m '+config.minTranscriptLength+' -a '+config.minAnchorLength+' -j '+config.minJunctionCoverage+' -c '+config.minReadCoverage+' -M '+config.maxMultiCoverage+' -g '+config.gapBetweenBundle+' -G ' + refgenPath + ' --fr -o '+soutputdir+sample+'_stringtie.gtf'
		cf.qsub_time(tissue, jobN, 300)
		cf.qsub_execute(tissue, queue, qsub_stringtie, soutputdir)
