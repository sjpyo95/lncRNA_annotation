import sys, os
import re
import config
import commands
regtoolsPath = '/export/home/sjpyo/new/packages/regtools/build/regtools'

#inputdir = sys.argv[1] + '/'
inputdirs = sys.argv[1]
inputdirs = inputdirs.split(',')
outputdir = sys.argv[2] + '/'
strand = sys.argv[3]	#0 = unstranded/XS, 1 = first-strand/RF, 2 = second-strand/FR
if not os.path.exists(outputdir): os.makedirs(outputdir)
#job = sys.argv[4]
server = sys.argv[4]
jobN = int(sys.argv[5])
#samples = filter(lambda x: not '.t' in x and not '.out' in x and not '.log' in x, os.listdir(inputdir))
for i in xrange(len(inputdirs)):
	inputdir = inputdirs[i]
#for i in xrange(len(samples)):
#	sample = samples[i]
#	sinputdir = inputdir + sample + '/2.mapping/'
	sinputdir = inputdir + '/2.mapping/'
	if 'rpds_star' in os.listdir(sinputdir):
		bamfile = sinputdir + '/rpds_star/'+filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(sinputdir+'/rpds_star/'))[0]
		bam = filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(sinputdir+'/rpds_star/'))[0].split('_Align')[0]
	elif 'cutadpt_star' in os.listdir(sinputdir):
		bamfile = sinputdir + '/cutadpt_star/'+filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(sinputdir+'/cutadpt_star/'))[0]
		bam = filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(sinputdir+'/cutadpt_star/'))[0].split('_Align')[0]
	else:
		bamfile = sinputdir + '/star/'+filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(sinputdir+'/star/'))[0]
		bam = filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(sinputdir+'/star/'))[0].split('_Align')[0]
		
#bamfiles = filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(inputdir))
#for i in xrange(len(bamfiles)):
#	bamfile =bamfiles[i]
#	bam = bamfile.split('_')[0]
#pre_samples = filter(lambda x: not '.t' in x and not '.log' in x and not '.OU' in x and not '.out' in x, os.listdir(inputdir))
#typedic = dict()
#for j in xrange(len(pre_samples)):
#	pre_sample = pre_samples[j]
#	celltype = 'all'
##	celltype = pre_sample.split('_')[0]
##	celltype = pre_sample.split('_rep')[0]
#	if not typedic.has_key(celltype):
#		typedic[celltype] = []
#	typedic[celltype].append(pre_sample)
#
#for celltype in typedic.keys():
#	samples = typedic[celltype]
#	print samples
##	soutputdir = outputdir + celltype + '/junction_reads/'
##	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
##samples = os.listdir(inputdir)
##samples = filter(lambda x: not '.t' in x, samples)
##print samples
#
#	for i in xrange(len(samples)):
#		sample = samples[i]
#		print sample
#		sinputdir = inputdir + sample + '/2.mapping/star/'
##		sinputdir = inputdir + sample + '/2.mapping/rpds/'
#		soutputdir = outputdir + '/junction_reads/'
#		if not os.path.exists(soutputdir): os.makedirs(soutputdir)
#
#		if os.path.exists(sinputdir):
##		print i, sample
#			bamfile = filter(lambda x: '.bam' in x and '.bai' not in x, os.listdir(sinputdir))[0]
##			bamfile = filter(lambda x: '.bam' in x and '.bai' not in x and 'rpds' in x, os.listdir(sinputdir))[0]
#			bam  = bamfile.split('.')[0]
	regtools = regtoolsPath + ' junctions extract -m 26 -M 1000000 -s '+ strand + ' -o ' + outputdir + bam + '_exj.bed ' + bamfile
#			print regtools + '\n'
	logfile = outputdir + bam+'.'
			#commands.getoutput(regtools)
	queue = config.queue(server, i)
	config.qsub_time('regtools', jobN, 100)
	config.qsub_execute('regtools', queue, regtools, logfile)
