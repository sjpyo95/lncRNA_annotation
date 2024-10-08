import sys, os
import config as cf
import commands
exjPath = '/home/sjpyo/new/codes/cafe/exj/exj.revised2.py'
inputdir = sys.argv[1]+ '/'
beddir = sys.argv[2]
outputdir = sys.argv[3]
species = sys.argv[4]
exj_type = sys.argv[5]
exj_level = sys.argv[6]

server = sys.argv[7]
jobN = int(sys.argv[8])
logdir = outputdir + 'logs/'
#if not os.path.exists(outputdir): os.makedirs(outputdir)
if not os.path.exists(logdir): os.makedirs(logdir)
#samples = filter(lambda x: not '.gtf' in x and not 'cafe' in x and not '.OU' in x, os.listdir(inputdir))

#for i in xrange(len(samples)):
#	sample = samples[i]
sinputdir = inputdir + '/assembly_chrN/'
#	sinputdir = inputdir + sample + '/'
inputGtfs = [sinputdir + x for x in filter(lambda x: '.gtf' in x, os.listdir(sinputdir))]
#	print inputGtfs
#	beddir = sinputdir + 'bedfiles/'
#	bedFile = [beddir + x for x in filter(lambda x: '.bed' in x, os.listdir(beddir))]
#	bedFile = ','.join(bedFile)
#	beddir = inputdir + sample + '/junction_reads/'
#	bedFile = beddir + filter(lambda x: '.bed' in x and 'merged' in x, os.listdir(beddir))[0]
#	bedFile = '/export/home/sjpyo/sophie/YSU/taco/N/allN_EJ_merged.bed'
if species == 'human' or species == 'hg19':
	chroms = cf.hm_chrs
	assembly = 'hg19'
elif species == 'mouse' or species == 'mm10':
	chroms = cf.ms_chrs
	assembly = 'mm10'
x = 0
soutputdir = outputdir + '/exj_update/'
if not os.path.exists(soutputdir): os.makedirs(soutputdir)
#if not os.path.exists(soutputdir+'/logs/'): os.makedirs(soutputdir+'/logs/')
noneed = []
for chrom in chroms:
	if chrom in noneed: continue
	bedFile = beddir + filter(lambda x: chrom + '.' in x, os.listdir(beddir))[0]
	inputGtf = filter(lambda x: chrom+'.' in x, inputGtfs)[0]
	command = '/usr/bin/python '+ exjPath+ ' '+ assembly+' '+inputGtf+' '+bedFile+' '+chrom+' '+soutputdir + ' ' + exj_type + ' ' + exj_level
#		commands.getoutput(command)
#		print command+'\n'
	queue = cf.queue(server, x)
	cf.qsub_time('exj', jobN, 300)
	logfile = logdir+chrom+'.exj.'
	cf.qsub_execute('exj'+chrom, queue, command, logfile)
	x += 1
