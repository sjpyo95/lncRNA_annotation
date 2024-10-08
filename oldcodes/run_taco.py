import sys, re, os
import config as cf

inputdir = sys.argv[1] + '/'
outputdir = sys.argv[2]
server = sys.argv[3]
jobN = int(sys.argv[4])
#if not os.path.exists(outputdir): os.makedirs(outputdir)

tacoPath = '/export/home/sjpyo/new/packages/taco-v0.7.3.Linux_x86_64/taco_run'

samples = os.listdir(inputdir)
samples = filter(lambda x: not '.t' in x and not '.out' in x and not '.OU' in x and not 'taco' in x, samples)

#typedic = dict()
#for j in xrange(len(pre_samples)):
#	pre_sample = pre_samples[j]
#	celltype = pre_sample.split('_')[0]
#	celltype = 'all'
#	if not typedic.has_key(celltype):
#		typedic[celltype] = []
#	typedic[celltype].append(pre_sample)
#print typedic
#exit()
#x = 0
#for celltype in typedic.keys():
gtflist = []
#	samples = typedic[celltype]
#	print samples
for i in xrange(len(samples)):
	sample = samples[i]
#	print sample
	sinputdir = inputdir + sample + '/3.stringtie/'
	if os.path.exists(sinputdir):
		gtfFile = os.listdir(sinputdir)
		if len(filter(lambda x: 'stringtie.gtf' in x and not 'chr' in x, gtfFile)) < 1:
			print sample + 'no GTF file'
			exit()
	
		else:
			gtfFile = filter(lambda x: 'stringtie.gtf' in x and not 'chr' in x, gtfFile)[0]
#			print sample
#		print gtfFile
		gtflist.append(sinputdir + gtfFile)
exit()
#	print gtflist
gtflistfile = open(inputdir + 'taco_gtf_dir.txt', 'w')
gtflistfile.write('\n'.join(gtflist))
gtflistfile.close()
#soutputdir = outputdir + celltype + '/'
#soutputdir = outputdir + '/'
gtfDir = inputdir+'taco_gtf_dir.txt'
taco = tacoPath + ' --gtf-expr-attr FPKM --filter-min-length 50 --filter-min-expr 0.1 --isoform-frac 0.1 -p 5 -o ' + outputdir + ' ' + gtfDir
queue = cf.queue(server, 0)
cf.qsub_time('taco', jobN, 200)
cf.qsub_execute('taco', queue,  taco, outputdir+'/taco_')
#x += 1
