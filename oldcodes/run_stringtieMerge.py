import sys, os
import config as cf
import commands
stringtiePath = '/share/apps/programs/stringtie/stringtie-1.3.4a.Linux_x86_64/stringtie'
inputdir = sys.argv[1] + '/'
outputdir = sys.argv[2]
server = sys.argv[3]
jobN = int(sys.argv[4])
clu = sys.argv[5]	#all, rep or single
if not os.path.exists(outputdir) : os.makedirs(outputdir)
pre_samples = os.listdir(inputdir)
pre_samples = filter(lambda x: not '.t' in x and not '.out' in x and not '.OU' in x and not '.log' in x, pre_samples)

typedic = dict()
for j in xrange(len(pre_samples)):
	pre_sample = pre_samples[j]
	if clu == 'all':
		celltype = 'all'
	elif clu == 'single':
		celltype = pre_sample
	else:
		celltype = pre_sample.split(clu)[0]
	if not typedic.has_key(celltype):
		typedic[celltype] = []
	typedic[celltype].append(pre_sample)
#print typedic
#exit()
x = 0
for celltype in typedic.keys():
	x = 0
	gtflist = []
	samples = typedic[celltype]
	print samples
	for i in xrange(len(samples)):
		sample = samples[i]
		if os.path.exists(inputdir+sample+'/3.assembly/rpds_stringtie/'):
			sinputdir = inputdir + sample + '/3.assembly/rpds_stringtie/'
		else:
			sinputdir = inputdir + sample + '/3.assembly/stringtie/'
		if os.path.exists(sinputdir):
			gtfFile = os.listdir(sinputdir)
			gtfFile = filter(lambda x: '.gtf' in x and not 'chr' in x, gtfFile)[0]
			gtflist.append(sinputdir + gtfFile)
		else:
			print 'Cannot find GTF file in your directory:', sinputdir
			exit()
		if clu == 'single':
			soutputdir = outputdir + sample + '/stringtieMerge/'
			if not os.path.exists(soutputdir): os.makedirs(soutputdir)
			outfile = soutputdir + sample + '_strMerge.gtf'
			comm = 'cp ' + sinputdir + gtfFile + ' > ' + outfile
			print comm
			print commands.getoutput(comm)
			
	
	if not clu == 'single':
		soutputdir = outputdir + celltype+'/stringtieMerge/'
		if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	
		merge_cmd = stringtiePath + ' --merge ' + ' '.join(gtflist) + ' -o ' + soutputdir +celltype+'_strMerge.gtf -m 50 -c 0 -F 0.1 -T 0 -f 0.1'
		queue = cf.queue(server, x)
		cf.qsub_time('strMerge', jobN, 100)
		cf.qsub_execute('strMerge_'+celltype, queue, merge_cmd, soutputdir+ 'strMerge.log')
		x += 1
