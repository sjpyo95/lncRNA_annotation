import sys, os
import config
gffcmpPath = '/home/sjpyo/new/packages/gffcompare-0.11.6.Linux_x86_64/gffcompare'

infile = sys.argv[1]
ref = sys.argv[2]
outputPrefix = sys.argv[3]#+'/'+assembly+'/'
if not os.path.exists('/'.join(outputPrefix.split('/')[:-1])): os.makedirs('/'.join(outputPrefix.split('/')[:-1]))
server = sys.argv[4]
jobN = int(sys.argv[5])
#samples = filter(lambda x: not '.log' in x and not '.t' in x and not '.OU' in x and not '.out' in x, os.listdir(inputdir))
#listfile = outputdir+assembly+'_gffcmp.txt'
#gtflist = open(listfile, 'w')

#for i in xrange(len(samples)):
#	sample = samples[i]
#	if assembly == 'stringtie':
#		sinputdir = inputdir + sample + '/3.assembly/stringtie/'
#		gtf = sinputdir + filter(lambda x: '.gtf' in x, os.listdir(sinputdir))[0]
#	elif assembly == 'rpds':
#		sinputdir = inputdir + sample + '/3.assembly/rpds_stringtie/'
#	else:
#		sinputdir = inputdir + sample + '/' + assembly + '/'
#	if os.path.exists(sinputdir):
#		print sinputdir
#		if assembly == 'stringtie' and assembly == 'rpds':
#			gtf = sinputdir + filter(lambda x: sample in x and '.gtf' in x  and not 'chr' in x, os.listdir(sinputdir))[0]
#		else:
#			gtf = sinputdir +  filter(lambda x: not 'chr' in x and '.gtf' in x, os.listdir(sinputdir))[0]

#	gtflist.write(gtf+'\n')
#gtflist.close()
#		soutputdir = outputdir + sample + '/'
#		if not os.path.exists(soutputdir): os.makedirs(soutputdir)
gffcmp = gffcmpPath + ' -R -M -r ' + ref + ' -o ' + outputPrefix +' '+ infile
queue = config.queue(server, 2)
config.qsub_time('gffcmp', jobN, 20)
config.qsub_execute('gffcmp', queue, gffcmp, '/'.join(outputPrefix.split('/')[:-1]))

