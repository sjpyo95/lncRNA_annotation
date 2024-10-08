import sys, os
import commands
inputdir = sys.argv[1]
allsamples = os.listdir(inputdir)
sampleName = sys.argv[2]
outputdir = sys.argv[3]
strand = sys.argv[4]	#rpds or star
samples = filter(lambda x: sampleName in x and not '.txt' in x, allsamples)
gtfs = []
for i in xrange(len(samples)):
	sample = samples[i]
	if strand == 'rpds':
		sinputdir = inputdir + sample + '/3.assembly/rpds_stringtie/'
	else:
		sinputdir = inputdir + sample + '/3.assembly/stringtie/'
	if os.path.exists(sinputdir):
		gtfs.append(sinputdir+filter(lambda x: '.gtf' in x, os.listdir(sinputdir))[0])
#if not sampleName +'_taco_gtf_dir.txt' in allsamples:
gtflistfile = open(inputdir + sampleName +'_taco_gtf_dir.txt', 'w')
gtflistfile.write('\n'.join(gtfs))
gtflistfile.close()

if os.path.exists(outputdir):
	print 'output directory already exist!'
	exit()
gtfdir = inputdir + sampleName +'_taco_gtf_dir.txt'
taco = '/export/home/sjpyo/packages/taco/taco-v0.7.3.Linux_x86_64/taco_run -v --gtf-expr-attr FPKM --filter-min-length 100 --filter-min-expr 0.1 --isoform-frac 0.01 --max-isoforms 10 -p 8 -o ' + outputdir + ' ' + gtfdir

print taco+'\n'
commands.getoutput(taco)

