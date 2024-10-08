import config as cf
import sys, os

rsemPath = '/home/sjpyo/new/packages/RSEM-master/rsem-calculate-expression'
inputdir = sys.argv[1] + '/'
ref = sys.argv[2]
strand = sys.argv[3]
outputdir = sys.argv[4]
rpds = sys.argv[5]
server = sys.argv[6]
jobN = int(sys.argv[7])

#samples = filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(inputdir))
samples = filter(lambda x: not '.bam' in x and not '.t' in x and not '.out' in x and not 'total' in x and not '.log' in x, os.listdir(inputdir))
#bamfiles = []
for i in xrange(len(samples)):
	sample = samples[i]
	if rpds == 'rpds':
		sinputdir = inputdir +sample + '/1.preprocess/rpds_fastqs/'
		fastq1 = sinputdir + filter(lambda x: 'rpds_1.fastq' in x, os.listdir(sinputdir))[0]
		fastq2 = sinputdir + filter(lambda x: 'rpds_2.fastq' in x, os.listdir(sinputdir))[0]
	else:
		sinputdir = inputdir + sample + '/1.preprocess/sickle/'
		fastq1 = sinputdir + filter(lambda x:'_1_trimmed.fastq' in x and not 'single' in x, os.listdir(sinputdir))[0]
		fastq2 = sinputdir + filter(lambda x:'_2_trimmed.fastq' in x and not 'single' in x, os.listdir(sinputdir))[0]
	soutputdir = outputdir + sample + '/'
	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	
	rsem = rsemPath + ' -p 8 --star --no-bam-output --paired-end --forward-prob ' + strand + ' ' + fastq1 + ' ' + fastq2 + ' ' + ref + ' ' + soutputdir + sample + '.rsem'
	queue = cf.queue(server, i)
	cf.qsub_time('rsem', jobN, 60)
	cf.qsub_execute('rsem', queue, rsem, soutputdir)
