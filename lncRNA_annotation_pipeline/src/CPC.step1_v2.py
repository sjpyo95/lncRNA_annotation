import sys, os
import time, commands
import module as mdl
import config as cf
species = sys.argv[1]	#human or mouse
inputgtf = sys.argv[2]
outputdir = sys.argv[3]
txnum = sys.argv[4]
server = sys.argv[5]
jobN = int(sys.argv[6])
if species == 'human': nib_dir = '/home/sjpyo/new/data/fasta/human/hg19/blat/';
if species == 'mouse': nib_dir = '/home/sjpyo/new/data/fasta/mouse/blat/';

gtfdic = mdl.getGtf(inputgtf)
#if species == 'human': chrs = map(lambda x: 'chr'+x, map(str, range(1, 23)) + ['X','Y'])
#elif species == 'mouse': chrs = map(lambda x: 'chr'+x, map(str, range(1, 20)) + ['X','Y'])

fadir = outputdir + '/fasta/'
if not os.path.exists(fadir): os.makedirs(fadir)
logdir = outputdir + '/logs/'
if not os.path.exists(logdir): os.makedirs(logdir)
outname = inputgtf.split('/')[-1][:-4]
i = 1; j = 1; tnum = 0
for chrom in gtfdic.keys():
	trxs = gtfdic[chrom]
	for trx in trxs:
		tnum += 1
		if i == 1:
			fafile = fadir + outname + '.' + str(j) + '.fa'
			outopen = open(fafile, 'w')
		trxid = trx.trxid(); coord = str(trx); sense = trx.sense()
		exons = trx.exons()
		exonseq = mdl.getSeq_exon(exons, sense, nib_dir)
		outopen.write('>' + trxid + '|' + coord + '\n' + exonseq + '\n')
		if i == int(txnum):
			fafile = fadir + outname + '.' + str(j) + '.fa'
			outfile = outputdir + outname + '.result.' + str(j) + '.txt'
			logfile = open(outputdir + outname + '.'+str(j)+ '.run_commandline.log','w')
			workdir = outputdir + 'tmp.' + str(j) + '/'
			if not os.path.exists(workdir): os.makedirs(workdir)
			evdfile = outputdir + outname + '.evd.' + str(j)
			run_predict = '/export/home/sjpyo/new/codes/lncRNA_annotation/cpc-0.9-r2/bin/run_predict.sh ' + fafile + ' ' + outfile + ' ' + workdir + ' ' + evdfile
			logfile.write(run_predict)
			logfile.close()
#			print run_predict + '\n'
			queue = cf.queue(server, j)
			cf.qsub_time('cpcPYO', jobN, 300)
			cf.qsub_execute('cpcPYO_'+str(tnum), queue, run_predict, logdir+str(tnum)+'.log')
#			commands.getoutput(run_predict)
			i = 0; j += 1
		i += 1
		
print tnum
fafile = fadir + outname + '.' + str(j) + '.fa'
outfile = outputdir + outname + '.result.' + str(j) + '.txt'
logfile = open(outputdir + outname + '.'+str(j)+ '.run_commandline.log','w')
workdir = outputdir + 'tmp.' + str(j) + '/'
if not os.path.exists(workdir): os.makedirs(workdir)
evdfile = outputdir + outname + '.evd.' + str(j)
run_predict = '/export/home/sjpyo/new/codes/lncRNA_annotation/cpc-0.9-r2/bin/run_predict.sh ' + fafile + ' ' + outfile + ' ' + workdir + ' ' + evdfile
logfile.write(run_predict)
logfile.close()
queue = cf.queue(server, j)
cf.qsub_time('cpcPYO', jobN, 300)
cf.qsub_execute('cpcPYO_last_'+str(tnum), queue, run_predict, logdir+str(tnum)+'.log')
#print run_predict + '\n'
#commands.getoutput(run_predict)

