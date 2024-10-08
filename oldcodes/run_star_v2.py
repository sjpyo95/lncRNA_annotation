#!/usr/bin/python
import sys, os
from subprocess import Popen, PIPE, STDOUT
import commands, time

import config  #mapping parameters

starPath = '/share/apps/programs/star/2.5.4/bin/Linux_x86_64/STAR'
#starPath = '/home/sjpyo/sjpyo/packages/STAR-2.7.3a/STAR/bin/Linux_x86_64/STAR'
bamIndex = '/usr/bin/python /home/sjpyo/new/codes/index_multiBAMs.py'
rseqcPath = '/home/sjpyo/new/packages/RSeQC-2.6.4/scripts/infer_experiment.py'

jobname = sys.argv[1]
ver = sys.argv[2]	#hg19 or mm10
#queue = sys.argv[3]
rpds = sys.argv[3]
server = sys.argv[4]
jobN = int(sys.argv[5])
#genome = 'hg19'

def intronLen (gtf, cutoff):
	import module as mdl
	import numpy as np

	gtf = mdl.getGtf(gtf)
	intronLen = []
	for chr in gtf.keys():
		trxs = gtf[chr]
		for i in xrange(len(trxs)):
			trx = trxs[i]
			introns = trx.introns()	
			intronLen += map(lambda x: abs(float(x.start())-float(x.end())), introns)
	maxlen = np.percentile(intronLen, 100-cutoff)
	minlen = np.percentile(intronLen, cutoff)
	return maxlen, minlen

if ver == 'hg19':
#fasta = '/home/sjpyo/sjpyo/data/v29/GRCh37.p13.genome.fa'
#	fasta = '/export/home/sjpyo/new/data/fasta/human/GRCh37.p13.genome.fa'
	fasta = '/home/sjpyo/new/data/fasta/human/hg19/hg19.fa'
	anno = '/home/sjpyo/new/data/anno/human/GENCODE/lift37/gencode.v34lift37.annotation.gtf' 
if ver == 'hg38':
	fasta = '/home/sjpyo/new/data/fasta/human/hg38/GRCh38.p10.genome.fa'
	anno = '/home/sjpyo/new/data/anno/human/GENCODE/GRCh38/gencode.v26.annotation.gtf'
#	anno = '/home/sjpyo/new/data/gencode.v29/gencode.v29lift37.annotation.gtf'
#elif ver == 'v19':
#	anno = '/export/home/sjpyo/data/anno/human/gencode/hg19/v19/gencode.v19.annotation.gtf'
#	anno = '/export/home/sjpyo/data/anno/human/gencode/hg19/v19/gencode.v19.annotation.chr22.gtf'
elif ver == 'mm10':
	fasta = '/export/home/sjpyo/new/data/fasta/mouse/GRCm38.p6.genome.fa'
	anno = '/export/home/sjpyo/new/data/anno/mouse/gencode.vM24.annotation.gtf'
if jobname == 'indexing':
#	if genome == 'hg19'
#		genomedir = '/home/sjpyo/new/data/star/'+ver+'/hg19_index/'
#	else:
	genomedir = '/home/sjpyo/new/data/star/'+ver+'/'
#	if os.path.exists(genomedir):
#		print 'You already generated genome indexes'
#		exit()
	if not os.path.exists(genomedir): os.makedirs(genomedir)
#	star = starPath + ' --runMode genomeGenerate --runThreadN ' + str(jobN)+ ' --genomeDir '+ genomedir + ' --genomeFastaFiles ' + fasta + ' --sjdbGTFfile ' + anno
	star = starPath + ' --runMode genomeGenerate --runThreadN 15 --genomeDir '+ genomedir + ' --genomeFastaFiles ' + fasta + ' --sjdbGTFfile ' + anno
#	print star + '\n'
#	commands.getoutput(star)
	queue = config.queue(server, 2)
	config.qsub_execute(jobname, queue, star, genomedir)

else: #jobname == 'mapping':
#	inputdir = sys.argv[6]+'/'
	inputdir = sys.argv[6]
	samples = os.listdir(inputdir)
	samples = filter(lambda x: not '.t' in x and not '.out' in x and not 'log' in x, samples)
	
##Define intron length
#	alignIntronMax, alignIntronMin = intronLen(anno, 0.1)
#	print 'alignIntronMax '  + ' : ' + str(alignIntronMax)
#	print 'alignIntronMin '  + ' : ' + str(alignIntronMin)
	genomedir = '/export/home/sjpyo/new/data/star/'+ver+'/'
	if not os.path.exists(genomedir) or len(os.listdir(genomedir)) < 1:
		print "Cannot find genome index directory!"
		exit()
	if jobname == 'mapping':
		print '--Start Alignment--'
		for i in xrange(len(samples)):
			sample = samples[i]
			if rpds == 'rpds':
				sinputdir = inputdir + sample + '/1.preprocess/rpds_fastqs/'
			else:
				if os.path.exists(inputdir+sample+'/1.preprocess/cutAdapt/'):
					sinputdir = inputdir+sample+'/1.preprocess/cutAdapt/'
				else:
					sinputdir = inputdir+sample+'/1.preprocess/sickle/'
#		sinputdir = inputdir + sample + '/'# + '/1.preprocess/sickle/'
			if os.path.exists(sinputdir):
				#print i, sample
				fastqs = os.listdir(sinputdir)
				if rpds == 'rpds':
					fastq1 = sinputdir + filter(lambda x: 'rpds_1.fastq' in x, fastqs)[0]
					fastq2 = sinputdir + filter(lambda x: 'rpds_2.fastq' in x, fastqs)[0]
				else:
					fastqs = [sinputdir + x for x in filter(lambda x: '.fastq' in x and not '_single' in x, fastqs)]
					print fastqs
					fastqs.sort()
					fastq1 = fastqs[0]; fastq2 = fastqs[1]
			
				if len(fastqs)>=2:
					fastqs = [sinputdir + x for x in fastqs]
#					print i, sample, fastq1,'\t',fastq2
					if rpds == 'rpds':
						soutputdir = inputdir+sample+'/2.mapping/rpds_star/'
					else:
						if os.path.exists(inputdir+sample+'/1.preprocess/cutAdapt/'):
							soutputdir = inputdir+sample+'/2.mapping/cutadpt_star/'
						else:
							soutputdir = inputdir+sample+'/2.mapping/star2/'
#							soutputdir = sys.argv[7]+sample+'/2.mapping/star/'
#					if not os.path.exists(soutputdir): os.makedirs(soutputdir)
#				star_comm = starPath + ' outReadsUnmapped Fastx --runThreadN 4 --outSAMtype BAM SortedByCoordinate --genomeDir '+ genomedir + ' --readFilesIn ' + ' '.join(fastqs) + ' --sjdbGTFfile ' + anno + ' --outFileNamePrefix ' + soutputdir + sample + '_'
#				queue=config.queue('biglab',i) ##
				#queue = 'biglab2'
#command line
#				star = starPath + ' --runThreadN 5 --runMode alignReads --genomeDir '+ genomedir + ' --readFilesIn ' + fastq1 + ' ' + fastq2 + ' --outFileNamePrefix ' + soutputdir + sample + '_ --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000'
					star = starPath+' --runMode alignReads --runThreadN '+config.runThread_star+' --sjdbOverhang '+config.sjdbOverhang+' --sjdbScore '+config.sjdbScore+' --outFilterMultimapScoreRange '+config.outFilterMultimapScoreRange+' --outFilterMultimapNmax '+config.outFilterMultimapNmax+' --outFilterMismatchNmax '+config.outFilterMismatchNmax+' --outFilterScoreMinOverLread '+config.outFilterScoreMinOverLread+' --outFilterMatchNminOverLread '+config.outFilterMatchNminOverLread+ ' --alignIntronMax '+config.alignIntronMax+' --alignMatesGapMax '+config.alignMatesGapMax+' --alignSJDBoverhangMin '+config.alignSJDBoverhangMin+' --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMstrandField intronMotif --genomeDir '+genomedir+' --readFilesIn '+fastq1 + ' ' + fastq2 +' --outFileNamePrefix '+soutputdir + sample + '_'
				queue = config.queue(server, i)
				config.qsub_time(jobname, jobN, 300)
#run STAR
				config.qsub_execute(jobname, queue, star, soutputdir)
			print star + '\n'
#			commands.getoutput(star)
		while commands.getoutput('qstat').count(jobname) > 0: time.sleep(10)
		print '--End Alignment--\n'
	if jobname == 'mapping' or jobname == 'bamindex':
		print '--Start BAM Index--'
#	for i in xrange(len(samples)):
#	sample = samples[i]
#	bamfile = soutputdir + filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(sinputdir))[0]
		indexing = bamIndex + ' ' + inputdir + ' ' + rpds + ' ' + server + ' ' + str(jobN)
		print indexing
		commands.getoutput(indexing)
		while commands.getoutput('qstat').count('index') > 0: time.sleep(10)
		print '--End BAM Index--\n'
	if jobname == 'mapping' or jobname == 'bamindex' or jobname == 'strand':
		print '--Determine Strandness--'
		if rpds == 'rpds':
			logfile = open(inputdir + 'strandness_RPDs_RSeQC.txt', 'w')
		else:
			logfile = open(inputdir + 'strandness_RSeQC.txt', 'w')
		logfile.write('sample\tfraction\tstrandness\n')
		unstranded = []
		for i in xrange(len(samples)):
			sample = samples[i]
			print sample
			if ver == 'hg19': ref = '/home/sjpyo/new/data/anno/human/GENCODE/lift37/hg19_GencodeCompV19.bed'
			else: ref = '/home/sjpyo/new/data/anno/mouse/mm10_Gencode_VM18.bed'
			if rpds == 'rpds': sinputdir = inputdir + sample + '/2.mapping/rpds_star/'
			else: sinputdir = inputdir + sample + '/2.mapping/star/'
			bamfile = sinputdir + filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(sinputdir))[0]
			rseqc = rseqcPath + ' -r ' + ref + ' -i ' + bamfile
			print rseqc
			output = commands.getoutput(rseqc).split('\n')
#			print output
#			print output[-2].split(':')[-1].strip()
			fraction = float(output[-2].split(':')[-1].strip())

			logfile.write(sample+'\t'+str(fraction)+'\t')
			if fraction > 0.8:
				logfile.write('fr-secondstrand\n')
			elif fraction < 0.2:
				logfile.write('fr-firststrand\n')
			else:
				unstranded.append(sample)
				logfile.write('Unstranded\n')
		logfile.close()
