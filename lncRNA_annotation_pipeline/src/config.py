#!/usr/bin/python
import sys, os
from subprocess import Popen, PIPE, STDOUT
import commands, time
import re
# queue name on biglab, cngc or kt(snu) servers #
def queue(server, i):
	if server=='biglab':
		if i%2==1: queue='biglab6'
#		elif i%24==2: queue='biglab4'
#		elif i%4==3: queue='biglab4'
#		elif i%7==4: queue='biglab4'
#		elif i%7==5: queue='biglab6'
#		elif i%7==6: queue='biglab6'
		else: queue='biglab6'
#		queue='biglab2'
#	elif server=='cngc':
#		queue='-lncpus=4'
#		queue='-lncpus=4 -lmem=10gb'
#	elif server=='kt':
#		queue='BIG_lab.q'
	return queue
hm_chrs = map(lambda x: 'chr' + x, map(str, range(1,23))+['X', 'Y', 'M'])
ms_chrs = map(lambda x: 'chr' + x, map(str, range(1, 20)) + ['X', 'Y'])
def qsub_execute(job, queue, command, logdir):
#	if server=='biglab':
	qsub = 'qsub -N '+job+' -lvnode='+queue+' -j oe -o '+logdir+'qsub.log -- '
#	elif server=='cngc':
#		qsub = 'qsub -N '+job+' '+queue+' -j oe -o '+logdir+'qsub.log -- '
#	elif server=='kt':
#		qsub = 'qsub -N '+job+' -q '+queue+' -b y -o '+logdir+'qsub.out.log -e '+logdir+'qsub.error.log '
	qsub+=command
	print '\n'+qsub+'\n'
	p=Popen(qsub, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()

def qsub_time(job, jobN, t):
	while commands.getoutput('qstat').count(job)>jobN: time.sleep(t)

# Parameters for RNA-seq analysis #

# library setting
def library(program, lib_type='unstranded'):
	if lib_type=='unstranded':
		if program=='stringtie':
			return ' '
		elif program=='featureCounts':
			return '0 '
	elif lib_type=='fr-firststrand':
		if program=='stringtie':
			return ' --rf '
		elif program=='featureCounts':
			return '1 '
	elif lib_type=='fr-secondstrand':
		if program=='stringtie':
			return ' --fr '
		elif program=='featureCounts':
			return '2 '

# run parameters
runThread_star='8'

##STAR mapping
# recommended parameters in ICGC

# splice junctions database
sjdbOverhang='100' #default: 100
sjdbScore='2' #default: 2

# output filtering
outFilterMultimapScoreRange='1' #default: 1
outFilterMultimapNmax='20' #default: 10
outFilterMismatchNmax='10' #default: 10
outFilterScoreMinOverLread='0.33' #default: 0.66
outFilterMatchNminOverLread='0.33' #default: 0.66

# alignments and seeding
alignIntronMax='500000' #default: 0
alignMatesGapMax='1000000' #default: 0
alignSJDBoverhangMin='1' #default: 3

##StringTie assembly

# run parameters
runThread_string='6'

# 
minIsoformFraction='0.05' #default: 0.1
minTranscriptLength='100' #default: 200
minAnchorLength='5' #default: 10
minJunctionCoverage='1' #default: 1
minReadCoverage='1' #default: 2.5
maxMultiCoverage='1.0' #default: 1.0
minSingleExonCoverage = '4.75' #default: 4.75

# 
gapBetweenBundle='50' #default: 50

##featureCounts counting reads

# annotation
formatAnno='GTF' #default: GTF
def feature(specify='default'):
	if specify=='default':
		return ' -t exon -g gene_id '
	else:
		return ' '+specify+' '

# reads and features
minOverlapLength='45' #default: 1
minOverlapFrac='0' #default: 0
def multiLoci(reads=False):
	if reads:
		return ' -M --fraction '
	else:
		return ''

# read filtering
minMapQual='0' #default: 0
def primary(alignments=False):
	if alignments:
		return ' --primary '
	else:
		return ' '

# paired-end reads
def paired(read_type='paired'):
	if read_type=='paired':
		return '-p -B -C '
	else:
		return ''

def tacoParse(gtfFile):
	gtf = open(gtfFile, 'r')
	lines = gtf.readlines(); gtf.close()
	for line in lines:
		line = line.strip()
		kind = line.split('\t')[2]
		trxid = re.search('.*transcript_id "(\S+)";', line).group(1)
		if kind == 'transcript':
			trxD[trxid] = [line]
		elif kind == 'exon':
			trxD[trxid].append(line)
	outputfile = open(sinputdir + sample + '_merged_transcripts.gtf', 'w')
	for trxid in trxD.keys():
		trx = trxD[trxid][0]; exon = trxD[trxid][1:]
		outputfile.write('\n'.join(trxD[trxid])+'\n')
	outputfile.close()	
