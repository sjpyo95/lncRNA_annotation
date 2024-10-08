import sys, os
import commands
samtools = '/share/apps/programs/samtools/1.7/bin/samtools'
faTonib = '/home/sjpyo/new/packages/blat/faToNib'
inputdir = sys.argv[1]
fasta = filter(lambda x: '.fa' in x and not '.fai' in x, os.listdir(inputdir))
for i in xrange(len(fasta)):
	fa = fasta[i]
	makefai = samtools + ' faidx ' + inputdir+fa
	makenib = faTonib + ' ' + inputdir+fa+ ' ' + inputdir+fa.split('.fa')[0]+'.nib'
	print '\n'+makefai
	commands.getoutput(makefai)
	print '\n'+makenib
	commands.getoutput(makenib)
