import sys, os
import commands

inputdir = sys.argv[1]
samples = filter(lambda x: not '.txt' in x, os.listdir(inputdir))

for sample in samples:
	sinputdir = inputdir + sample + '/1.preprocess/rpds_fastqs/'
	fastq1files = sinputdir + sample + '_chr*_1.fastq'
	fastq2files = sinputdir + sample + '_chr*_2.fastq'
	print 'cat ' + fastq1files + ' > ' + sinputdir + sample + '.rpds_1.fastq'
	print commands.getoutput('cat ' + fastq1files + ' > ' + sinputdir + sample + '.rpds_1.fastq')
	print 'cat ' + fastq2files + ' > ' + sinputdir + sample + '.rpds_2.fastq'
	print commands.getoutput('cat ' + fastq2files + ' > ' + sinputdir + sample + '.rpds_2.fastq')
