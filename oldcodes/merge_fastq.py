import sys, os

inputdir = sys.argv[1]
samples = filter(lambda x: not '.txt' in x and not '.out' in x, os.listdir(inputdir))
#print samples
def readFastq(fastq):
	infile = open(fastq, 'r')
	lines = infile.readlines(); infile.close()
	print len(lines)
	reads = []
	for i in range(0,len(lines), 4):
		read = lines[i:i+4]
		info = read[0]; seq = read[1]; info2 = read[2]; qual = read[3]
		reads.append(read)
	return reads
	
def catFastq(reads):
	newRead = []
	for i in xrange(len(reads)):
		read = reads[i]
		info = read[0].strip().split(' '); seq = read[1].strip(); info2 = read[2].strip().split(' '); qual = read[3].strip()

		ninfo = ' '.join([info[0].split('.')[0]+'.'+str(i+1), str(i+1), info[2]])
		ninfo2 = ' '.join([info2[0].split('.')[0]+'.'+str(i+1), str(i+1), info2[2]])

		newRead.append('\n'.join([ninfo, seq, ninfo2, qual]))
#		print len(newRead)
	return newRead

for i in xrange(len(samples)):
	sample = samples[i]
	sinputdir = inputdir + sample + '/'
	fastqs = filter(lambda x: 'fastq' in x and 'SRR' in x, os.listdir(sinputdir))
	fastq1s = sorted(filter(lambda x: '_1' in x, fastqs)); fastq2s = sorted(filter(lambda x: '_2' in x, fastqs))
	print fastq1s
	print fastq2s
	readslist1 = []
	readslist2 = []
	for fq in fastq1s:
		print fq
		reads = readFastq(sinputdir + fq)
#		print len(reads)
		readslist1 += reads
	newReads1 = catFastq(readslist1)
	for fq in fastq2s:
		print fq
		reads = readFastq(sinputdir + fq)
		readslist2 += reads
	newReads2 = catFastq(readslist2)
	nfastq1 = open(sinputdir+sample+'_1.fastq', 'w')
	nfastq2 = open(sinputdir+sample+'_2.fastq','w')
	nfastq1.write('\n'.join(newReads1))
	nfastq2.write('\n'.join(newReads2))
	nfastq1.close(); nfastq2.close()
	

