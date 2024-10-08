import sys, os

inputdir = sys.argv[1] + '/'
outputdir = sys.argv[2] + '/'
if not os.path.exists(outputdir): os.makedirs(outputdir)
samples = os.listdir(inputdir)
samples = filter(lambda x: '.t' not in x and '.out' not in x and '.log' not in x, samples)
outputfile = open(outputdir + 'star_quality_table.txt', 'w')
outputfile.write('Library_type\tUniquely_mapped_reads\tMultiple_loci_reads\tMismatched_reads(unmapped)\tToo_short_reads(unmapped)\tOther_unmapped_reads\n')

def parse_starLog(infile):
	infile = open(infile, 'r')
	lines = infile.readlines()
	infile.close()
	for line in lines:
		line = line.strip().split('|')
		title = line[0].strip(); value = line[-1].strip().split('%')[0]

		if title == 'Uniquely mapped reads %':
			uq_reads = value
		elif title == '% of reads mapped to multiple loci':
			multi_reads = value
		elif title == '% of reads unmapped: too many mismatches':
			un_mis_reads = value
		elif title == '% of reads unmapped: too short':
			un_short_reads = value
		elif title == '% of reads unmapped: other':
			un_other_reads = value
	return [uq_reads, multi_reads, un_mis_reads, un_short_reads, un_other_reads]

for i in xrange(len(samples)):
	sample = samples[i]

	sinputdir = inputdir + sample + '/2.mapping/rpds_star/'
	print sinputdir
	logFile = filter(lambda x: 'Log.final.out' in x, os.listdir(sinputdir))[0]
	logInfo = parse_starLog(sinputdir+logFile)
	outputfile.write(sample + '\t')
	for i in xrange(len(logInfo)):
		outputfile.write(logInfo[i] + '\t')
	outputfile.write('\n')
outputfile.close()
