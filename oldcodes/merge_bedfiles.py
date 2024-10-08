import sys, os
inputdir = sys.argv[1]
outfile = sys.argv[2]
#if outputdir == 'same':
#	outputdir = inputdir

def readingBed(bedfiles):
	junctions = dict()
	for i in xrange(len(bedfiles)):
		bedfile = bedfiles[i]
		print bedfile
		bed = open(bedfile, 'r')
		lines = bed.readlines(); bed.close()
		for i in xrange(len(lines)):
			line = lines[i].strip().split('\t')
			chrom = line[0]
			rstart = int(line[1]); rend = int(line[2])
			count = int(line[4])
			sense = line[5]
			if sense == '?': continue
			block = line[-2]
#			if sense != '+' and sense != '-': continue
			start = rstart + int(block.split(',')[0])+1
			end = rend - int(block.split(',')[1])
			if not junctions.has_key(chrom):
				junctions[chrom] = dict()
			intron = [chrom, str(start), str(end), sense]
			if not junctions[chrom].has_key(str(intron)):
				junctions[chrom][str(intron)] = [intron, count]
			else:
#				print intron
				junctions[chrom][str(intron)][1] += count
	return junctions

def writeBed(junctions, outputfile):
	output = open(outputfile, 'w')
#	junNum = 0
	for chrom in junctions.keys():
		jdic = junctions[chrom]
		for intron in jdic.keys():
#			junNum += 1
			info = jdic[intron][0]
			count = jdic[intron][1]
#			juncName = 'JUNC' + str(junNum)
			output.write('\t'.join(info[:-1])+'\t' + str(count)+'\t'+info[-1]+'\n')
	output.close()

bedfiles = [inputdir + x for x in filter(lambda x: '.bed' in x, os.listdir(inputdir))]
junctions = readingBed(bedfiles)
#outfile = outputdir+bedfiles[0].split('/')[-1][:-4].split('_rep')[0]+'_merged.bed'
writeBed(junctions, outfile)
