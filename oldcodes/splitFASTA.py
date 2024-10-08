import sys, os
import config as cf

infile = open(sys.argv[1], 'r')
specie = sys.argv[2]
outputdir = sys.argv[3]
if not os.path.exists(outputdir): os.makedirs(outputdir)
lines = infile.read(); infile.close()
lines = lines.split('>')
splitfa = dict()
for line in lines:
	chrom = line.split('\n')[0].split(' ')[0]
	seq = '\n'.join(line.split('\n')[1:])
	splitfa[chrom] = seq
if specie == 'human': chroms = cf.hm_chrs
elif specie == 'mouse': chroms = cf.ms_chrs
for chr in chroms:
	outfile = open(outputdir + chr+'.fa', 'w')
	outfile.write('>'+chr+'\n'+splitfa[chr])
	outfile.close()

