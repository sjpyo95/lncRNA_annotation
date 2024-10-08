import sys, os
import module as mdl
import config as cf
fafile = sys.argv[1]
outfile = open(sys.argv[2], 'w')
specie = sys.argv[3]
if specie == 'human':
	chrs = cf.hm_chrs
elif specie == 'mouse':
	chrs = cf.ms_chrs
seqdic = mdl.getSeq_dict(fafile)
for chr in chrs:
	seq = seqdic[chr]
	outfile.write(chr+'\t'+str(len(seq))+'\n')
outfile.close()

