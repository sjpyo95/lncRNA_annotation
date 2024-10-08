import sys, os
import commands
import config as cf

stringtiePath = '/share/apps/programs/stringtie/stringtie-1.3.4a.Linux_x86_64/stringtie'
#refgenPath = '/export/home/sjpyo/data/anno/human/gencode/hg19/gencode.v19.annotation.gtf'

bamfile = sys.argv[1]
strand = sys.argv[2]	#--fr or --rf
refgenPath = sys.argv[3]
outputdir = sys.argv[4]
if not os.path.exists(outputdir): os.makedirs(outputdir)
sample = bamfile.split('/')[-1][:-4]
if 'Aligned.sortedByCoord.out' in sample:
	sample = sample.split('Aligned.sortedByCoord.out')[0]
stringtie = stringtiePath+' '+ bamfile+' -p '+ cf.runThread_string + ' -v -f '+cf.minIsoformFraction+' -m '+cf.minTranscriptLength+' -a '+cf.minAnchorLength+' -j '+cf.minJunctionCoverage+' -c '+cf.minReadCoverage+' -M '+cf.maxMultiCoverage+' -g '+cf.gapBetweenBundle+' -s '+ cf.minSingleExonCoverage + ' -G ' + refgenPath + ' '+ strand + ' -o '+outputdir+sample+'stringtie.gtf'

print stringtie + '\n'
commands.getoutput(stringtie)
