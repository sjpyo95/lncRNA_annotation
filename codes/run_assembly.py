import sys, os
import module as mdl
import config as cf
stringtiePath = '/share/apps/programs/stringtie/stringtie-1.3.4a.Linux_x86_64/stringtie'
refgenePath = '/home/sjpyo/new/data/anno/human/GENCODE/lift37/gencode.v34lift37.annotation.gtf'
inputdir = sys.argv[1]
outputdir = sys.argv[2]
if not os.path.exists(outputdir): os.makedirs(outputdir)

library = cf.library('stringtie', sys.argv[3]) #unstranded, fr-firststrand, fr-secondstrand
jobN = int(sys.argv[4])

bamfiles = [inputdir + s for s in filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(inputdir))]
print 'Total samples: ', len(bamfiles)
print 'Input directory: ' + inputdir
print '<< Parameters >>'
print ' -p '+ cf.runThread_string
print ' -f '+cf.minIsoformFraction
print ' -m '+cf.minTranscriptLength
print ' -a '+cf.minAnchorLength
print ' -j '+cf.minJunctionCoverage
print ' -j '+cf.minJunctionCoverage
print ' -c'+cf.minReadCoverage
print ' -M '+cf.maxMultiCoverage
print ' -g '+cf.gapBetweenBundle
print ' -s '+ cf.minSingleExonCoverage
print ' -G ' + refgenePath
print library
print ' -o '+outputdir+'SAMPLE_NAME.stringtie.gtf'

for i in range(len(bamfiles)):
	bamfile = bamfiles[i]
	sample = bamfile.split('/')[-1].split('_Aligned.sortedByCoord.out')[0]

	stringtie_CM = stringtiePath + ' '+ bamfile+' -p '+ cf.runThread_string + ' -f '+cf.minIsoformFraction+' -m '+cf.minTranscriptLength+' -a '+cf.minAnchorLength+' -j '+cf.minJunctionCoverage+' -c '+cf.minReadCoverage+' -M '+cf.maxMultiCoverage+' -g '+cf.gapBetweenBundle+' -s '+ cf.minSingleExonCoverage + ' -G ' + refgenePath + ' '+ library + ' -o '+outputdir+sample+'.stringtie.gtf'
	queue = cf.queue('biglab', i)
	cf.qsub_time('PYO_job', jobN, 60)
	cf.qsub_execute('PYO_job', queue, stringtie_CM, outputdir+sample)
	

