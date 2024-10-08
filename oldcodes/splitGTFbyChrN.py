import sys, os
from subprocess import Popen, PIPE, STDOUT
import commands, time
#inputdir = sys.argv[1] + '/'
inputFile = sys.argv[1]
specie = sys.argv[2]
outputdir = sys.argv[3] + '/assembly_chrN/'#+ '/'+assemblyFile[:-4]+'_chrN/'#inputdir+assemblyFile[:-4] +'_chrN/'
#tissue = sys.argv[2]
#queue = sys.argv[3]
#jobN = int(sys.argv[4])
#ingtfFile = sys.argv[1] + '/'
#gtfFile = open(ingtfFile, 'r')
#outputdir = ingtfFile.split('/')[:-1]
#outDir = '/'.join(outputdir) + '/chrN/'
#gtflines = gtfFile.readlines();gtfFile.close()
#create gtf files of each chromosomes

#def qsub_execute(job, queue, command, outputdir):
#	qsub = 'qsub -N '+job+' -j oe -lvnode='+queue+' -o '+outputdir+'qsub.log -- '+command
#	print qsub
#	p=Popen(qsub, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
#	p.wait()

def splitGtf(hm_chrs, inputFile, outDir):
	gtfFile = open(inputFile, 'r')
	gtflines = gtfFile.readlines();gtfFile.close()
	for i in xrange(len(hm_chrs)):
		chrom = hm_chrs[i]
		outfile = open(outDir + '/'+inputFile.split('/')[-1][:-4] + '_'+chrom + '.gtf', 'w')
		for line in gtflines:
			chrN = line.strip().split('\t')[0]
			if chrN == chrom:
				outfile.write(line)
	outfile.close()
if specie == 'human': chrs = map(lambda x: 'chr' + x, map(str, range(1, 23)) + ['X', 'Y'])
elif specie == 'mouse': chrs = map(lambda x: 'chr' + x, map(str, range(1, 20)) + ['X', 'Y'])
#assemblyFiles = os.listdir(inputdir)
#assemblyFiles = filter(lambda x: '.gtf' in x, assemblyFiles)
#for i in xrange(len(assemblyFiles)):
assemblyFile = inputFile.split('/')[-1]
if not os.path.exists(outputdir) : os.makedirs(outputdir)
splitGtf(chrs, inputFile, outputdir)

#	for i in xrange(len(hm_chrs)):
#		chrN = hm_chrs[i]
#			command = 'grep -w ' + chrN + ' ' + sinputdir + assemblyFile + ' > ' + soutputdir + assemblyFile[:-4] + '.' + chrN + '.gtf'
#			print '\n' + command + '\n'
   		#while commands.getoutput('qstat').count(tissue)>jobN: time.sleep(30)
#			commands.getoutput(command)
#				splitGtf(hm_chrs, sinputdir+assembly,soutputdir)
