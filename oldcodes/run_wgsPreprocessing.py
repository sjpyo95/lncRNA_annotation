import sys, os
import argparse
import commands
##Run inside the node
##ssh biglab#

fastqc = '/share/apps/programs/fastqc/0.11.8/fastqc'
sickle = '/share/apps/programs/sickle/1.33/sickle'
java = '/share/apps/programs/java/jdk1.8.0_231/bin/java'
picard = '/home/sjpyo/new/packages/picard/build/libs/picard.jar'

def fastqc(f1, f2, targetDir, readType):
#a command line
	fastqFiles = [f1, f2]	
	if readType == 'paired'or readType == 'single':
		outputDir = targetDir + '1.preprocess/fastqc/'
		if not os.path.exists(outputDir): os.makedirs(outputDir)
#		fastqFiles = map(lambda x: targetDir + x, fastqFiles)
	elif readType == 'trimmed_paired':
		outputDir = targetDir + '../fastqc-sickle/'
		if not os.path.exists(outputDir): os.makedirs(outputDir)
#		fastqFiles = map(lambda x: targetDir + x, fastqFiles)

	fastqcCmd = '/share/apps/programs/fastqc/0.11.8/fastqc -o ' + outputDir + ' --extract -f fastq -t 4 -q '
	if readType == 'paired' or readType == 'trimmed_paired':    fastqcCmd += fastqFiles[0] + ' ' + fastqFiles[1]
	elif readType == 'single':  fastqcCmd += fastqFiles[0]
	jobNameBase = False
	print '\n' + fastqcCmd + '\n'
#run the command line
	commands.getoutput(fastqcCmd)

def sickle(f1, f2, targetDir, prefix):
	fastqFiles = [f1,f2]
	if prefix:
		name = prefix
	if not prefix:
		name = fastqFiles[0].split('.fa')[0].split('R1')
#file sorting
	outputDir = targetDir + '1.preprocess/sickle/'
	if not os.path.exists(outputDir): os.makedirs(outputDir)
#	fastqFiles = [targetDir + x for x in fastqFiles]
	fastqFiles.sort()
#sickle defaults
	quality = int(20); length = int(20)
	quality_type = 'sanger'
	readType = 'paired'
#a command line
	if readType == 'paired': sickle = '/share/apps/programs/sickle/1.33/sickle pe -f ' + fastqFiles[0] + ' -r ' + fastqFiles[1]
	elif readType == 'single': sickle = '/share/apps/programs/sickle/1.33/sickle se -f ' + fastqFiles[0]
	sickle += ' -t ' + quality_type + ' -q ' + str(quality) + ' -l ' + str(length)
	if readType == 'paired':
		sickle += ' -o ' + outputDir + name + '_R1_trimmed.fastq'
		sickle += ' -p ' + outputDir + name + '_R2_trimmed.fastq'
		sickle += ' -s ' + outputDir + name + '_single_trimmed.fastq'
	elif readType == 'single':
		sickle += ' -o ' + outputDir + fastqFiles[0].split('.fa')[0] + '_trimmed.fastq'
	jobNameBase = False
	print '\n' + sickle + '\n'
	commands.getoutput(sickle)

def fastqToSam(f1, f2, targetDir, prefix):
	outputDir = targetDir+'1.preprocess/uBam/'
	if not os.path.exists(outputDir): os.makedirs(outputDir)
	tmpdir = targetDir + '/tmp/'
	if prefix:
		name = prefix
	if not prefix:
		name = ' '.join(targetDir.split('/')).split()[-1]
	fastqToSam_cmm = java + ' -Xmx2g -Djava.io.tmpdir='+tmpdir+' -jar '+picard + ' FastqToSam F1='+f1+' F2='+f2+' O='+outputDir+name+'_unaligned_read_pair.bam SM=sample_'+name+' RG=rg_'+name+' TMP_DIR='+tmpdir
	print '\n'+fastqToSam_cmm+'\n'
	commands.getoutput(fastqToSam_cmm)

def markAdapters(uBAM, targetDir, prefix):
	outputDir = targetDir + '/1.preprocess/adapterMark/'
	if not os.path.exists(outputDir): os.makedirs(outputDir)
	tmpdir = targetDir+'/tmp/'
	if prefix:
		 name = prefix
	if not prefix:
		name = ' '.join(targetDir.split('/')).split()[-1]
	markAdapters = java + ' -jar '+picard + ' MarkIlluminaAdapters I='+ uBAM+ ' O= ' + outputDir+ name+ '_markAdapters.bam M='+ outputDir+name+'_markAdapters._metrics.txt'
	print '\n'+markAdapters+'\n'
	commands.getoutput(markAdapters)

def samToFastq(markBam, targetDir, prefix):
	outputDir = targetDir + '/1.preprocess/samToFastq/'
	if not os.path.exists(outputDir): os.makedirs(outputDir)
	tmpdir = targetDir+'/tmp/'
	if prefix:
		name = prefix
	if not prefix:
		name = ' '.join(targetDir.split('/')).split()[-1]
	samToFastq = java + ' -jar '+picard + ' SamToFastq I='+markBam+' FASTQ='+ outputDir + name +'_interleaved.fastq CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true'
	print '\n'+samToFastq+'\n'
	commands.getoutput(samToFastq)

def __main__(args):
	fastqc(args.f1, args.f2, args.outputdir, 'paired')
	print 'FastQC of raw data finished.'
	sickle(args.f1, args.f2, args.outputdir, args.prefix)
	print 'Sickle finished.'
	trimF1 = args.outputdir + '/1.preprocess/sickle/' + filter(lambda x: 'trimmed.fastq' in x and 'R1' in x and not 'single' in x, os.listdir(args.outputdir + '/1.preprocess/sickle/'))[0]
	trimF2 = args.outputdir + '/1.preprocess/sickle/' + filter(lambda x: 'trimmed.fastq' in x and 'R2' in x and not 'single' in x, os.listdir(args.outputdir+ '/1.preprocess/sickle/'))[0]
	fastqc(trimF1, trimF2, args.outputdir, 'trimmed_paired')
	print 'trimmed reads FastQC finished.'

	tmpdir = args.outputdir+'/tmp/'
	if not os.path.exists(tmpdir): os.makedirs(tmpdir)

	fastqToSam(trimF1, trimF2, args.outputdir, args.prefix)
	print 'Unaligned BAM (uBAM) file is produced.'
	ubam = args.outputdir+'1.preprocess/uBam/' + filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(args.outputdir+'1.preprocess/uBam/'))[0]
	markAdapters(ubam, args.outputdir, args.prefix)
	print 'Mark Adapters by picard MarkIlluminaAdapters'
	markBam = args.outputdir + '1.preprocess/uBam/' + filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(args.outputdir + '1.preprocess/uBam/'))[0]
	samToFastq(markBam, args.outputdir, args.prefix)
	print 'Bam file to single interleaved FASTQ file'

if __name__=='__main__':
	parser = argparse.ArgumentParser(description = 'WGS fastq files preprocessing')
	parser.add_argument('-f1', help = 'foward reads fastq file', required = True)
	parser.add_argument('-f2', help = 'reverse reads fastq file', required = True)
	parser.add_argument('-o', '--outputdir', help = 'the directory where preprocessed files will be produced', required = True)
	parser.add_argument('-p', '--prefix', help = 'header name for all producing files', required = False)
#	parser.add_argument('-s', '--step', help = 'The step that you want to start (e.g. fastqc, sickle, fastqc-sickle, fastqToSam, etc.)', required = False)
	args = parser.parse_args()
	__main__(args)
