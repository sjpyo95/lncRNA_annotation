import sys, os
import argparse
import commands, time

def run_rpds(inputdir , specie, stranded, jobN):
	rpdsPath = '/home/sjpyo/new/codes/map_assembly_annotation_pipeline/run_rpds_v3.py'
	runRPDs = 'python '+rpdsPath+' '+inputdir+' '+specie+' '+stranded+' '+inputdir+' biglab '+jobN
	print commands.getoutput(runRPDs)
def mergeRPDS_fa(sample, inputdir):
	catfastq1 = 'cat ' + inputdir + '*_chr*_1.fastq > '+ inputdir + sample+'.rpds_1.fastq'
	catfastq2 = 'cat ' + inputdir + '*_chr*_2.fastq > '+ inputdir + sample+'.rpds_2.fastq'
	print commands.getoutput(catfastq1)
	print commands.getoutput(catfastq2)
	
def run_star(inputdir, sp, rpds, jobN):
	starPath = '/home/sjpyo/new/codes/run_star_v2.py'
	if sp == 'human': specie = 'hg19'
	elif sp == 'mouse': specie = 'mm10'
	runSTAR = 'python '+starPath+' mapping ' + specie + ' ' + rpds + ' biglab ' + jobN+' '+inputdir
	print runSTAR
	print commands.getoutput(runSTAR)
	


def check_strand_specific(rseqc_output):
	infile = open(rseqc_output, 'r')
	lines = infile.readlines(); infile.close()
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if not line.startswith('sample'):
			tmp = line.split('\t')
			sample = tmp[0]; frac = float(tmp[1]); strand = tmp[2]
			if frac < 0.2 or strand == 'fr-firststrand':
				st = 'fr-firststrand'
			elif frac > 0.8 or strand == 'fr-secondstrand':
				st = 'fr-secondstrand'
			else:
				st = 'unstranded'
	return st

def run_stringtie(inputdir, rpds, st, anno, jobN):
	stringtiePath = '/home/sjpyo/new/codes/run_stringtie_v2.py'
	if st == 'fr-firststrand': stS = '--rf'
	elif st == 'fr-secondstrand': stS = '--fr'
	runStringtie = 'python '+stringtiePath+ ' '+ inputdir + ' ' + rpds + ' ' + stS + ' ' + anno + ' biglab '+jobN+' '+inputdir
	print commands.getoutput(runStringtie)

def run_stringtieMerge(inputdir, outputdir, jobN, cluster, specie):
	stringmergePath = '/home/sjpyo/new/codes/run_stringtieMerge.py'
	stringmerge = '/usr/bin/python '+stringmergePath + ' '+ inputdir + ' ' + outputdir + '/assembly/ biglab '+ jobN+ ' '+cluster
	print commands.getoutput(stringmerge)

def run_regtools(inputdir, outputdir, st, jobN):
	regtoolsPath = '/home/sjpyo/new/codes/run_regtools.py'
	mergebedPath = '/home/sjpyo/new/codes/merge_bedfiles.py'
	splitBedPath = '/home/sjpyo/new/codes/splitBEDinChrN.py'
	if st == 'fr-firststrand': stN = '1'
	elif st == 'fr-secondstrand': stN = '2'
	routputdir = outputdir + '/regtools/'
	regtools = '/usr/bin/python '+regtoolsPath + ' ' + inputdir + ' '+ routputdir+ ' ' + stN + ' biglab '+ jobN
	print regtools
	print commands.getoutput(regtools)
	while commands.getoutput('qstat').count('regtools') > 0: time.sleep(5)
	mergefile = routputdir+inputdir.split('/')[-2]+'.merge_junction.bed'
	mergeBed = '/usr/bin/python '+ mergebedPath + ' ' + routputdir+' '+ mergefile
	print commands.getoutput(mergeBed)
	splitBed = '/usr/bin/python '+ splitBedPath + ' ' + mergefile + ' ' + routputdir+'/chrN/'
	print commands.getoutput(splitBed)
	beddir = routputdir+'/chrN/'
	return beddir

def run_exj_update(chrdir, cafeoutdir, beddir, specie, jobN):
	outputdir = cafeoutdir + '/cafe/'
	if specie == 'human': sp = 'hg19'
	elif specie == 'mouse': sp = 'mm10'
	exjPath = '/home/sjpyo/new/codes/run_exj_revised2.py'
	exj ='/usr/bin/python '+exjPath+' '+ chrdir + ' ' + beddir + ' ' + outputdir + ' ' + sp + ' read 2 biglab '+ jobN
	print commands.getoutput(exj)
	return outputdir

def run_tss_cps_update(inputdir, specie, jobN):
	if specie == 'human': sp = 'hg19'
	elif specie == 'mouse': sp = 'mm10'
	endsupPath = '/home/sjpyo/new/codes/run_tss_cps_filter.py'
	endsup = '/usr/bin/python '+endsupPath+ ' '+ inputdir + ' ' + inputdir + ' ' + sp + ' tss near biglab ' +jobN
	print commands.getoutput(endsup)

def run_lncRNA_filtration(inputdir, outputdir, specie, jobN):
	lncFiltPath = '/home/sjpyo/new/codes/map_assembly_annotation_pipeline/run_lncRNAfilter.py'
	lncFilt = 'python '+ lncFiltPath+' '+inputdir+' '+outputdir+' '+specie+' biglab '+jobN
	print lncFilt
	commands.getoutput(lncFilt)

def run_CPC(specie, inputgtf, outputdir, txnum, jobN):
	cpcPath = '/home/sjpyo/new/codes/lncRNA_annotation/CPC.step1_v2.py'
	cpc = 'python '+cpcPath+' '+specie+' '+inputgtf+' '+outputdir+' '+str(txnum)+' biglab '+jobN
	print cpc
	print commands.getoutput(cpc)
#	faDir = outputdir+'/fasta/'
#	return faDir

def get_TrxNum(gtf):
	commL = '''awk -F '\t' '$3=="transcript"' ''' + gtf + ' | wc -l'
	return commands.getoutput(commL)

def run_CPAT(faDir, outputdir, jobN, specie):
	cpatPath = '/home/sjpyo/new/codes/lncRNA_annotation/run_cpat.py'
	cpat = 'python '+cpatPath+' '+faDir+' '+outputdir+' '+jobN+' '+specie
	print commands.getoutput(cpat)

def get_CPC_CPAT_result(cpcdir, cpatdir, pTrx, outputdir, sp):
##Get CPC result
	cpcresults = [cpcdir + '/' + x for x in filter(lambda x: 'result' in x and '.txt' in x, os.listdir(cpcdir))]
	print commands.getoutput('cat ' + ' '.join(cpcresults) + ' > ' + cpcdir + '/putative_transcripts.result.all.txt')
	cpc_allresults = cpcdir + '/putative_transcripts.result.all.txt'
	getcpcPath = '/home/sjpyo/new/codes/map_assembly_annotation_pipeline/getcpcLnc.py'
	getcpc = 'python ' + getcpcPath + ' ' + cpc_allresults + ' ' + pTrx + ' ' + outputdir
	print getcpc
	print commands.getoutput(getcpc)
##Get CPAT result
	pTrx2 = outputdir + '/' + filter(lambda x: '.gtf' in x and 'cpc' in x, os.listdir(outputdir))[0]
	print commands.getoutput('cat ' + cpatdir+'/*CPAT_result.txt > ' + cpatdir + '/putative_transcripts.all.CPAT_result.txt')
	cpat_allresults = cpatdir + '/putative_transcripts.all.CPAT_result.txt'
	getcpatPath = '/home/sjpyo/new/codes/map_assembly_annotation_pipeline/getcpatLnc.py'
	getcpat = 'python ' + getcpatPath + ' ' + cpat_allresults + ' ' + pTrx2 + ' ' + outputdir+' '+sp
	print getcpat
	print commands.getoutput(getcpat)

def novel_filter(novel, anno, outfile, logfile):
	filterPath = '/home/sjpyo/new/codes/map_assembly_annotation_pipeline/novel_filter.py'
	novelfilt = 'python ' + filterPath + ' ' + novel + ' ' + anno + ' ' + outfile + ' ' + logfile
	print novelfilt
	print commands.getoutput(novelfilt)

def __main__(args):
	jobN = args.jobN
	step = args.step
	if not step:
		step = 'star'
##STEP1 : STAR Alignment
	if step == 'star':
		print '-- Run STAR Alignment --'
		run_star(args.inputdir, args.specie, 'no', jobN)
		while commands.getoutput('qstat').count('mapping') > 0: time.sleep(5)
		print '----- STAR Alignment Done -----'
		st = check_strand_specific(args.inputdir+'/strandness_RSeQC.txt')
##STEP2 : RPDs

	rpds = args.rpds
	if step == 'star' or step == 'rpds':
		print '-- Strandness identification --'
		st = check_strand_specific(args.inputdir+'/strandness_RSeQC.txt')
		if st == 'unstranded':
			if not args.strandBAM:
				print 'Strand-specific BAM file is required for RPDs: Use -str/--strandBAM option'
				exit()
			else:
				print '-- Run RPDs --'
				run_rpds(args.inputdir, args.specie, args.strandBAM, jobN)
				while commands.getoutput('qstat').count('rpds') > 0: time.sleep(5)
				rpds = 'rpds'
			run_star(args.inputdir, args.specie, 'rpds', jobN)
			while commands.getoutput('qstat').count('rpds') > 0: time.sleep(5)
		else:
			rpds = 'no'
	
			while commands.getoutput('qstat').count('mapping')>0: time.sleep(5)
		print '----- Strandness identification Done -----'

	if args.strand:
		st = args.strand
	else:
		if rpds == 'rpds':
			st = check_strand_specific(args.inputdir+'/strandness_RPDs_RSeQC.txt')
		else:
			st = check_strand_specific(args.inputdir+'/strandness_RSeQC.txt')


##STEP3 : StringTie

	if step == 'star' or step == 'rpds' or step == 'stringtie':
		if st == 'unstranded':
			print 'The data is unstranded. Please run from RPDs for strand prediction.'
			exit()
		print '-- Run StringTie --'
		run_stringtie(args.inputdir, rpds, st, args.anno, jobN)
		while commands.getoutput('qstat').count('stringtie') > 0: time.sleep(5)

		print '---- StringTie Done -----'

##STEP4 : StringTie-merge

	if step == 'star' or step == 'rpds' or step == 'stringtie' or step == 'merge':
		print '-- Run StringTie-merge --'
		cluster = args.cluster
		if not cluster:
			cluster = 'single'
		run_stringtieMerge(args.inputdir, args.outputdir, jobN, cluster, args.specie)
		while commands.getoutput('qstat').count('strMerge') > 0: time.sleep(5)
		print '----- StringTie-merge Done -----'

	minputdirs = [args.outputdir+ '/assembly/' + x for x in filter(lambda x: not '.txt' in x and not '.out' in x, os.listdir(args.outputdir+'/assembly/'))]

##STEP5 : Regtools

	for i in xrange(len(minputdirs)):
		minputdir = minputdirs[i]
		sample = minputdir.split('/')[-1]
		if step == 'star' or step == 'rpds' or step == 'stringtie' or step == 'merge' or step == 'junction':
			print '-- Run Regtools --'
			if sample != 'all':
				reginputdir = ','.join([args.inputdir+ x for x in filter(lambda x: not '.txt' in x and sample in x and len(x.split(sample)[0]) == 0, os.listdir(args.inputdir))])
			else:
				reginputdir = ','.join([args.inputdir + x for x in filter(lambda x: not '.txt' in x, os.listdir(args.inputdir))])
			print reginputdir
			beddir = run_regtools(reginputdir, minputdir, st, jobN)
			print '----- Regtools Done -----'

##STEP6 : CAFE module (exj, TSS, CPS update)

	if step == 'star' or step == 'rpds' or step == 'stringtie' or step == 'merge' or step == 'junction' or step == 'cafe' or step == 'exj':
		for i in xrange(len(minputdirs)):
			minputdir = minputdirs[i]
			sample = minputdir.split('/')[-1]

			###STEP6-1 : split merged GTF file by chromosomes
			splitPath = '/home/sjpyo/new/codes/splitGTFbyChrN.py'
			gtf = minputdir + '/stringtieMerge/'+filter(lambda x: '.gtf' in x, os.listdir(minputdir+'/stringtieMerge/'))[0]
			splitGtf = 'python '+splitPath + ' ' + gtf + ' '+ args.specie + ' ' +  minputdir + '/stringtieMerge/'
			print splitGtf
			print commands.getoutput(splitGtf)
			beddir = minputdir + '/regtools/chrN/'
			chrdir = minputdir + '/stringtieMerge/'

			###STEP6-2 : Exon junction update
			print '-- Run Exon junction update --'
			cafeout = run_exj_update(chrdir, minputdir, beddir, args.specie, jobN)
		while commands.getoutput('qstat').count('exj') > 0: time.sleep(5)
	if step == 'star' or step == 'rpds' or step == 'stringtie' or step == 'merge' or step == 'junction' or step == 'cafe' or step == 'exj' or step == 'end':
		for i in xrange(len(minputdirs)):
			minputdir = minputdirs[i]
			cafeout = minputdir+'/cafe/'

			###STEP6-3 : TSS CPS update
			print '-- Run TSS CPS udate --'
			run_tss_cps_update(cafeout, args.specie, jobN)
		while commands.getoutput('qstat').count('filter') > 0: time.sleep(5)
		print '----- CAFE Done -----'
	samples = filter(lambda x: not '.out' in x and not '.txt' in x, os.listdir(args.outputdir+'/assembly/'))

##STEP7 : lncRNA annotation
	
	###STEP7-1 : known genes filtration
	if step == 'star' or step == 'rpds' or step == 'stringtie' or step == 'merge' or step == 'junction' or step == 'cafe' or step == 'exj' or step == 'end' or step == 'lncRNA':
		print '-- Run lncRNA annotation --'
		print '-- Run lncRNA filtration --'
		samples = filter(lambda x: not '.out' in x and not '.txt' in x, os.listdir(args.outputdir+'/assembly/'))
		for i in xrange(len(samples)):
			sample = samples[i]
			print sample
			linputdir = args.outputdir + '/assembly/'+sample+'/cafe/' + filter(lambda x: 'filter' in x, os.listdir(args.outputdir + '/assembly/'+sample+'/cafe/'))[0]
			loutputdir = args.outputdir + '/lncRNA_annotation/'+sample+'/'
			run_lncRNA_filtration(linputdir, loutputdir, args.specie, jobN)
		while commands.getoutput('qstat').count('lncF') >0: time.sleep(50)
	
	###STEP7-2 : CPC
	if step == 'star' or step == 'rpds' or step == 'stringtie' or step == 'merge' or step == 'junction' or step == 'cafe' or step == 'exj' or step == 'end' or step == 'lncRNA' or step == 'cpc':
		print '-- Run CPC --'
		for i in xrange(len(samples)):
			sample = samples[i]
			loutputdir = args.outputdir+'/lncRNA_annotation/'+sample+'/'
			putativeGtf = loutputdir + filter(lambda x: 'putative_transcripts' in x and '.gtf' in x, os.listdir(loutputdir))[0]
			trxN = get_TrxNum(putativeGtf)
			txnum = 300
			coutputdir = loutputdir + '/CPC/'
			run_CPC(args.specie, putativeGtf, coutputdir, txnum, jobN)
		while commands.getoutput('qstat').count('cpcsj') > int(0): time.sleep(5)
	
	###STEP7-3 : CPAT
	if step == 'star' or step == 'rpds' or step == 'stringtie' or step == 'merge' or step == 'junction' or step == 'cafe' or step == 'exj' or step == 'end' or step == 'lncRNA' or step == 'cpc' or step == 'cpat':
		print '-- Run CPAT --'
		for i in xrange(len(samples)):
			sample = samples[i]
			loutputdir = args.outputdir + '/lncRNA_annotation/'+sample+'/'
			putativeGtf = loutputdir + filter(lambda x: 'putative_transcripts' in x and '.gtf' in x, os.listdir(loutputdir))[0]
			faDir = loutputdir + '/CPC/fasta/'

			cpoutputdir = loutputdir+'/CPAT/'
			run_CPAT(faDir, cpoutputdir, jobN, args.specie)
			while commands.getoutput('qstat').count('cpat') > int(0): time.sleep(5)
	if step == 'star' or step == 'rpds' or step == 'stringtie' or step == 'merge' or step == 'junction' or step == 'cafe' or step == 'exj' or step == 'end' or step == 'lncRNA' or step == 'cpc' or step == 'cpat' or step == 'novel':
		for i in xrange(len(samples)):
			sample = samples[i]
			print sample
			loutputdir = args.outputdir + '/lncRNA_annotation/'+sample+'/'
			pTrx = loutputdir + filter(lambda x: 'putative_transcripts' in x and '.gtf' in x, os.listdir(loutputdir))[0]
			cpcdir = loutputdir + '/CPC/'
			cpatdir = loutputdir+'/CPAT/'
			outputdir = loutputdir+ '/final/'
			if not os.path.exists(outputdir): os.makedirs(outputdir)
			get_CPC_CPAT_result(cpcdir, cpatdir, pTrx, outputdir, args.specie)
			novel = outputdir + filter(lambda x: 'cpc' in x and 'cpat' in x and '.gtf' in x, os.listdir(outputdir))[0]
			outfile = outputdir+ sample +'.novel_lncRNA.gtf'
			logfile = outputdir + sample + '.novel.type.log'
			if args.specie == 'human':
				pcglncanno = '/home/sjpyo/new/data/anno/human/GENCODE/lift37/gencode.v34lift37.lncRNA+pcgs.gtf'
			elif args.specie == 'mouse':
				pcglncanno = '/home/sjpyo/new/data/anno/mouse/gencode.vM24.lncRNAs+pcgs.gtf'

			novel_filter(novel, pcglncanno, outfile, logfile)		

		print '----- lncRNA annotation Done -----'
	
	
	print '****ALL PROCESS COMPLETED****'

if __name__=="__main__":
	parser = argparse.ArgumentParser(description = '')
	parser.add_argument('-i', '--inputdir', help = 'the directory where raw data are', required = True)
	parser.add_argument('-st', '--step', help = 'Step name to start the process\nStarts from mapping (STAR) as a defualt (e.g. star, rpds, stringtie, merge, junction, cafe or lncRNA)', required = False)
	parser.add_argument('-o', '--outputdir', help = 'output directory', required = True)
	parser.add_argument('-s', '--specie', help = 'human or mouse', required = True)
	parser.add_argument('-j', '--jobN', help = 'number of jobs that will run simultaneously (multiprocessing)', required = True)
	parser.add_argument('-b', '--strandBAM', help = 'stranded BAM file for RPDs', required = False)
	parser.add_argument('-a', '--anno', help = 'reference annotation GTF file', required = True)
	parser.add_argument('-c', '--cluster', help = 'prefix to separate in merging step ("single" as defualt). If "single", no merging process.', required = False)
	parser.add_argument('-r', '--rpds', help = 'rpds information if the step starts from StringTie or after ("rpds" or "no", "no" as default)', required = False)
	parser.add_argument('-l', '--strand', help = 'strand-specific library type if the data is strand-specific (fr-firststrand or fr-secondstrand)', required = False)
#	parser.add_argument('-pr', '--jobprefix', help = 'prefix of job name', required = False)
	args = parser.parse_args()
	__main__(args)