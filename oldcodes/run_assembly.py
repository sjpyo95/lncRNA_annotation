import sys, os
import config as cf
import commands, time
stringtiePath = '/home/sjpyo/new/codes/run_stringtie_v2.py'
stringmergePath = '/home/sjpyo/new/codes/run_stringtieMerge.py'
splitGtfPath = '/home/sjpyo/new/codes/assembly/splitGTFbyChrN.py'
regtoolsPath = '/home/sjpyo/new/codes/run_regtools.py'
mergebedPath = '/home/sjpyo/new/codes/merge_bedfiles.py'
splitBedPath = '/home/sjpyo/new/codes/splitBEDinChrN.py'
exjPath = '/home/sjpyo/new/codes/run_exj_revised2.py'
endsupPath = '/home/sjpyo/new/codes/run_tss_cps_filter.py'

inputdir = sys.argv[1]
specie = sys.argv[2]	#hg19 or mm10
strandness = sys.argv[3]	# rf(fr-firststrand) or fr(fr-secondstrand)
anno = sys.argv[4]
outputdir = sys.argv[5]
jobN = sys.argv[6]
startpoint = sys.argv[7]
rpds = sys.argv[8]
if not os.path.exists(outputdir) : os.makedirs(outputdir)

####StringTie####
if startpoint == 'stringtie':
	print '-- Run StringTie --\n'

	stringtie = '/usr/bin/python '+stringtiePath+ ' '+ inputdir + ' ' + rpds + ' ' + strandness + ' ' + anno + ' biglab ' + jobN + ' ' + inputdir
	print stringtie
	print commands.getoutput(stringtie)
	while commands.getoutput('qstat').count('stringtie') > 0: time.sleep(5)
	print '--StringTie Finished--\n'

####StringTie-merge####
if startpoint == 'stringtie' or startpoint == 'merge':
	print '-- Run Stringtie-Merge --\n'
	prefix = inputdir.split('/')[-2]
	stringmerge = '/usr/bin/python '+stringmergePath + ' '+ inputdir + ' ' + outputdir + ' biglab '+ jobN+ ' _rep '
	print stringmerge
	print commands.getoutput(stringmerge)
	while commands.getoutput('qstat').count('strMerge') > 0: time.sleep(5)

if startpoint == 'stringtie' or startpoint == 'merge' or startpoint == 'split':
	print '(split GTF by chrN)\n'
	minputdir = [outputdir + x for x in filter(lambda x: '.out' not in x and '.txt' not in x, os.listdir(outputdir))]
	for mdir in minputdir:
		gtf = mdir + '/stringtieMerge/'+filter(lambda x: '.gtf' in x, os.listdir(mdir+'/stringtieMerge/'))[0]
		splitGtf = '/usr/bin/python '+splitGtfPath+ ' ' + gtf+ ' '+specie+ ' ' + mdir+'/stringtieMerge/'
		print splitGtf
		print commands.getoutput(splitGtf)

	print '\n-- StringTie-merge Finished --\n'

####Regtools####
if startpoint == 'regtools' or startpoint == 'stringtie' or startpoint == 'merge' or startpoint == 'split':
	print '\n-- Run Regtools --\n'
	if strandness == 'rf':  strandN = '1'
	elif strandness == 'fr': strandN = '2'
	routputdir = [outputdir + x for x in filter(lambda x: '.out' not in x and '.txt' not in x, os.listdir(outputdir))]
	for rdir in routputdir:
		name = rdir.split('/')[-1]
		rinputdir = ','.join([inputdir + x for x in filter(lambda x: name in x and len(x.split(name)[0]) == 0, os.listdir(inputdir))])
		regtools = '/usr/bin/python '+regtoolsPath + ' ' + rinputdir + ' '+ rdir+ '/regtools/ ' + strandN + ' biglab '+ jobN
		print regtools
		print commands.getoutput(regtools)
	while commands.getoutput('qstat').count('regtools') > 0: time.sleep(5)
	for rdir in routputdir:
		name = rdir.split('/')[-1]
		mergefile = rdir+'/regtools/'+name+'.merge_junction.bed'
		mergeBed = '/usr/bin/python '+ mergebedPath + ' ' + rdir+'/regtools/ '+ mergefile
		print mergeBed+'\n'
		commands.getoutput(mergeBed)
		splitBed = '/usr/bin/python '+ splitBedPath + ' ' + mergefile + ' ' + rdir+'/regtools/chrN/'
		print commands.getoutput(splitBed)
	print '-- Regtools Finished --\n'

####Exon junction update####
if startpoint == 'exj' or startpoint == 'regtools' or startpoint == 'stringtie' or startpoint == 'merge' or startpoint == 'split':
	print '\n-- Run Exon Junction Update (CAFE) --\n'
	soutputdir = [outputdir + x for x in filter(lambda x: '.out' not in x and '.txt' not in x, os.listdir(outputdir))]
	for out in soutputdir:
		minputdir = out + '/stringtieMerge/'
		beddir = out + '/regtools/chrN/'
		exj = '/usr/bin/python '+exjPath+' '+ minputdir + ' ' + beddir + ' ' + out + '/cafe/ ' + specie + ' read 2 biglab '+ jobN
		print exj
		print commands.getoutput(exj)
	while commands.getoutput('qstat').count('exj') > 0: time.sleep(5)

####End support filtering####
print '\n-- End Support Filtration --\n'
soutputdir = [outputdir + x for x in filter(lambda x: '.out' not in x and '.txt' not in x, os.listdir(outputdir))]
for out in soutputdir:
	endsup = '/usr/bin/python '+endsupPath+ ' '+ out + '/cafe/ ' + out + '/cafe/ ' + specie + ' tss near biglab ' + jobN
	print endsup + '\n'
	print commands.getoutput(endsup)
