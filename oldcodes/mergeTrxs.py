import sys, os
import module as mdl
import config as cf
import commands
inputdir = sys.argv[1]
outputdir = sys.argv[2]
if not os.path.exists(outputdir): os.makedirs(outputdir)
bedtoolsPath = '/share/apps/programs/bedtools2/2.29.0/bin/bedtools'

def get_index_positions(list_of_elems, element):
	''' Returns the indexes of all occurrences of give element in
	    the list- listOfElements '''
	index_pos_list = []
	index_pos = 0
	while True:
		try:
			# Search for item in list from indexPos to the end of list
			index_pos = list_of_elems.index(element, index_pos)
			# Add the index position in list
			index_pos_list.append(index_pos)
			index_pos += 1
		except ValueError as e: break
	return index_pos_list

def identical(trx, otherTrx):
	exons = trx.exons()
	otherexons = otherTrx.exons()
	exonsites = set(); othexonsites = set()
	for exon in exons:
		exonsites.add((exon.start(), exon.end()))
	for oexon in otherexons:
		othexonsites.add((oexon.start(), oexon.end()))
	if exonsites == othexonsites: return otherTrx
	else: return False
	

def overlapTrx(trx, otherTrx):
	if trx.chr() != otherTrx.chr(): return False
	elif not(trx.sense() == '.' or \
	otherTrx.sense() == '.' or \
	trx.sense() == otherTrx.sense()): return False
	elif trx.start() > otherTrx.end() or otherTrx.start() > trx.end(): return False
	else: return True

def overlapTrxs(trxs, chr):
	ptrxs = mdl.comp(filter(lambda x: x.sense() == '+', trxs), lambda x : (int(x.start())), False)
	mtrxs = mdl.comp(filter(lambda x: x.sense() == '-', trxs), lambda x : (int(x.start())), False)
	preov = dict()
	x = 1
## (+) strand
	for i in range(0,len(ptrxs)-1):
		trx = ptrxs[i]
		ntrx = ptrxs[i+1]
		if overlapTrx(trx,ntrx):
			if i == 0:
				newGeneid = 'BIG-imm-'+chr+'-'+str(x)
				trx.setGeneid(newGeneid)
				ntrx.setGeneid(newGeneid)
				preov[newGeneid] = [trx,ntrx]
				x += 1
			else:
				if trx.geneid() == 'BIG-imm-'+chr+'-'+str(x-1):
					ntrx.setGeneid(trx.geneid())
					preov[trx.geneid()].append(ntrx)
				else:					
					newGeneid = 'BIG-imm-'+chr+'-'+str(x)
					trx.setGeneid(newGeneid)
					ntrx.setGeneid(newGeneid)
					if not preov.has_key(newGeneid):
						preov[newGeneid] = []
					preov[newGeneid]+=[trx,ntrx]
					x += 1
		else:
			if trx.geneid() == 'BIG-imm-'+chr+'-'+str(x-1):
				continue
			else:	
				newGeneid = 'BIG-imm-'+chr+'-'+str(x)
				trx.setGeneid(newGeneid)
				preov[newGeneid] = [trx]
				x += 1
## (-) strand
	for i in range(0,len(mtrxs)-1):
		trx = mtrxs[i]
		ntrx = mtrxs[i+1]
		if overlapTrx(trx,ntrx):
			if i == 0:
				newGeneid = 'BIG-imm-'+chr+'-'+str(x)
				trx.setGeneid(newGeneid)
				ntrx.setGeneid(newGeneid)
				preov[newGeneid] = [trx,ntrx]
				x += 1
			else:
				if trx.geneid() == 'BIG-imm-'+chr+'-'+str(x-1):
					ntrx.setGeneid(trx.geneid())
					preov[trx.geneid()].append(ntrx)
				else:					
					newGeneid = 'BIG-imm-'+chr+'-'+str(x)
					trx.setGeneid(newGeneid)
					ntrx.setGeneid(newGeneid)
					if not preov.has_key(newGeneid):
						preov[newGeneid] = []
					preov[newGeneid]+=[trx,ntrx]
					x += 1
		else:
			if trx.geneid() == 'BIG-imm-'+chr+'-'+str(x-1):
				continue
			else:	
				newGeneid = 'BIG-imm-'+chr+'-'+str(x)
				trx.setGeneid(newGeneid)
				preov[newGeneid] = [trx]
				x += 1

	ovDic = dict()
	for geneid in preov.keys():
		trxs = preov[geneid]
		y = 1
		if not ovDic.has_key(geneid):
			ovDic[geneid] = []
#		trxs = mdl.comp(trxs, lambda x : (int(x.start())), False)
#		print geneid
#		print map(lambda x: x.trxid(), trxs)
		checklist = []
		iden = []
		for i in range(0,len(trxs)-1):
			if trxs[i].trxid() in checklist: continue
#			print trxs[i].trxid()
			idenCheck = [identical(trxs[i], trxs[j]) for j in range(i+1,len(trxs))]
			identrxs = filter(lambda x: x != False, idenCheck)
#			print map(lambda x: x.attri(), identrxs)
			checklist += map(lambda x: x.trxid(), identrxs)
			iden.append(identrxs)
#			print iden
		y = 0
		for i in range(len(iden)):
			idtrxs = iden[i]
			repreTrx = iden[i][0]
			y += 1
			newTrxid = geneid + '.' + str(y)
			repreTrx.setTrxid(newTrxid)
			repreTrx.setTrxtype(repreTrx.genetype())
			repreTrx.setGenetype('')
			endsup = repreTrx.attri().split('|')[0]
			types = list(set(map(lambda x: x.attri().split('|')[-1], idtrxs)))
			typeNum = len(types)

			repreTrx.setAttri(endsup+'|'+str(typeNum)+ 'types:' + ','.join(types))
			ovDic[geneid].append(repreTrx)
#		print map(lambda x: str(x), ovDic[geneid])
#	print ovDic	
	return ovDic


#		print geneid
#		print ovDic[geneid]
#			if identical(trxs[i],trxs[i+1]):
			#	print trxs[i].trxid(), trxs[i+1].trxid()
			#	print y
#				exit()
#				newTrxid = geneid + '.' + str(y)
#				sample1 = set(trxs[i].attri().split('|')[-1].split(',')); sample2 = set(trxs[i+1].attri().split('|')[-1].split(','))
#				trxs[i].setTrxid(newTrxid)
#				attri = list(sample1|sample2)
#				trxs[i].setAttri(trxs[i].attri().split('|')[0]+'|'+','.join(attri))
#				trxs[i+1] = trxs[i]
#				ovDic[geneid].append(trxs[i])
#			else:
#				ovDic[geneid].append(trxs[i])
#				y+=1
#		ovDic[geneid] = list(set(ovDic[geneid]))
#		print ovDic
#	return ovDic

def getGeneloc(trxs):
	minstart = min([int(x.start()) for x in trxs])
	maxend = max([int(x.end()) for x in trxs])
	return minstart, maxend

#def getIntersectBed(bedfile):
#	bed = open(bedfile, 'r')
#	lines = bed.readlines(); bed.close()
#	beddic = dict()
#	for i in xrange(len(lines)):
#		line = lines[i].strip().split('\t')
#		chr = line[0]; geneid = line[3]
#		tmp = line[8]
#		trxid = tmp.split('transcript_id')[0].split('"')[0]
#		celltype = tmp.split('attribute')[0].split('"')[0]
#		geneid = line[3]
#		ovtmp = line[-1]
#		ovtrxid = ovtmp.split('transcript_id')[0].split('"')[0]
#		ovtype = ovtmp.split('attribute')[0].split('"')[0]
#		if not beddic.has_key(geneid):
#			beddic[geneid] = []
#		beddic[geneid].append([trxid,ovtrxid])
#	return beddic

def mergeGtf(gtflist):
	mgtfdic = dict()
	for gtf in gtflist:
		gtfdic = mdl.getGtf(gtf)
		for chr in gtfdic.keys():
			trxs = gtfdic[chr]
			if not mgtfdic.has_key(chr):
				mgtfdic[chr] = []
			mgtfdic[chr] += trxs
	return mgtfdic

samples = filter(lambda x: not '.txt' in x and not '.log' in x and not '.out' in x, os.listdir(inputdir))
gtfall = dict()
bedlist = []
for i in xrange(len(samples)):
	sample = samples[i]
	sinputdir = inputdir +'/'+sample+'/final/'
	novelfile = sinputdir + filter(lambda x: '.gtf' in x and 'novel_lncRNA' in x and not 'gene' in x and not 'trxs' in x, os.listdir(sinputdir))[0]
####Rename every transcripts geneid and transcript id (different transcripts from different samples can use identical id)
	rename_comm = 'python /home/sjpyo/new/codes/addGenelines.py ' + novelfile + ' NOVEL ' + sample + ' ' + sample + ' ' + sinputdir + sample + '.novel_lncRNA.gene.gtf novel'
#	print rename_comm
#	print commands.getoutput(rename_comm)
	newnovelfile = sinputdir + filter(lambda x: '.gtf' in x and 'novel_lncRNA' in x and 'gene' in x, os.listdir(sinputdir))[0]
	gtf = mdl.getGtf(newnovelfile)
	for chr in gtf.keys():
		trxs = gtf[chr]
		if not gtfall.has_key(chr):
			gtfall[chr] = []
		gtfall[chr] += trxs

newgtfdic = dict()
for chr in gtfall.keys():
	print chr
	trxs = gtfall[chr]
	ovdic = overlapTrxs(trxs,chr)
	ovtrxs = sum(ovdic.values(),[])
#	newTrxs = filter(lambda x: 'BIG' in x.trxid(), ovtrxs)
	newgtfdic[chr] = ovtrxs

outfile = outputdir + '/monaco_immune.novel_lncRNA.merge.gene.gtf'
mdl.writeGtf2(newgtfdic, outfile)
#mdl.writeGtf(newgtfdic,outfile)
#	trxs_comm = '''awk -F '\t' '$3!="gene"' ''' + newnovelfile + ' > ' + sinputdir + sample + '.novel_lncRNA.trxs.gtf'
#	print trxs_comm
#	print commands.getoutput(trxs_comm)
#	novelTrx = sinputdir + sample + '.novel_lncRNA.trxs.gtf'
#
#	gtf2bed_comm = '/home/sjpyo/new/packages/bedops/bin/convert2bed -d -i gtf < ' + novelTrx + ' > ' + sinputdir + sample + '.novel_lncRNA.bed'
#	print gtf2bed_comm
#	print commands.getoutput(gtf2bed_comm)
#	bedfile = sinputdir + sample + '.novel_lncRNA.bed'
#	trxbed = sinputdir + sample + '.novel_lncRNA.trx.bed'
#	print trxbed
#	print commands.getoutput('''awk -F '\t' '$8=="transcript"' ''' + bedfile + ' > ' + trxbed)
#	bedlist.append(trxbed)
#for j in xrange(len(bedlist)-1):
#	sample = bedlist[j].split('/')[-1].split('.novel')[0]
#	bedtools_comm = bedtoolsPath + ' intersect -a ' + bedlist[j] +' -b ' + ' '.join(bedlist[j+1:]) + ' -wa -wb -s -names ' + ' '.join(map(lambda x: x.split('/')[-1].split('.novel')[0], bedlist[j+1:])) + ' > ' + outputdir + sample + '.novel.ovlp.bed'
#	print bedtools_comm
#	print commands.getoutput(bedtools_comm)
#print commands.getoutput('cat ' + outputdir + '/* > ' + outputdir + '/all.novel.ovlp.bed')
#interTrxs = dict()
#for bed in bedlist:
#	beddic = getIntersectBed(bed)
#	for geneid in beddic.keys():
#		trxs = beddic[geneid]
#		if not interTrxs.has_key(geneid):
#			interTrxs[geneid] = []
#		interTrxs[geneid] += trxs
#mgtf = mergeGtf(gtflist)
#
