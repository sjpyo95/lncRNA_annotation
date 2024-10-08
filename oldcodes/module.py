#!/usr/bin/python
import sys, os
import copy, string
from subprocess import Popen, PIPE, STDOUT
import cafe_dataset as data
#import re
#flag list
senseFlag = ['0', '256', '83', '163', '81', '161', '65', '129', '67', '131', '89', '137']
antisenseFlag = ['16', '272', '99', '147', '97', '145', '113', '177', '115', '179', '73', '153']

#flag table
flag_table = open('/home/bhyou/cafe/dataset/info/flag_table')
flag_dic = dict() #flag dictionary
for line in flag_table:
	if not line.startswith('flag'):
		line = line.split('\t'); flag_dic[line[0]] = []
		line2 = line[1].split('+')
		for i in range(len(line2)):
			flag_dic[line[0]].append(line2[i].rstrip())

class gene():
	def __init__(self, chr, geneid, start, end, sense, genetype, symb, attri):
		self.__chr = chr
		self.__geneid = geneid
		self.__start = start
		self.__end = end
		self.__sense = sense
		self.__genetype = genetype
		self.__symb = symb
		self.__attri = attri
	def chr(self): return self.__chr
	def geneid(self): return self.__geneid
	def start(self): return self.__start
	def end(self): return self.__end
#	def len(self): return abs(self.__end - self.__start + 1)
	def len(self): return abs(self.__end - self.__start)
	def sense(self): return self.__sense
	def genetype(self): return self.__genetype
	def symb(self): return self.__symb
	def attri(self): return self.__attri
class transcript():
	def __init__(self, chr, geneid, trxid, start, end, exons, sense, genetype, trxtype, symb, attri):
		self.__chr = chr
		self.__geneid = geneid
		self.__trxid = trxid
		self.__start = start
		self.__end = end
		self.__exons = exons
		self.__sense = sense
		self.__genetype = genetype
		self.__trxtype = trxtype
		self.__symb = symb
		self.__attri=attri
		for ex in self.__exons: ex.chgSenseTo(sense)
		self.__use = False
		self.__exj = False
		self.__fpPos = -1
		self.__tpPos = -1
		self.__overlap = 0
		self.__introns = []
		for n in range(len(exons)):
			if n < len(exons) - 1: self.__introns.append(exon(chr, exons[n].trxid(), int(exons[n].end()) + 1, int(exons[n+1].start() - 1), sense, exons[n].symb(), exons[n].attri()))
	def chr(self): return self.__chr
	def geneid(self): return self.__geneid
	def trxid(self): return self.__trxid
	def start(self): return self.__start
	def end(self): return self.__end
#	def len(self): return abs(self.__end - self.__start + 1)
	def len(self): return abs(self.__end - self.__start)
	def exons(self): return self.__exons
	def sense(self): return self.__sense
	def genetype(self): return self.__genetype
	def trxtype(self): return self.__trxtype
	def symb(self): return self.__symb
	def attri(self): return self.__attri
	def introns(self): return self.__introns
	def exonNum(self): return len(self.__exons)
	def exonLen(self): return map(lambda x: x.len(), self.__exons)
	def coord(self): return self.__chr + ':' + str(self.__start) + '-' + str(self.__end)
	def coords(self): return [self.__start, self.__end]
	def __str__(self): return self.__chr + '(' + self.__sense + '):' + '-'.join(map(str, self.coords()))
	def setGeneid(self, geneid): self.__geneid = geneid
	def setTrxid(self, trxid): self.__trxid = trxid
	def setAttri(self, attri): self.__attri = attri
	def setTrxtype(self, trxtype): self.__trxtype = trxtype
	def setGenetype(self, genetype): self.__genetype = genetype
	def setUse(self, flag): self.__use = flag
	def setExj(self, flag): self.__exj = flag
	def getUse(self): return self.__use
	def getExj(self): return self.__exj
	def setFpPos(self, pos): self.__fpPos = pos
	def setTpPos(self, pos): self.__tpPos = pos
	def getFpPos(self): return self.__fpPos
	def getTpPos(self): return self.__tpPos
	def chgSense(self, sense): self.__sense = sense
	def setExons(self, exons): self.__exons = exons
	def setOverlap(self, overlap): self.__overlap = self.__overlap + overlap
	def getOverlap(self): return self.__overlap
	def update(self, geneid, trxid, start, end, exons1, exons2):
		self.__geneid = geneid
		self.__trxid = trxid
		self.__start = start
		self.__end = end
		self.__exons = exons1 + exons2
	def overlaps(self, otherExon): #sense overlap
		if self.chr() != otherExon.chr(): return False
		elif not(self.sense() == '.' or \
			otherExon.sense() == '.' or \
			self.sense() == otherExon.sense()): return False
		elif self.start() > otherExon.end() or otherExon.start() > self.end(): return False
		else: return True
	def overlaps2(self, otherExon): #overlap
		if self.chr() != otherExon.chr(): return False
		elif self.start() > otherExon.end() or otherExon.start() > self.end(): return False
		else: return True
	def overlaps3(self, otherExon):	#antisense overlap
		if self.chr() != otherExon.chr(): return False
		elif self.sense() == otherExon.sense(): return False
		elif self.start() > otherExon.end() or otherExon.start() > self.end(): return False
		else: return True

	def toString(self):
		lines = []
		line = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; gene_type \"%s\";' % (self.chr(), 'CAFE', 'transcript', self.start(), self.end(), '.', self.sense(), '.', self.geneid(), self.trxid(), self.genetype())
		lines.append(line+'\n')
		exons = self.exons()
		for exon in exons:
			line = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; gene_type \"%s\";' % (self.chr(), 'CAFE', 'exon', exon.start(), exon.end(), '.', self.sense(), '.', self.geneid(), self.trxid(), self.genetype())
			lines.append(line+'\n')
		return lines

class exon():
	def __init__(self, chr, trxid, start, end, sense, symb, attri):
		self.__chr = chr
		self.__trxid = trxid
		self.__start = start
		self.__end = end
		self.__sense = sense
		self.__overlap = 0
		self.__symb = symb
		self.__attri=attri
	def chr(self): return self.__chr
	def trxid(self): return self.__trxid
	def start(self): return self.__start
	def end(self): return self.__end
	def sense(self): return self.__sense
	def symb(self): return self.__symb
	def attri(self): return self.__attri
	def len(self): return abs(self.__end - self.__start + 1)
	def newStart(self, pos):
		self.__start = pos
#		if self.__sense == '+': self.__start = pos
#		elif self.__sense == '-': self.__end = pos
	def newEnd(self, pos):
		self.__end = pos
#		if self.__sense == '+': self.__end = pos
#		elif self.__sense == '-': self.__start = pos
	def chgSense(self):
		if self.__sense == '+': self.__sense = '-'
		elif self.__sense == '-': self.__sense = '+'
	def chgSenseTo(self, newSense): self.__sense = newSense
	def coord(self): return self.__chr + ':' + str(self.__start) + '-' + str(self.__end)
	def coords(self): return [self.__start, self.__end]
	def __str__(self): return self.__chr + '(' + self.__sense + '):' + '-'.join(map(str, self.coords()))
	def overlaps(self, otherExon): #sense overlap
		if self.chr() != otherExon.chr(): return False
		elif not(self.sense() == '.' or \
			otherExon.sense() == '.' or \
			self.sense() == otherExon.sense()): return False
		elif self.start() > otherExon.end() or otherExon.start() > self.end(): return False
		else: return True
	def overlaps2(self, otherExon): #overlap
		if self.chr() != otherExon.chr(): return False
		elif self.start() > otherExon.end() or otherExon.start() > self.end(): return False
		else: return True
	def overlaps3(self, otherExon):	#antisense overlap
		if self.chr() != otherExon.chr(): return False
		elif self.sense() == otherExon.sense(): return False
		elif self.start() > otherExon.end() or otherExon.start() > self.end(): return False
		else: return True
	def identical(self, otherExon):
		if self.__chr == otherExon.chr() and self.__sense == otherExon.sense() and self.__start == otherExon.start() and self.__end == otherExon.end(): return True
	def inclusive(self, otherExon):
                if self.__chr == otherExon.chr() and self.__sense == otherExon.sense() and self.__start <= otherExon.start() and self.__end >= otherExon.end(): return True
	def overlapped(self, otherExon, ratio):
		if self.__chr == otherExon.chr() and self.__sense == otherExon.sense() and (min(self.__end, otherExon.end()) - max(self.__start, otherExon.start()))/float(max(self.__end, otherExon.end()) - min(self.__start, otherExon.start())) >= ratio: return True
		else: return False
	def setOverlap(self, flag): self.__overlap = flag
	def getOverlap(self): return self.__overlap

class bam():
	def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tag):
		self.__qname = qname
		self.__flag = flag
		self.__rname = rname
		self.__pos = pos
		self.__mapq = mapq
		self.__cigar = cigar
		self.__rnext = rnext
		self.__pnext = pnext
		self.__tlen = tlen
		self.__seq = seq
		self.__qual = qual
		self.__tag = tag
	def qname(self): return self.__qname
	def flag(self): return self.__flag
	def rname(self): return self.__rname
	def pos(self): return self.__pos
	def mapq(self): return self.__mapq
	def cigar(self): return self.__cigar
	def rnext(self): return self.__rnext
	def pnext(self): return self.__pnext
	def tlen(self): return self.__tlen
	def seq(self): return self.__seq
	def qual(self): return self.__qual
	def tag(self): return self.__tag
	def setTag(self, tag): self.__tag = tag
	def setreadRatio(self, fReadN, rReadN):
		self.__fReadN = fReadN
		self.__rReadN = rReadN
		self.__tReadN = self.__fReadN + self.__rReadN
		if self.__tReadN == 0: self.__readRatio = 0.5
		else: self.__readRatio = round((self.__fReadN / self.__tReadN), 3)
	def readRatio(self): return self.__readRatio
	def setreadProb(self, fProb, rProb):
		self.__fProb = float(fProb)
		self.__rProb = float(rProb)
		if self.__fProb == 0 and self.__rProb == 0: self.__readProb = 1
		elif self.__fProb == 0: self.__readProb = -1111 ##check
		elif self.__rProb == 0: self.__readProb = 1111 ##check
		else: self.__readProb = round((self.__fProb / self.__rProb), 3)
	def readProb(self): return self.__readProb
	def readfProb(self): return self.__fProb
	def readrProb(self): return self.__rProb

def cmp0(a1, a2): return cmp(a1.start(), a2.start())

def getMatrix(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	header = lines[0].strip()
	samples = header.split('\t')[2:]

#	print samples
	geneDic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		id = line[0]; symb = line[1]; exp = line[2:]
#		print exp
		if not geneDic.has_key((id,symb)):
			geneDic[(id, symb)] = dict()
		for x in xrange(len(exp)):
			sample = samples[x]; express = exp[x]
			geneDic[(id,symb)][sample] = express
	return geneDic
def comp(sampleList, stand, torf): return sorted(sampleList, key=stand, reverse = torf)

def writeMatrix(matdic, filename):
	outfile = open(filename, 'w')
	samples = matdic.values()[0].keys()
	outfile.write('ID\tsymbol\t'+ '\t'.join(samples)+'\n')
	for idsym in matdic.keys():
		express = matdic[idsym]
		outfile.write('\t'.join(idsym))
		for sample in samples:
			outfile.write('\t'+express[sample])
		outfile.write('\n')	
	outfile.close()

def getGene(filename):
	filein = open(filename)
	lines = filein.readlines(); filein.close()
	gtfChr = dict()
	preChr = ''; preGeneID = ''; preTrxID = ''; preStart = -1; preEnd = -1; preSense = ''; preGenetype = ''; preSymb = ''; preAttri = ''
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('#'):
			line = line.strip().split('\t')
			if not line[2] == 'gene': continue
			chrom = line[0]; temp = line[-1].strip()
			geneid = temp.split('gene_id')[1].split(';')[0].split('"')[1]
			preChr = chrom; preStart = int(line[3].split(' ')[0]); preEnd = int(line[4].split(' ')[0]); preSense = line[6].split(' ')[0]
			if temp.find('gene_type') > -1:
				genetype = temp.split('gene_type')[1].split(';')[0].split('"')[1]
			else:
				genetype = ""
			if temp.find('gene_name') > -1:
				symb = temp.split('gene_name')[1].split(';')[0].split('"')[1]
			else:
				symb = ""
			if temp.find('attribute') > -1:
				attri = temp.split('attribute')[1].split(';')[0].split('"')[1]
			else:
				attri = ''
			if not(gtfChr.has_key(chrom)): gtfChr[chrom] = []
			preGeneID = geneid; preGenetype = genetype; preSymb = symb; preAttri = attri
			ngene = gene(preChr, preGeneID, preStart, preEnd, preSense, preGenetype, preSymb, preAttri)
			if preChr != '': gtfChr[preChr].append(ngene)
	return gtfChr

def getGtf(filename):
	filein = open(filename)
	lines = filein.readlines()
	#lines = lines.split('\n'); filein.close()
	gtfChr = dict(); exons = []
	preChr = ''; preGeneID = ''; preTrxID = ''; preStart = -1; preEnd = -1; preSense = ''; preGenetype = ''; preTrxtype = ''; preSymb = ''; preAttri = ''
#	for i in xrange(len(lines[:-1])):
	for i in xrange(len(lines)):
		line = lines[i].strip()
		if len(line) != 0 and not line.startswith('#'):
			line = line.strip().split('\t')
			if line[2] == 'gene': continue
			#print line
#			chrom = line[0].split(' ')[0]; temp = line[8].split('"')
			chrom = line[0]; temp = line[-1].strip()
			
#			geneid = temp[1]; trxid = temp[3]
#			geneid = re.search('.*gene_id "(\S+)";', temp).group(1)
#			trxid = re.search('.*transcript_id "(\S+)";', temp).group(1)
			#print temp, geneid, trxid
#			geneid = line[-1].split(';')[0].split('"')[1]
#			trxid = line[-1].split(';')[1].split('"')[1]
			geneid = temp.split('gene_id')[1].split(';')[0].split('"')[1]
			trxid = temp.split('transcript_id')[1].split(';')[0].split('"')[1]

			if temp.find('gene_type') > -1:
				genetype = temp.split('gene_type')[1].split(';')[0].split('"')[1]
			else:
				genetype = ""
			if temp.find('transcript_type') > -1:
				trxtype = temp.split('transcript_type')[1].split(';')[0].split('"')[1]
			else:
				trxtype = ""
			if temp.find('gene_name') > -1:
				symb = temp.split('gene_name')[1].split(';')[0].split('"')[1]
			else:
				symb = ""
			if temp.find('attribute') > -1:
				attri = temp.split('attribute')[1].split(';')[0].split('"')[1]
			else:
				attri = ''
			if not(gtfChr.has_key(chrom)): gtfChr[chrom] = []
			if i > 0 and line[2] == 'transcript':
				#print line[2]
				exons.sort(cmp0)
				ntranscript = transcript(preChr, preGeneID, preTrxID, preStart, preEnd, exons, preSense, preGenetype, preTrxtype, preSymb, preAttri)
				if preChr != '': gtfChr[preChr].append(ntranscript)
				exons = []
			elif line[2] == 'exon':
				#print line[2]
				nexon = exon(chrom, trxid, int(line[3]), int(line[4]), line[6], symb, attri); exons.append(nexon)
			preGeneID = geneid; preTrxID = trxid; preGenetype = genetype; preTrxtype = trxtype; preSymb = symb; preAttri = attri
			if line[2] == 'transcript': preChr = chrom; preStart = int(line[3].split(' ')[0]); preEnd = int(line[4].split(' ')[0]); preSense = line[6].split(' ')[0]
			if len(lines) - 2 == i:
				exons.sort(cmp0)
				ntranscript = transcript(preChr, preGeneID, preTrxID, preStart, preEnd, exons, preSense, preGenetype, preTrxtype, preSymb, preAttri)
				if preChr != '': gtfChr[preChr].append(ntranscript)
	return gtfChr

#def getGtf2(filename):
#        filein = open(filename)
#        lines = filein.read(); lines = lines.split('\n'); filein.close()
#        gtfChr = dict(); exons = []
#        preChr = ''; preGeneID = ''; preTrxID = ''; preStart = -1; preEnd = -1; preSense = ''
##       for i in xrange(len(lines[:-1])):
#        for i in xrange(len(lines[:])):
#                line = lines[i]
#                if len(line) != 0 and not line.startswith('#'):
#                        line = line.split('\t')
##                       chrom = line[0]; temp = line[-1].split('\"')
#						chrom = line[0]; temp = line[-1].strip()
##                       geneid = temp[1]; trxid = temp[3]
#						geneid = temp.strip().split('gene_id')[1].split(';')[0].split('"')[1]
#						trxid = temp.strip().split('transcript_id')[1].split(';')[0].split('"')[1]
##						geneid = re.search('.*gene_id "(\S+)";', temp).group(1)
##						trxid = re.search('.*transcript_id "(\S+)";', temp).group(1)
#						if not(gtfChr.has_key(chrom)): gtfChr[chrom] = []
#						if i > 0 and line[2] == 'transcript':
#								exons.sort(cmp0)
#                                ntranscript = transcript(preChr, preGeneID, preTrxID, preStart, preEnd, exons, preSense)
#                                if preChr != '': gtfChr[preChr].append(ntranscript)
#                                exons = []
#                        elif line[2] == 'exon':
#                                nexon = exon(chrom, int(line[3]), int(line[4]), line[6]); exons.append(nexon)
#                        preGeneID = geneid; preTrxID = trxid
#                        if line[2] == 'transcript': preChr = chrom; preStart = int(line[3]); preEnd = int(line[4]); preSense = line[6]
#                        if len(lines) - 2 == i:
#                                exons.sort(cmp0)
#                                ntranscript = transcript(preChr, preGeneID, preTrxID, preStart, preEnd, exons, preSense)
#                                if preChr != '': gtfChr[preChr].append(ntranscript)
#        return gtfChr

def writeGtf(gtfD, outfilename):
	fileout = open(outfilename, 'w')
#	fileout = open(outfilename, 'a')
	lines = []
	for chr in gtfD.keys():
		gtf = gtfD[chr]
		for trx in gtf:
#			line = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\";' % (chr, 'CorrectedCuff', 'transcript', trx.start(), trx.end(), '.', trx.sense(), '.', trx.geneid(), trx.trxid())	
			line = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; gene_type \"%s\"; transcript_type \"%s\"; attribute \"%s\";' % (chr, 'CAFE', 'transcript', trx.start(), trx.end(), '.', trx.sense(), '.', trx.geneid(), trx.trxid(), trx.symb(), trx.genetype(), trx.trxtype(), trx.attri())	
			lines.append(line)
			exons = trx.exons()
			for ex in exons:
#				line = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\";' % (chr, 'CorrectedCuff', 'exon', ex.start(), ex.end(), '.', trx.sense(), '.', trx.geneid(), trx.trxid())
				line = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; gene_type \"%s\"; transcript_type \"%s\"; attribute \"%s\";' % (chr, 'CAFE', 'exon', ex.start(), ex.end(), '.', trx.sense(), '.', trx.geneid(), trx.trxid(), trx.symb(), trx.genetype(), trx.trxtype(), trx.attri())
				lines.append(line)
	fileout.write('\n'.join(lines))
	fileout.write('\n')
	fileout.close()

def getGeneloc(trxs):
	minstart = min([int(x.start()) for x in trxs])
	maxend = max([int(x.end()) for x in trxs])
	return minstart, maxend

def writeGtf2(gtfD, outfilename):
	fileout = open(outfilename, 'w')
	lines = []
	geneDic = dict()
	for chr in gtfD.keys():
		gtf = gtfD[chr]
		for trx in gtf:
			geneid = trx.geneid()
			if not geneDic.has_key(geneid):
				geneDic[geneid] = []
			geneDic[geneid].append(trx)
	for geneid in geneDic.keys():
		trxs = geneDic[geneid]
		attributes = map(lambda x: x.attri().split('|'), trxs)
		endsup = attributes[0][0]; typeNum = len(list(set(sum(map(lambda x: x[-1].split('types')[-1].split(','), attributes), []))))
		geneattri = endsup+'|'+str(typeNum)+ 'types'
		geneStart, geneEnd = getGeneloc(trxs)
		genetype = ''
		geneline = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; gene_type \"%s\"; gene_name \"%s\"; attribute \"%s\";' % (trxs[0].chr(), 'BIG', 'gene', geneStart, geneEnd, '.', trxs[0].sense(), '.', trxs[0].geneid(), genetype, trxs[0].symb(), geneattri)
		lines.append(geneline)
		for trx in trxs:
			line = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; gene_type \"%s\"; transcript_type \"%s\"; attribute \"%s\";' % (trx.chr(), 'BIG', 'transcript', trx.start(), trx.end(), '.', trx.sense(), '.', trx.geneid(), trx.trxid(), trx.symb(), trx.genetype(), trx.trxtype(), trx.attri())	
			lines.append(line)
			exons = trx.exons()
			for ex in exons:
#				line = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\";' % (chr, 'CorrectedCuff', 'exon', ex.start(), ex.end(), '.', trx.sense(), '.', trx.geneid(), trx.trxid())
				line = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; gene_type \"%s\"; transcript_type \"%s\"; attribute \"%s\";' % (trx.chr(), 'BIG', 'exon', ex.start(), ex.end(), '.', trx.sense(), '.', trx.geneid(), trx.trxid(), trx.symb(), trx.genetype(), trx.trxtype(), trx.attri())
				lines.append(line)
	fileout.write('\n'.join(lines))
	fileout.write('\n')
	fileout.close()
		

def bothStrand(trx):
	ftrx = transcript(trx.chr(), trx.geneid(), trx.trxid() + '+', trx.start(), trx.end(), trx.exons(), '+')
	rtrx = transcript(trx.chr(), trx.geneid(), trx.trxid() + '-', trx.start(), trx.end(), trx.exons(), '-')
	return ftrx, rtrx	

def flagCorrection(read):
	senseflag_list = ['99', '147', '97', '145', '67', '131', '65', '129']
	antisenseflag_list = ['83', '163', '81', '161', '115', '179', '113', '177']
	nflag_dic = {'83':'99', '163':'147', '81':'97', '161':'145', '115':'67', '179':'131', '113':'65', '177':'129'}

	if int(read.flag()) >= 256: #the alingment is not primary
		nflag = int(read.flag()) - 256
	else: 
		nflag = int(read.flag())

	if float(read.readRatio()) >= 0.99: #forward strand read
		if str(nflag) == '0': nflag = '73'
		elif '8' in flag_dic[str(nflag)]: #the mate is unmapped
			if (str(nflag) == '73') or (str(nflag) == '153'): pass
			else:
				if str(nflag) == '89': nflag = '73'
#				elif str(nflag) == '137': nflag = '153'
				else: nflag = '153'
		else: #the mate is mapped
			if str(nflag) in senseflag_list: pass
			else:
				if str(nflag) in antisenseflag_list: nflag = nflag_dic[str(nflag)]
	elif float(read.readRatio()) <= 0.01: #reverse strand read
		if str(nflag) == '16': nflag = '89'
		elif '8' in flag_dic[str(nflag)]: #the mate is unmapped
			if (str(nflag) == '89') or (str(nflag) == '137'): pass
			else:
				if str(nflag) == '73': nflag = '89'
#				elif str(nflag) == '153': nflag = '137'
				else: nflag = '137'
		else: #the mate is mapped
			if str(nflag) in antisenseflag_list: pass
			else:
				if str(nflag) in senseflag_list:
					for flag_key in nflag_dic.keys():
						if str(nflag) == nflag_dic[flag_key]: nflag = flag_key
	else: pass

	ntlen = int(read.tlen()) #read length
	if '16' in flag_dic[str(nflag)]: ntlen = -abs(ntlen)
	else: ntlen = abs(ntlen)
	
	if int(read.flag()) >= 256: nflag = int(nflag) + 256
	nread = bam(read.qname(), int(nflag), read.rname(), read.pos(), read.mapq(), read.cigar(), read.rnext(), read.pnext(), ntlen, read.seq(), read.qual(), read.tag())
	return nread

def checkDirection(bamfile, locus, libType):
	senseReads = []; antisenseReads = []
	commandBase = 'samtools view ' + bamfile + ' ' + locus.coord()
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	output = p.stdout.read(); output = output.split('\n')
	#excluding junction hits when determing a direction
	mappingReads = filter(lambda x: x.split('\t')[5].find('N') < 0, output[:-1])
	if len(mappingReads) > 0:
		for mappingRead in mappingReads:
			mappingRead = mappingRead.split('\t')
			if int(mappingRead[1]) >= 256: flag = str(int(mappingRead[1]) - 256) 
			else: flag = str(mappingRead[1])
			if libType == 'fr': #fr-firststrand
				if flag in senseFlag: senseReads.append(mappingRead)
				elif flag in antisenseFlag: antisenseReads.append(mappingRead)
			else: #fr-secondstrand
				if flag in antisenseFlag: senseReads.append(mappingRead)
				elif flag in senseFlag: antisenseReads.append(mappingRead)
	return senseReads, antisenseReads

def checkDirectionTx(sp_bamfile, locus, libType = 'fr'):
	senseReads = []; antisenseReads = []
	exons = locus.exons()
	for exon in exons:
		esenseReads, eantisenseReads = checkDirection(sp_bamfile, exon, libType)
		senseReads += esenseReads; antisenseReads += eantisenseReads
	totalReadN = len(senseReads) + len(antisenseReads)
	if totalReadN > 0: strandRatio = round((len(senseReads)/float(totalReadN)), 3)
	else: strandRatio = 'null'
	return strandRatio, senseReads, antisenseReads

def checkDirectionEx(sp_bamfile, locus, libType = 'fr'):
	senseReads = []; antisenseReads = []
	senseReads, antisenseReads = checkDirection(sp_bamfile, locus, libType)
	totalReadN = len(senseReads) + len(antisenseReads)
	if totalReadN > 0: strandRatio = round((len(senseReads)/float(totalReadN)), 3)
	else: strandRatio = 'null'
	return strandRatio, senseReads, antisenseReads

def getPos(read):
	pos = int(read[3]); cigar = read[5]; last = 0
	N = [0,0]; I = [0,0]; D = [0,0]
	tags = ['N', 'I', 'D', 'M']

	if int(read[1]) >= 256:
		flag = int(read[1]) - 256 #the alignment is not primary
	else:
		flag = int(read[1])

	if '16' in flag_dic[str(flag)]:	#reverse strand
		length = len(read[9])
		for z in range(len(cigar)):
			if cigar[z] in tags:
				if cigar[z] == 'N':
					n = int(cigar[last+1:z])
					N.append(n)
				elif cigar[z] == 'D':
					d = int(cigar[last+1:z])
					D.append(d)
				elif cigar[z] == 'I':
					i = int(cigar[last+1:z])
					I.append(i)
				last = z
			else: pass
		pos = pos + length - 1 + sum(N) + sum(D) - sum(I)
	return pos

def getPos_sp(reads, libType = 'fr'): #stranded RNA-seq reads
	readPos = dict()
	for read in reads:
		read = read.split('\t'); pos = getPos(read); sense = ''
		if int(read[1]) >= 256: flag = str(int(read[1]) - 256)
		else: flag = str(read[1])
		
		if libType == 'fr': #fr-firststrand
			if flag in senseFlag: sense = '+'
			elif flag in antisenseFlag: sense = '-'
		else: #fr-secondstrand
			if flag in antisenseFlag: sense = '+'
			elif flag in senseFlag: sense = '-'

		if not readPos.has_key(pos):
			readPos[pos] = dict() #readPos => {pos:dict()}
		if not readPos[pos].has_key(sense):
			readPos[pos][sense] = []
			readPos[pos][sense].append(read)
		else:
			readPos[pos][sense].append(read)
	return readPos

def getPos_np(reads): #unstranded RNA-seq reads
	readPos = dict()
	for read in reads:
		read = read.split('\t'); pos = getPos(read); pair = ''
		if int(read[1]) >= 256: flag = str(int(read[1]) - 256)
		else: flag = str(read[1])

		if '64' in flag_dic[str(flag)]: pair = 'first'
#		elif '128' in flag_dic[str(flag)]: pair = 'second'
		else: pair = 'second'

		if not readPos.has_key(pos):
			readPos[pos] = dict() #readPos => {pos:dict()}
		if not readPos[pos].has_key(pair):
			readPos[pos][pair] = []
			readPos[pos][pair].append(read)
		else:
			readPos[pos][pair].append(read)
	return readPos

#def getAdjacentPos(npPos, spReads_list): #distance == 0
#	pos_dict = dict()
#	spReads_list = filter(lambda x: int(x) == int(npPos), spReads_list) #distance != 0
#	for x in range(len(spReads_list)):
#		distance = abs(int(npPos) - int(spReads_list[x]))
#		if not pos_dict.has_key(distance): pos_dict[distance] = []
#		pos_dict[distance].append(spReads_list[x])
#	return pos_dict

def getAdjacentPos(npPos, spReads_list): #distance == 0
	nspReads_list = []; pos_dict = dict()
#	spReads_list = filter(lambda x: int(x) == int(npPos), spReads_list) #distance != 0
	distance = abs(int(npPos) - int(spReads_list[0]))
	for x in range(1, len(spReads_list)):
		pre_distance = abs(int(npPos) - int(spReads_list[x]))
		if pre_distance < distance: distance = pre_distance
		else: break
	for y in range(max(0, x-10), min(x+10, len(spReads_list))):
		ndistance = abs(int(npPos) - int(spReads_list[y]))
		if not pos_dict.has_key(ndistance): pos_dict[ndistance] = []
		pos_dict[ndistance].append(spReads_list[y])
	return pos_dict

def getAdjacentPos2(npPos, spReads_list):
	pos_dict = getAdjacentPos(npPos, spReads_list)
	pos_list = pos_dict.keys(); pos_list.sort()
	return pos_dict[pos_list[0]], pos_list[0]

def getAdjacentPos3(npPos, snReads, spReads_list):
	nspReads_list = []
	for spRead in spReads_list:
		if not spRead in snReads: nspReads_list.append(spRead)
	pos_dict = getAdjacentPos(npPos, nspReads_list)
	pos_list = pos_dict.keys(); pos_list.sort()
	return pos_list[0], pos_dict[pos_list[0]]

def getAdjacentPos4(snReads_dict, spReads_list, x):
	nsnReads_dict = dict(); nsnReads_dict[x] = []
	flag = -1
	for snReads in snReads_dict[x-1]:
		distance, apReads_list = getAdjacentPos3(snReads[x-1], snReads, spReads_list)
		for snRead in apReads_list:
			nsnReads = copy.copy(snReads); nsnReads.append(snRead); nsnReads_dict[x].append(nsnReads)
		if distance > 50: flag = 1
	if flag < 0: return nsnReads_dict
	else: return snReads_dict

def return_sense(sense_list, sense_list2, x):
	nsense_list = []
	if x == 0:
		for sense in sense_list2[0]:
			nsense_list.append([sense])
	else:
		for sense1 in sense_list:
			for sense2 in sense_list2[x]:
				nsense_list.append(sense1 + [sense2])
	return nsense_list

def makeSam(read):
	read = read.split('\t')
	qname = read[0]; flag = read[1]; rname = read[2]; pos = read[3]; mapq = read[4]; cigar = read[5]
	rnext = read[6]; pnext = read[7]; tlen = read[8]; seq = read[9]; qual = read[10]; tag = read[11:]
	nread = bam(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tag)
	return nread

def writeSamHeader(outputFile, database): #header for the SAM
#	outputFile = open(outputFile, 'a')
	if database == 'hg19':
                outputFile.write('@HD\tVN:1.0\tSO:coordinate\n')
                outputFile.write('@SQ\tSN:chr1\tLN:249250621\n')
                outputFile.write('@SQ\tSN:chr10\tLN:135534747\n')
                outputFile.write('@SQ\tSN:chr11\tLN:135006516\n')
                outputFile.write('@SQ\tSN:chr12\tLN:133851895\n')
                outputFile.write('@SQ\tSN:chr13\tLN:115169878\n')
                outputFile.write('@SQ\tSN:chr14\tLN:107349540\n')
                outputFile.write('@SQ\tSN:chr15\tLN:102531392\n')
                outputFile.write('@SQ\tSN:chr16\tLN:90354753\n')
                outputFile.write('@SQ\tSN:chr17\tLN:81195210\n')
                outputFile.write('@SQ\tSN:chr18\tLN:78077248\n')
                outputFile.write('@SQ\tSN:chr19\tLN:59128983\n')
                outputFile.write('@SQ\tSN:chr2\tLN:243199373\n')
                outputFile.write('@SQ\tSN:chr20\tLN:63025520\n')
                outputFile.write('@SQ\tSN:chr21\tLN:48129895\n')
                outputFile.write('@SQ\tSN:chr22\tLN:51304566\n')
                outputFile.write('@SQ\tSN:chr3\tLN:198022430\n')
                outputFile.write('@SQ\tSN:chr4\tLN:191154276\n')
                outputFile.write('@SQ\tSN:chr5\tLN:180915260\n')
                outputFile.write('@SQ\tSN:chr6\tLN:171115067\n')
                outputFile.write('@SQ\tSN:chr7\tLN:159138663\n')
                outputFile.write('@SQ\tSN:chr8\tLN:146364022\n')
                outputFile.write('@SQ\tSN:chr9\tLN:141213431\n')
                outputFile.write('@SQ\tSN:chrM\tLN:16571\n')
                outputFile.write('@SQ\tSN:chrX\tLN:155270560\n')
                outputFile.write('@SQ\tSN:chrY\tLN:59373566\n')
#               outputFile.write('@PG\tID:TopHat\tVN:2.0.6\tCL:strandEstimation\n')
	elif database == 'mm9':
                outputFile.write('@HD\tVN:1.0\tSO:coordinate\n')
                outputFile.write('@SQ\tSN:chr1\tLN:197195432\n')
                outputFile.write('@SQ\tSN:chr10\tLN:129993255\n')
                outputFile.write('@SQ\tSN:chr11\tLN:121843856\n')
                outputFile.write('@SQ\tSN:chr12\tLN:121257530\n')
                outputFile.write('@SQ\tSN:chr13\tLN:120284312\n')
                outputFile.write('@SQ\tSN:chr14\tLN:125194864\n')
                outputFile.write('@SQ\tSN:chr15\tLN:103494974\n')
                outputFile.write('@SQ\tSN:chr16\tLN:98319150\n')
                outputFile.write('@SQ\tSN:chr17\tLN:95272651\n')
                outputFile.write('@SQ\tSN:chr18\tLN:90772031\n')
                outputFile.write('@SQ\tSN:chr19\tLN:61342430\n')
                outputFile.write('@SQ\tSN:chr2\tLN:181748087\n')
                outputFile.write('@SQ\tSN:chr3\tLN:159599783\n')
                outputFile.write('@SQ\tSN:chr4\tLN:155630120\n')
                outputFile.write('@SQ\tSN:chr5\tLN:152537259\n')
                outputFile.write('@SQ\tSN:chr6\tLN:149517037\n')
                outputFile.write('@SQ\tSN:chr7\tLN:152524553\n')
                outputFile.write('@SQ\tSN:chr8\tLN:131738871\n')
                outputFile.write('@SQ\tSN:chr9\tLN:124076172\n')
                outputFile.write('@SQ\tSN:chrM\tLN:16299\n')
                outputFile.write('@SQ\tSN:chrX\tLN:166650296\n')
                outputFile.write('@SQ\tSN:chrY\tLN:15902555\n')
                outputFile.write('@PG\tID:TopHat\tVN:2.0.6\tCL:strandEstimation\n')
#	outputFile.close()

def writeSam(read):
        readLine = read.qname() + '\t' + `int(read.flag())` + '\t' + read.rname() + '\t' + `int(read.pos())` + '\t' + `int(read.mapq())` + '\t' + read.cigar() + '\t' + read.rnext() + '\t' + `int(read.pnext())` + '\t' + `int(read.tlen())` + '\t' + read.seq() + '\t' + read.qual() + '\t' + '\t'.join(read.tag())
        return readLine

def samToBam(samfilename):
        bamfilename = '.'.join(samfilename.split('.')[:-1]) + '.bam'
        commandBase = 'samtools view -Sb ' + samfilename + ' > ' + bamfilename
        print commandBase
        p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        p.wait()
        os.remove(samfilename)

def samtoolsCat(bamfilelist, outputfilename):
        commandBase = 'samtools cat -o ' + outputfilename + ' ' + ' '.join(bamfilelist)
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()
	samtoolsSort(outputfilename)
	samtoolsIndex(outputfilename)

def samtoolsSort(outputfilename):
	sortedFile = '.'.join(outputfilename.split('.')[:-1]) + 'Sorted'
#	commandBase = 'samtools sort ' + outputfilename + ' ' + sortedFile
	commandBase = 'samtools sort -@ 4 ' + outputfilename + ' ' + sortedFile
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()
	os.remove(outputfilename)
	os.renames(sortedFile+'.bam', outputfilename)

def samtoolsIndex(bamfilename):
	commandBase = 'samtools index ' + bamfilename
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()

def samtoolsRetrieve(bamfilename, coord):
        commandBase = 'samtools view ' + bamfilename + ' ' + coord
        p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        output = p.stdout.read()
        return output

def getBam(bamfile, coord, libType = 'np'): #distinguish between mappingReads and exonJunctionReads
	commandBase = 'samtools view ' + bamfile + ' ' + coord
	print commandBase
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	output = p.stdout.read(); output = output.split('\n')
	mappingReads = filter(lambda x: x.split('\t')[5].find('N') < 0, output[:-1])
	exonJunctionReads = filter(lambda x: x.split('\t')[5].find('N') > 0, output[:-1])
	return mappingReads, exonJunctionReads

def getBamDict(bamfile, coord, libType = 'np'): #distinguish between mappingReads and exonJunctionReads
	mappingReads_dic = dict(); exonJunctionReads_dic = dict(); readPos_dic = dict()
	commandBase = 'samtools view ' + bamfile + ' ' + coord
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	output = p.stdout.read(); output = output.split('\n'); output = output[:-1]
	for read in output:
		nread = makeSam(read)
		if read.split('\t')[5].find('N') < 0: #mappingReads
			if not mappingReads_dic.has_key(nread.qname()):
				mappingReads_dic[nread.qname()] = dict()
			mappingReads_dic[nread.qname()][read] = ''

			read = read.split('\t'); pos = getPos(read); info = ''
			if int(read[1]) >= 256: flag = str(int(read[1]) - 256)
			else: flag = str(read[1])

			if libType == 'np': #unstranded RNA-seq reads
				if '64' in flag_dic[str(flag)]: info = 'first'
				else: info = 'second'
			else: #stranded RNA-seq reads
				if libType == 'fr': #fr-firststrand
					if flag in senseFlag: info = '+'
					elif flag in antisenseFlag: info = '-'
				else: #fr-secondstrand
					if flag in antisenseFlag: info = '+'
					elif flag in senseFlag: info = '-'

			if not readPos_dic.has_key(pos):
				readPos_dic[pos] = dict() #readPos_dic => {pos:dict()}
			if not readPos_dic[pos].has_key(info):
				readPos_dic[pos][info] = []
			readPos_dic[pos][info].append(read)
		else: #exonJunctionReads
			if 'XS:A:+' in read.split('\t')[11:]: nread.setreadRatio(1,0) #forward strand
			else: nread.setreadRatio(0,1) #reverse strand
			nread = flagCorrection(nread)
			if not exonJunctionReads_dic.has_key(nread.qname()):
				exonJunctionReads_dic[nread.qname()] = dict()
			exonJunctionReads_dic[nread.qname()][read] = nread
	return mappingReads_dic, exonJunctionReads_dic, readPos_dic

def realJunctions(junctionReads, coord):
	coord = coord.split(':')
	start = coord[1].split('-')[0]; end = coord[1].split('-')[1]
	tags = ['N', 'I', 'D', 'M']; njunctionsReads = []
	for junctionRead in junctionReads:
		N = [0,0]; I = [0,0]; D = [0,0]
		junctionread = junctionRead.split('\t'); readStart = int(junctionread[3])
		readLength = len(junctionread[9]); cigar = junctionread[5]

		last = 0
		for i in range(len(cigar)):
			if cigar[i] in tags:
				if cigar[i]=='N':
					n = int(cigar[last+1:i])
					N.append(n)
				elif cigar[i]=='D':
					d = int(cigar[last+1:i])
					D.append(d)
				elif cigar[i]=='I':
					ii = int(cigar[last+1:i])
					I.append(ii)
				last = i
			else: pass
		readEnd = int(readStart) + sum(N) + sum(D) - sum(I) + int(readLength) - 1
		if (int(readStart) < int(start)) and (int(readEnd) > int(end)): pass
		else: njunctionsReads.append(junctionRead)
	return njunctionsReads

def cmp1(a1, a2): return cmp(a1[0], a2[0])
#def updateJunctions(bamfile, exon1, exon2, libType = 'sp'):
def updateJunctions(exon1, exon2, junctions, sense, libType = 'sp'):
	if exon1.end() >= exon2.start():
		print exon1, exon2
		print 'updateJunctions Error: start is bigger than end'
		return False, [exon1, exon2]
#	coord = exon1.chr() + ':' + str(exon1.end()-1) + '-' + str(exon1.end()+1)
#	output = samtoolsRetrieve(bamfile, coord); output = output.split('\n')
	coord = exon(exon1.chr(), exon1.start()+1, exon2.end()-1, sense)

	def keySort(a1, a2): return cmp(int(a1.split('.')[0]), int(a2.split('.')[0]))
	def getIntrons(lines):
		intronHits = dict()
		for line in lines:
			line = line.split('\t')
			if len(line[5].split('M')[0].split('S')) == 2:
				fMatchN = int(line[5].split('M')[0].split('S')[1])
			else:
				fMatchN = int(line[5].split('M')[0])
			interval = int(line[5].split('M')[1].split('N')[0])
			fiveIntron  = int(line[3]) + fMatchN
			threeIntron = int(line[3]) + fMatchN + interval - 1
			index = str(fiveIntron) + '.' + str(threeIntron)
#			if line[17][-1] == exon1.sense(): #/// for only a same direction
			bamFlag = int(line[1])
			if bamFlag >= 256:
				bamFlag = bamFlag - 256
			flagVals =[128,64,32,16,8,4,2]
			checkVals =[]
			#print bamFlag
			for i in flagVals:
				if bamFlag >= i:
					#print bamFlag, i
					bamFlag = bamFlag - i
					checkVals.append(i)
				else:
					continue
			#print checkVals
			if 16 in checkVals:
				bamSense = '-'
			else:
				bamSense = '+'
			#print bamSense, exon1.sense()
			if bamSense == exon1.sense() : #/// for only a same direction
				if not(intronHits.has_key(index)): intronHits[index] = 0
				intronHits[index] += 1
		#/// choose only introns with more than 2 events
		intronHitA = filter(lambda x: intronHits[x] > 1, intronHits.keys())
		## checking overlapping introns. If it is, choose intron with most hits
		intronHitA2 = []
		for i in xrange(len(intronHitA)):
			intron = intronHitA[i]
			fIntron, tIntron = map(int, intron.split('.'))
			flag = True
			for j in xrange(len(intronHitA)):
				intron2 = intronHitA[j]
				fIntron2, tIntron2 = map(int, intron2.split('.'))
				if i != j:
					if min(tIntron, tIntron2) - max(fIntron, fIntron2) >= 0:
						if intronHits[intron] < intronHits[intron2]: flag = False
			if flag: intronHitA2.append(intron)
		intronHitA2 = map(lambda x: map(int, x.split('.')), intronHitA2)
		intronHitA2.sort(cmp1)
		#print intronHitA2
		return intronHitA2

	## extracting junction hits
	junctionHits = filter(lambda x: x.split('\t')[5].find('N') > 0 and x.split('\t')[5].find('D') < 0 and x.split('\t')[5].find('I') < 0, output[:-1])
	introns = getIntrons(junctionHits)
	if len(introns) > 0:
		exons = []
		#/// excluding 5' most intron beyond start pos
#		introns = filter(lambda x: exon1.start() < x[0] and x[1] < exon2.end(), introns)
		introns = filter(lambda x: (exon1.start() < x[0] <= exon1.end() + 100) and (exon2.start() - 100 <= x[1] < exon2.end()), introns)
		for i in xrange(len(introns)):
			if i == 0: nexon = exon(exon1.chr(), exon1.start(), introns[i][0]-1, exon1.sense())
			else: nexon = exon(exon1.chr(), introns[i-1][1]+1, introns[i][0]-1, exon1.sense())
			exons.append(nexon)
			if i == len(introns) - 1:
				nexon = exon(exon1.chr(), introns[i][1]+1, exon2.end(), exon1.sense())
				exons.append(nexon)
		#/// excluding 3' most intron beyond end pos
		if len(exons) > 0: return True, exons
		else: return False, [exon1, exon2]
	else: return False, [exon1, exon2]

def readingBed(bedFiles):
	bedFiles = bedFiles.split(',')
	junctions = dict()
	for bedfile in bedFiles:
		for line in open(bedfile, 'r'):
			if not line.startswith('#'):
				line = line.strip().split('\t')
				chrom = line[0]; start = int(line[1]); end = int(line[2])
				hits = int(line[4]); sense = line[5]; block = line[10]
				start = start+int(block.split(',')[0])+1 ##
				end = end-int(block.split(',')[1]) ##
				if not junctions.has_key(chrom): junctions[chrom] = dict()
				intron = exon(chrom, start, end, sense)
				if not junctions[chrom].has_key(str(intron)): junctions[chrom][str(intron)] = [intron, hits]
				else:
					junctions[chrom][str(intron)][1] += hits
	return junctions

def readingAnno(anno, type):
	file = open(anno)
	lines = file.readlines(); file.close()
	lines = filter(lambda x: x[0] != '#', lines)

	siteD = dict()
	for line in lines:
		line = line.split('\t')
		if type == 'cage':
			chr = line[0]; sense = line[5]; tagN = float(line[4])
			if sense == '+': site = int(line[6])+1
			elif sense == '-': site = int(line[7])
		else: #type == 'polya':
			chr = line[0]; sense = line[5]; tagN = float(line[3])
			if sense == '+': site = int(line[2])
			elif sense == '-': site = int(line[1])+1
#		if int(tagN) >= 2:
		if not (siteD.has_key(chr)): siteD[chr] = []
		siteD[chr].append((site, tagN, sense))
	return siteD

def getSeq(locus, f):
	if f[-1]=='/': f = f[:-1]
	chr = locus.chr()
	start = locus.start()-1; end = locus.end(); sense = '+'
	commandbase = '/home/sjpyo/new/packages/blat/nibFrag ' + f + '/' + chr + '.nib ' + str(start) + ' ' + str(end) + ' ' + str(sense) + ' stdout'
	p=Popen(commandbase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	b=p.stdout.read().split('\n')
	return ''.join(b[1:]).upper()

def getSeq_exon(exons, sense, f):
	seq = ''
	for exon in exons:
		seq += getSeq(exon, f)
	if sense == '-':
		seq = reverseComp(seq)
	return seq

def getSeq_dict(fa):
	seq_dict = dict(); flag = 0
	inFile = open(fa, 'r')
	while 1:
		line = inFile.readline()
		if not line:
			seq_dict[pre_chr] = ''.join(seq)
			break
		if line.startswith('>'):
			if flag > 0:
				seq_dict[pre_chr] = ''.join(seq)
			pre_chr = line[1:].strip().split(' ')[0]
			seq = []; flag = 1
		else:
			seq.append(line.strip().upper())
	return seq_dict

def getSeq_dictbyChr(fa, chr):
	seq_dict = dict(); flag = 0
	inFile = open(fa, 'r')
	while 1:
		line = inFile.readline()
		if not line:
			seq_dict[pre_chr] = ''.join(seq)
			break
		if line.startswith('>'):
			if flag > 0 and pre_chr in [chr]:
				seq_dict[pre_chr] = ''.join(seq)
				break
			pre_chr = line[1:].strip()
			seq = []; flag = 1
		else:
			seq.append(line.strip().upper())
	return seq_dict

def reverseString(seq):
        li = []
        for i in seq: li.append(i)
        li.reverse()
        return ''.join(li)

def reverseComp(seq):
        comp = string.maketrans('ATCG', 'TAGC')
        return reverseString(seq).translate(comp)

def get1stCluster(tags):
	cluster = []; ntags = []
	toppos, read = tags[0]
	cluster.append((int(toppos), int(read)))
	for tag in tags:
		if -10 <= int(tag[0]) - toppos <= 10:
			if not int(toppos) == tag[0]: cluster.append(tag)
		else: ntags.append(tag)
	return cluster, ntags

def cmp2(a1, a2): return cmp(a1[1], a2[1])
def get5CAGESite(firstExon, lastExon, flankingLen, exons, cSiteD, fpseqThreshold, type):
	start = -1; end = -1; sense = firstExon.sense(); cexons = []; cexons += exons
	if sense == '+':
		start = firstExon.start() - flankingLen; end = lastExon.start()
		if end == firstExon.start(): end = lastExon.end()
		if len(exons) > 0:
			cexons[0] = exon(firstExon.chr(), start, firstExon.end(), sense)
			if not len(exons) == 1: cexons = cexons[:-1]
		else: cexons.append(exon(firstExon.chr(), start, firstExon.end(), sense))
	else: #sense == '-':
		start = lastExon.end(); end = firstExon.end() + flankingLen
		if start == firstExon.end(): start = lastExon.start()
		if len(exons) > 0:
			cexons[-1] = exon(firstExon.chr(), firstExon.start(), end, sense)
			if not len(exons) == 1: cexons = cexons[1:]
		else: cexons.append(exon(firstExon.chr(), firstExon.start(), end, sense))

	FpTags = []
	for csite in cSiteD[firstExon.chr()]:
		for cexon in cexons:
			if (csite[0] > cexon.start()) and (csite[0] < cexon.end()) and (csite[2] == sense):
				FpTags.append(csite)

	fpTags = dict()
	for Tag in FpTags:
		newTag = 0
		if fpTags.has_key(Tag[0]):
			newTag = int(fpTags[Tag[0]]) + int(Tag[1])
			fpTags[Tag[0]] = newTag
		else: fpTags[Tag[0]] = Tag[1]
	fpTags = fpTags.items(); fpTags.sort(cmp2); fpTags.reverse()

	clusters = []
	while len(fpTags) > 0:
		(cluster, fpTags) = get1stCluster(fpTags)
		clusters.append(cluster)

	distance = dict(); cfpEnds = []
	# Representative csite reads should be more than 2
	clusters = filter(lambda x: x[0][1] >= 2, clusters)
	# Representative csite reads should be more than 10% of maximum reads
	if len(clusters) > 0: maxRead = clusters[0][0][1]
	clusters = filter(lambda x: x[0][1]/float(maxRead) >= fpseqThreshold, clusters)
	if type == 'all':
		totalRead = 0
		for cluster in clusters:
			for site in cluster: totalRead += site[1]
		nclusters = []
		for cluster in clusters:
			ncluster = []
			for nsite in cluster:
#				if float(nsite[1])/float(totalRead) >= 0.05: ncluster.append(nsite)
				if float(nsite[1])/float(totalRead) >= 0.1: ncluster.append(nsite)
			if len(ncluster) > 0: nclusters.append(ncluster)
		if len(nclusters) > 0:  fpEnds = map(lambda x: x[0][0], nclusters)
		else: fpEnds = []
	elif type == 'major':
		if len(clusters) > 0: fpEnds = [clusters[0][0][0]]
		else: fpEnds = []
	else: #type == 'near':
		nfpEnds = map(lambda x: x[0][0], clusters)
		if len(nfpEnds) > 0:
			dfpEnds = dict()
			if sense == '+': trxEnd = firstExon.start()
			elif sense == '-': trxEnd = firstExon.end()
			for nfpEnd in nfpEnds:
				dis = abs(trxEnd-nfpEnd)
				if not dfpEnds.has_key(dis): dfpEnds[dis] = []
				dfpEnds[dis].append(nfpEnd)
			dfpEnd = dfpEnds.keys(); dfpEnd.sort()
			fpEnds = dfpEnds[dfpEnd[0]]
		else: fpEnds = []
	if len(fpEnds) > 0:
		for cluster in clusters:
			if cluster[0][0] in fpEnds:
				for fsite in cluster: cfpEnds.append(fsite[0])
	return fpEnds, distance
#	return fpEnds, cfpEnds, distance

def get3IPSite(firstExon, lastExon, flankingLen, exons, cSiteD, tpseqThreshold, type):
	start = -1; end = -1; sense = lastExon.sense(); cexons = []; cexons += exons
	if sense == '+':
		start = firstExon.end(); end = lastExon.end() + flankingLen
		if start == lastExon.end(): start = firstExon.start()
		cexons[-1] = exon(lastExon.chr(), lastExon.start(), end, sense)
		if not len(exons) == 1: cexons = cexons[1:]
	else: #sense == '-':
		start = lastExon.start() - flankingLen; end = firstExon.start()
		if end == lastExon.start(): end = firstExon.end()
		cexons[0] = exon(lastExon.chr(), start, lastExon.end(), sense)
		if not len(exons) == 1: cexons = cexons[:-1]

	TpTags = []
        for csite in cSiteD[lastExon.chr()]:
                for cexon in cexons:
                        if (csite[0] > cexon.start()) and (csite[0] < cexon.end()) and (csite[2] == sense):
                                TpTags.append(csite)

        tpTags = dict()
        for Tag in TpTags:
                newTag = 0
                if tpTags.has_key(Tag[0]):
                        newTag = int(tpTags[Tag[0]]) + int(Tag[1])
                        tpTags[Tag[0]] = newTag
                else: tpTags[Tag[0]] = Tag[1]
        tpTags = tpTags.items(); tpTags.sort(cmp2); tpTags.reverse()

        clusters = []
        while len(tpTags) > 0:
                (cluster, tpTags) = get1stCluster(tpTags)
                clusters.append(cluster)

	distance = dict(); ctpEnds = []
	# Representative csite reads should be more than 2
	clusters = filter(lambda x: x[0][1] >= 2, clusters)
	# Representative csite reads should be more than 10% of maximum reads
	if len(clusters) > 0: maxRead = clusters[0][0][1]
	clusters = filter(lambda x: x[0][1]/float(maxRead) >= tpseqThreshold, clusters)
	if type == 'all':
		totalRead = 0
		for cluster in clusters:
			for site in cluster: totalRead += site[1]
		nclusters = []
		for cluster in clusters:
			ncluster = []
			for nsite in cluster:
#				if float(nsite[1])/float(totalRead) >= 0.05: ncluster.append(nsite)
				if float(nsite[1])/float(totalRead) >= 0.1: ncluster.append(nsite)
			if len(ncluster) > 0: nclusters.append(ncluster)
		if len(nclusters) > 0: tpEnds = map(lambda x: x[0][0], nclusters)
		else: tpEnds = []
	elif type == 'major':
		if len(clusters) > 0: tpEnds = [clusters[0][0][0]]
		else: tpEnds = []
	else: #type == 'near':
		ntpEnds = map(lambda x: x[0][0], clusters)
		if len(ntpEnds) > 0:
			dtpEnds = dict()
			if sense == '+': trxEnd = lastExon.end()
			elif sense == '-': trxEnd =  lastExon.start()
			for ntpEnd in ntpEnds:
				dis = abs(trxEnd-ntpEnd)
				if not dtpEnds.has_key(dis): dtpEnds[dis] = []
				dtpEnds[dis].append(ntpEnd)
			dtpEnd = dtpEnds.keys(); dtpEnd.sort()
			tpEnds = dtpEnds[dtpEnd[0]]
		else: tpEnds = []
	if len(tpEnds) > 0:
		for cluster in clusters:
			if cluster[0][0] in tpEnds:
				for tsite in cluster: ctpEnds.append(tsite[0])
	return tpEnds, distance
#	return tpEnds, ctpEnds, distance

def filterNoneTrxs(trxs):
        ntrxs = []
        for trx in trxs:
                exons = trx.exons()
                flag = True
                for exon in exons:
                        if exon.len() <= 2: flag = False; break ##check
                if flag: ntrxs.append(trx)
        return ntrxs

def overlappedTrx(trx1, trx2): #checking trx1 -> trx2
        if trx1.sense() == trx2.sense() and (min(trx2.end(), trx1.end()) - max(trx1.start(), trx2.start())) / (float)(trx1.len()) >= 1 and (trx2.start() != trx1.start() or trx2.end() != trx1.end()): return True

def checkingProperTrx(trxs):
        ntrxs = []
        for trx in trxs:
                exons = trx.exons()
                if len(exons) == 1:
                        if exons[0].start() < exons[0].end(): ntrxs.append(trx)
                else:
                        preEnd = -1; flag = True
                        for exon in exons:
                                if exon.start() <= preEnd: flag = False; break
                                elif exon.start() >= exon.end(): flag = False; break
                                preEnd = exon.end()
                        if flag: ntrxs.append(trx)
        return ntrxs

def makingProperTrx(trxs):
	for trx in trxs:
		exons = trx.exons()
		for exon in exons:
			preStart = exon.start(); preEnd = exon.end()
			exon.newStart(min(preStart, preEnd)); exon.newEnd(max(preStart, preEnd))
	return trxs

def makeIndex(trx):
        return trx.chr() + '_' + trx.sense() + '_' + '_'.join(map(lambda x: str(x.start()), trx.exons())) + '_' + '_'.join(map(lambda x: str(x.end()), trx.exons()))

def makeIndex2(trx):
	return trx.chr() + '_' + trx.sense() + '_' + '_'.join(map(lambda x: str(x.start()), trx.introns())) + '_' + '_'.join(map(lambda x: str(x.end()), trx.introns()))

def filterSameTrxs(trxs):
	ntrxs = dict(); geneNum = len(trxs)
	for i in xrange(geneNum):
		trx = trxs[i]; index = makeIndex(trx)
		if not (ntrxs.has_key(index)): ntrxs[index] = []
#		else: print 'SameTrx: ', trx.trxid(), trx.coord()
		ntrxs[index].append(trx)
	ntrxs2 = []
	for trxA in ntrxs.values(): ntrxs2.append(trxA[0])
#	for trxA in ntrxs.values():
#		ntrxs2.append(trxA[0])
#		if len(trxA) > 1: print trxA[0].coord()
#	print 'filterTranscriptNum: ', geneNum, ' ', len(ntrxs2)
	return ntrxs2

def filterIncluTrxs(trxs):
	ntrxs = []; trxNum = len(trxs)
	for i in xrange(len(trxs)):
		ntrx = trxs[i]; exons = ntrx.exons()
#		trxsToBeChecked = trxs[max(0, i-20):i]
#		trxsToBeChecked += trxs[i+1:min(trxNum, i+20)]
		trxsToBeChecked = trxs[max(0, i-50):i]
		trxsToBeChecked += trxs[i+1:min(trxNum, i+50)]
		ans = False
		for trx in trxsToBeChecked:
			exons2 = trx.exons()
			# the parent trx should have one more exon to avoid messup
#			if (len(exons2) >= len(exons)): ##
			if (len(exons2) > len(exons)): ##
				inc = 0
				for ex in exons:
					for ex2 in exons2:
						if ex2.inclusive(ex): inc += 1
				if inc == len(exons): ans = True; break
		if not ans: ntrxs.append(ntrx)
#		else: print 'IncluTrx: ', ntrx.trxid(), ntrx.coord()
	return ntrxs

def filledIntrons(trxs):
	ntrxs = []; k = 1
	trxNum = len(trxs)
	for i in xrange(trxNum):
		trx = trxs[i]; exons = trx.exons()
#		if not len(exons)-1 == len(trx.introns()):
#			trx = transcript(trx.chr(), trx.geneid(), trx.trxid(), trx.start(), trx.end(), trx.exons(), trx.sense())
		if len(exons) > 1:
			intronLen = map(lambda x: x.len(), trx.introns())
			if min(intronLen) <= 60: ##
#				print trx.coord()
				nexons = []; nexons.append(exons[0])
				for x in xrange(1, len(exons)):
					if 0 < exons[x].start() - nexons[-1].end() <= 60: ##
						nexon = exon(nexons[-1].chr(), nexons[-1].start(), exons[x].end(), nexons[-1].sense())
						nexons[-1] = nexon
					else: nexons.append(exons[x])
				ntrx = transcript(trx.chr(), "Big-Gene-" + str(k), "Big-Trx-" + str(k) + '.1', trx.start(), trx.end(), nexons, trx.sense())
				ntrxs.append(ntrx); k+=1
			else:
				trx.setGeneid("Big-Gene-" + str(k)); trx.setTrxid("Big-Trx-" + str(k) + '.1')
				ntrxs.append(trx); k+=1
		else:
			trx.setGeneid("Big-Gene-" + str(k)); trx.setTrxid("Big-Trx-" + str(k) + '.1')
			ntrxs.append(trx); k+=1
	return ntrxs

def clusterOverlappedTrx(trxs):
	newGtf = dict(); loci = []
	trxNum = len(trxs); k = 1
	for i in xrange(trxNum):
		trx = trxs[i]
		if i == 0:
			newGtf['loci.' + str(k)] = []
			newGtf['loci.' + str(k)].append(trx)
			loci.append('loci.' + str(k))
		else:
			exons = trx.exons(); flag = False; cflag = False
			for x in xrange(len(loci)):
				for trx2 in newGtf[loci[x]]:
					if not trx.start() > trx2.end(): cflag = True
					exons2 = trx2.exons()
					for exon in exons:
						for exon2 in exons2:
							if exon.overlaps(exon2): flag = True
			if not flag and not cflag: loci = []
			if flag: newGtf[loci[x]].append(trx)
			else:
				k += 1
				newGtf['loci.' + str(k)] = []
				newGtf['loci.' + str(k)].append(trx)
				loci.append('loci.' + str(k))
	return newGtf

def renameTrx(loci, trxs):
	trxNum = len(trxs)
	for i in xrange(trxNum):
		trx = trxs[i]
		if trx.trxid().find('_TSS') > 0: tmp = '_TSS' + trx.trxid().split('_TSS')[-1]
		elif trx.trxid().find('_CPS') > 0: tmp = '_CPS' + trx.trxid().split('_CPS')[-1]
		else: tmp = ''
		trx.setGeneid("Big-Gene-Chr" + trx.chr()[3:] + '_' + str(loci.split('.')[-1]))
		trx.setTrxid("Big-Trx-Chr" + trx.chr()[3:] + '_' + str(loci.split('.')[-1]) + '.' + str(i+1) + tmp)
	return trxs

def renameTrx2(loci, trxs):
        trxNum = len(trxs)
        for i in xrange(trxNum):
                trx = trxs[i]
                trx.setGeneid("Gencode-Gene-Chr" + trx.chr()[3:] + '_' + str(loci.split('.')[-1]))
                trx.setTrxid("Gencode-Trx-Chr" + trx.chr()[3:] + '_' + str(loci.split('.')[-1]) + '.' + str(i+1))
        return trxs

def nf(inputfilename):
        ninputfilename = inputfilename + 'n'
        commandBase = 'awk NF ' + inputfilename + ' > ' + ninputfilename
        print commandBase
        p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        p.wait()
        os.remove(inputfilename)
        os.renames(ninputfilename, inputfilename)

def echo(inputfilename):
        commandBase = 'echo >> ' + inputfilename
        print commandBase
        p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        p.wait()

def cat(inputfilelist, outputfilename):
        commandBase = 'cat ' + ' '.join(inputfilelist) + ' >> ' + outputfilename
        print commandBase
        p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        p.wait()

def cufflinks(database, bamfilename, outputdir):
	if database == 'hg19': intronmin = 61; intrommax = 265006
	elif database == 'mm9': intronmin = 52; intrommax = 240764
	commandBase = 'cufflinks -L BIG_cNp --min-intron-length ' + str(intronmin) + ' --max-intron-length ' + str(intrommax) + ' --library-type fr-secondstrand --no-update-check -o ' + outputdir + ' ' + bamfilename
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()

def cuffcompare(database, inputfilename, outputdir, option="default"):
	if database == 'hg19':
		protein = data.hm_gencode_protein#; lncRNA = hm_gencode_lincRNA
#		protein_rpkm_1 = data.hela_protein_rpkm_1; lncRNA_rpkm_01 = data.hela_lncRNA_rpkm_01
#	elif database == 'mm9':
#		protein_rpkm_1 = data.mes_protein_rpkm_1; lncRNA_rpkm_01 = data.mes_lncRNA_rpkm_01
	outprefix = outputdir + '/protein'#; outprefix2 = outputdir + '/lncRNA/'
	commandBase = 'bsub -J cuffcompare cuffcompare -T -r ' + protein
#	commandBase = 'bsub -J cuffcompare cuffcompare -r ' + protein_rpkm_1
#	commandBase2 = 'bsub -J cuffcompare cuffcompare -r ' + lncRNA_rpkm_01
	if option == 'R': commandBase += ' -R'; outprefix += 'R'
	elif option == 'Q': commandBase += ' -Q'; outprefix += 'Q'
#	if option == 'R': commandBase += ' -R'; commandBase2 += ' -R '; outprefix += 'R'; outprefix2 += 'R'
#	elif option == 'Q': commandBase += ' -Q'; commandBase2 += ' -Q '; outprefix += 'Q'; outprefix2 += 'Q'
	commandBase += ' -o ' + outprefix + ' ' + inputfilename
#	commandBase2 += ' -o ' + outprefix2 + ' ' + inputfilename
	print commandBase
#	print commandBase, '\n', commandBase2
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
#	p2=Popen(commandBase2, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()#; p2.wait()
