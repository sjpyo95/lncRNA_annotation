##### SJ_module.py
##### "SJ module - variables, classes, functions"
##### Made by Sangjin Lee
##### Last modified: 2020-09-02
import os, sys
from subprocess import Popen, PIPE, STDOUT
from collections import defaultdict
import random
import time
import re

#***************  variables  ***************

hm_chrList = list(map(lambda x: 'chr'+x, list(map(str, range(1,23)))+['X','Y','M']))
ms_chrList = list(map(lambda x: 'chr'+x, list(map(str, range(1,20)))+['X','Y']))

annoM24 = '/lustre/export/home/sjlee/data/ref_seq/mouse/M24/gencode.vM24.annotation.gtf'
faM24 = '/lustre/export/home/sjlee/data/ref_seq/mouse/M24/GRCm38.p6.genome.fa'
annoM24pc = '/home/sjpyo/new/data/anno/mouse/gencode.vM24.protein_coding_genes.gtf'
annoM24lnc = '/home/sjpyo/new/data/anno/mouse/gencode.vM24.long_noncoding_RNAs.gtf'


#***************  classes  ***************

def chrOrderKey(chr):
    if not chr.startswith(('chr','Chr','CHR')): chr = 'chr' + chr
    if chr.lower() == 'chrx': chr = 'chr100'
    elif chr.lower() == 'chry': chr = 'chr101'
    elif chr.lower() == 'chrm': chr = 'chr102'
    return int(chr.strip('chr'))

class exon():
    def __init__(self, gene_id, trx_id, exon_number, chr, start, end, strand,
                 gene_type='', gene_name='', attribute='', source='.', score='.', frame='.'):
        self.__gene_id = gene_id
        self.__trx_id = trx_id
        self.__exon_number = exon_number
        if not chr.startswith('chr'): chr = 'chr' + chr
        self.__chr = chr
        self.__start = int(start)
        self.__end = int(end)
        self.__strand = strand
        self.__gene_type = gene_type
        self.__gene_name = gene_name
        self.__attribute = attribute
        self.__source = source
        self.__score = score
        self.__frame = frame
    def geneID(self): return self.__gene_id
    def trxID(self): return self.__trx_id
    def exonNum(self):
        if self.__exon_number != '.': return (int(self.__exon_number))
        return self.__exon_number
    def chr(self): return self.__chr
    def start(self):
        if self.__start != '.': return (int(self.__start))
        return self.__start
    def end(self):
        if self.__end != '.': return (int(self.__end))
        return self.__end
    def strand(self): return self.__strand
    def geneType(self): return self.__gene_type
    def geneName(self): return self.__gene_name
    def attribute(self): return self.__attribute
    def source(self): return self.__source
    def score(self): return self.__score
    def frame(self): return self.__frame
    def set_geneID(self, newGene_id): self.__gene_id = newGene_id
    def set_trxID(self, newTrx_id): self.__trx_id = newTrx_id
    def set_exonNum(self, newExon_number): self.__exon_number = newExon_number
    def set_chr(self, newChr): self.__chr = newChr
    def set_start(self, newStart): self.__start = newStart
    def set_end(self, newEnd): self.__end = newEnd
    def set_strand(self, newStrand): self.__strand = newStrand
    def set_geneType(self, newGene_type): self.__gene_type = newGene_type
    def set_geneName(self, newGeneName): self.__gene_name = newGeneName
    def set_attribute(self, newAttribute): self.__attribute = newAttribute
    def set_source(self, newSource): self.__source = newSource
    def set_score(self, newScore): self.__score = newScore
    def set_frame(self, newFrame): self.__frame = newFrame
    def update(self, gene_id='.', trx_id='.', exon_number='.', chr='.', start='.', end='.', strand='.',
               gene_type='.', gene_name='.', attribute='.', source='.', score='.', frame='.'):
        if gene_id != '.': self.set_geneID(gene_id)
        if trx_id != '.': self.set_trxID(trx_id)
        if exon_number != '.': self.set_exonNum(exon_number)
        if chr != '.': self.set_chr(chr)
        if start != '.': self.set_start(start)
        if end != '.': self.set_end(end)
        if strand != '.': self.set_strand(strand)
        if gene_type != '.': self.set_geneType(gene_type)
        if gene_name != '.': self.set_geneName(gene_name)
        if attribute != '.': self.set_attribute(attribute)
        if source != '.': self.set_source(source)
        if score != '.': self.set_score(score)
        if frame != '.': self.set_frame(frame)
    def len(self): return abs(self.__end - self.__start)
    def coords(self):
        return
    def flipStrand(self):
        if self.__strand == '+': self.__strand = '-'
        elif self.__strand == '-': self.__strand = '+'
    def strandedOverlap(self, otherExon):
        overlap = False
        if self.chr() != otherExon.chr(): return False
        if self.strand() == '.' or otherExon.strand() == '.': return False
        if self.strand() != otherExon.strand(): return False
        if otherExon.start() <= self.start() and self.start() < otherExon.end(): overlap = True
        if otherExon.start() < self.end() and self.end() <= otherExon.end(): overlap = True
        return overlap
    def unstrandedOverlap(self, otherExon):
        overlap = False
        if self.chr() != otherExon.chr(): return False
        # if not (self.strand() == '.' and otherExon.strand() == '.'): return False
        # if otherExon.start() <= self.start() and self.start() < otherExon.end(): overlap += 1
        # if otherExon.start() < self.end() and self.end() <= otherExon.end(): overlap += 1
        if otherExon.start() <= self.start() and self.start() < otherExon.end(): overlap = True
        if otherExon.start() < self.end() and self.end() <= otherExon.end(): overlap = True
        return overlap
    def isidentical(self, otherExon):
        identical = False
        selfExon = [self.chr(), self.start(), self.end(), self.strand()]
        otherExon = [otherExon.chr(), otherExon.start(), otherExon.end(), otherExon.strand()]
        if selfExon == otherExon:
            identical = True
        return identical

    def PSJ_inclusive(self, otherExon):
        if self.strand() == '.':
            if self.chr() == otherExon.chr() and self.start() <= otherExon.start() and self.end() >= otherExon.end(): return True
        else:
            if self.chr() == otherExon.chr() and self.strand() == otherExon.strand() and self.start() <= otherExon.start() and self.end() >= otherExon.end(): return True


def exonOrderKey_num(exon):
    return exon.exonNum(), exon.start(), exon.end()

def exonOrderKey_coords(exon):
    chr = exon.chr()
    if chr.lower() == 'chrx': chr = 'chr100'
    elif chr.lower() == 'chry': chr = 'chr101'
    elif chr.lower() == 'chrm': chr = 'chr102'
    chrOrder = int(chr.strip('chr'))
    if exon.strand() == '+': strand = 0
    elif exon.strand() == '-': strand = 1
    elif exon.strand() == '.': strand = 2
    else: strand = 3
    return chrOrder, exon.start(), exon.end(), strand

class transcript():
    def __init__(self, gene_id, trx_id, chr, strand, exons,
                 gene_type='', gene_name='', attribute='', cpc_score='', cpat_score='', gtf_count='', tss='', cps='',
                 source='.', score='.', frame='.'):
        self.__gene_id = gene_id
        self.__trx_id = trx_id
        if not chr.startswith('chr'): chr = 'chr' + chr
        self.__chr = chr
        self.__strand = strand
        self.__gene_type = gene_type
        self.__gene_name = gene_name
        self.__attribute = attribute
        self.__cpc_score = cpc_score
        self.__cpat_score = cpat_score
        self.__source = source
        self.__score = score
        self.__frame = frame
        exons = sorted(exons, key=exonOrderKey_coords)
        exonNum = 0
        for exon in exons:
            exonNum += 1
            exon.update(gene_id=self.__gene_id,trx_id=self.__trx_id,exon_number=exonNum,strand=self.__strand,
                        gene_type=self.__gene_type,gene_name=self.__gene_name,attribute=self.__attribute)
        self.__exons = exons
        self.set_introns()
        self.__start = min(map(lambda x: int(x.start()), exons))
        self.__end = max(map(lambda x: int(x.end()), exons))
        self.__gtf_count = gtf_count
        self.__skip = 0
        tssPattern = 'TSS\.\d+'
        tssMatches = re.findall(tssPattern, self.__trx_id)
        if len(tssMatches) >= 1:
            self.__tss = tssMatches[0]
        else:
            self.__tss = tss
            if self.__tss != '' and not self.__tss.startswith('TSS.'): self.__tss = 'TSS.' + self.__tss
        cpsPattern = 'CPS\.\d+'
        cpsMatches = re.findall(cpsPattern, self.__trx_id)
        if len(cpsMatches) >= 1:
            self.__cps = cpsMatches[0]
        else:
            self.__cps = cps
            if self.__cps != '' and not self.__cps.startswith('CPS.'): self.__cps = 'CPS.' + self.__cps
        # self.__use = False
        # self.__exj = False
        # self.__fpPos = -1
        # self.__tpPos = -1
        # self.__overlap = 0
    def set_introns(self):
        self.__introns = []
        exons = sorted(self.exons(), key=exonOrderKey_coords)
        for i in range(len(exons)-1):
            intron = exon(gene_id=exons[i].geneID(), trx_id=exons[i].trxID(), exon_number=exons[i].exonNum(),
                     chr=exons[i].chr(), start=exons[i].end()+1, end=exons[i + 1].start()-1,
                     strand=exons[i].strand(), gene_type=exons[i].geneType(), gene_name=exons[i].geneName(),
                     attribute=exons[i].attribute, source=exons[i].source(), score=exons[i].score(),
                     frame=exons[i].frame())
            self.__introns.append(intron)
    def exons(self): return self.__exons
    def geneID(self): return self.__gene_id
    def trxID(self): return self.__trx_id
    def chr(self): return self.__chr
    def start(self):
        if self.__start != '.': return (int(self.__start))
        return self.__start
    def end(self):
        if self.__end != '.': return (int(self.__end))
        return self.__end
    def strand(self): return self.__strand
    def geneType(self): return self.__gene_type
    def geneName(self): return self.__gene_name
    def attribute(self): return self.__attribute
    def introns(self): return self.__introns
    def cpcScore(self):
        if self.__cpc_score != '': return(float(self.__cpc_score))
        return self.__cpc_score
    def cpatScore(self):
        if self.__cpat_score != '': return (float(self.__cpat_score))
        return self.__cpat_score
    def source(self): return self.__source
    def score(self): return self.__score
    def frame(self): return self.__frame
    def tss(self): return self.__tss
    def cps(self): return self.__cps
    def gtfCount(self): return self.__gtf_count
    def skip(self): return self.__skip
    def set_exons(self, newExons):
        self.__exons = newExons
        self.set_introns()
    def set_geneID(self, newGene_id):
        self.__gene_id = newGene_id
        exons = sorted(self.exons(), key=exonOrderKey_num)
        for exon in exons:
            exon.set_geneID(newGene_id)
        self.set_exons(exons)
    def set_trxID(self, newTrx_id):
        self.__trx_id = newTrx_id
        exons = sorted(self.exons(), key=exonOrderKey_num)
        for exon in exons:
            exon.set_trxID(newTrx_id)
        self.set_exons(exons)
    def set_chr(self, newChr):
        self.__chr = newChr
        exons = sorted(self.exons(), key=exonOrderKey_num)
        for exon in exons:
            exon.set_chr(newChr)
        self.set_exons(exons)
    def set_start(self, newStart):
        self.__start = newStart
    def set_end(self, newEnd):
        self.__end = newEnd
    def set_strand(self, newStrand):
        self.__strand = newStrand
        exons = sorted(self.exons(), key=exonOrderKey_num)
        for exon in exons:
            exon.set_strand(newStrand)
        self.set_exons(exons)
    def set_geneType(self, newGene_type):
        self.__gene_type = newGene_type
        exons = sorted(self.exons(), key=exonOrderKey_num)
        for exon in exons:
            exon.set_geneType(newGene_type)
        self.set_exons(exons)
    def set_geneName(self, newGeneName):
        self.__gene_name = newGeneName
        exons = sorted(self.exons(), key=exonOrderKey_num)
        for exon in exons:
            exon.set_geneName(newGeneName)
        self.set_exons(exons)
    def set_attribute(self, newAttribute):
        self.__attribute = newAttribute
        exons = sorted(self.exons(), key=exonOrderKey_num)
        for exon in exons:
            exon.set_attribute(newAttribute)
        self.set_exons(exons)
    def set_cpcScore(self, newCpcScore):
        self.__cpc_score = newCpcScore
    def set_cpatScore(self, newCpatScore):
        self.__cpat_score = newCpatScore
    def set_source(self, newSource):
        self.__source = newSource
        exons = sorted(self.exons(), key=exonOrderKey_num)
        for exon in exons:
            exon.set_source(newSource)
        self.set_exons(exons)
    def set_score(self, newScore):
        self.__score = newScore
        exons = sorted(self.exons(), key=exonOrderKey_num)
        for exon in exons:
            exon.set_score(newScore)
        self.set_exons(exons)
    def set_frame(self, newFrame):
        self.__frame = newFrame
        exons = sorted(self.exons(), key=exonOrderKey_num)
        for exon in exons:
            exon.set_frame(newFrame)
        self.set_exons(exons)
    def set_tss(self, newTSS):
        if newTSS != '' and not newTSS.startswith('TSS.'): newTSS = 'TSS.' + newTSS
        self.__tss = newTSS
    def set_cps(self, newCPS):
        if newCPS != '' and not newCPS.startswith('CPS.'): newCPS = 'TSS.' + newCPS
        self.__cps = newCPS
    def update(self, gene_id='.', trx_id='.', chr='.', strand='.',
               gene_type='.', gene_name='.', attribute='.', cpc_score='.', cpat_score='.',
               source='.', score='.', frame='.', exons='.'):
        if gene_id != '.': self.set_geneID(gene_id)
        if trx_id != '.': self.set_trxID(trx_id)
        if chr != '.': self.set_chr(chr)
        if strand != '.': self.set_strand(strand)
        if gene_type != '.': self.set_geneType(gene_type)
        if gene_name != '.': self.set_geneName(gene_name)
        if attribute != '.': self.set_attribute(attribute)
        if cpc_score != '.': self.set_cpcScore(cpc_score)
        if cpat_score != '.': self.set_cpatScore(cpat_score)
        if source != '.': self.set_source(source)
        if score != '.': self.set_score(score)
        if frame != '.': self.set_frame(frame)
        if exons != '.' and isinstance(exons, list):
            if len(exons) > 0:
                self.set_tss('')
                self.set_cps('')
                if '_TSS.' in self.trxID():
                    self.set_trxID(self.trxID().split('_TSS.')[0])
                if '_CPS.' in self.trxID():
                    self.set_trxID(self.trxID().split('_CPS.')[0])
                print("TSS, CPS reset to ''")
                exons = sorted(exons, key=exonOrderKey_coords)
                exNum = 0
                for exon in exons:
                    exNum += 1
                    exon.update(gene_id=self.geneID(),trx_id=self.trxID(),exon_number=exNum,
                                gene_type=self.geneType(),gene_name=self.geneName(),attribute=self.attribute())
                self.set_start(min(map(lambda x: x.start(), exons)))
                self.set_end(max(map(lambda x: x.end(), exons)))
                self.set_exons(exons)
            else:
                raise Exception("'exons' given as an empty list")
    def len(self): return abs(self.__end - self.__start)
    def coords(self):
        return
    def exonCount(self): return len(self.__exons)
    def exonLen(self): return list(map(lambda x: x.len(), self.__exons))
    def exonLenTotal(self): return sum(map(lambda x: x.len(), self.__exons))
    def strandedOverlaps(self, otherTrx):
        overlaps = 0
        if self.chr() != otherTrx.chr(): return False
        if self.strand() == '.' or otherTrx.strand() == '.': return False
        if self.strand() != otherTrx.strand(): return False
        selfExons = self.exons()
        otherExons = otherTrx.exons()
        for selfExon in selfExons:
            hit = 0
            for otherExon in otherExons:
                x = 0
                if selfExon.strandedOverlap(otherExon):
                    x = 1
                    overlaps += 1
                if x == 1: hit = x
                if hit > x: break
        return overlaps
    def isidentical(self, otherTrx):
        identical = False
        selfExons = list(map(lambda x: [x.chr(), x.start(), x.end(), x.strand()], sorted(self.exons(), key=exonOrderKey_coords)))
        otherExons = list(map(lambda x: [x.chr(), x.start(), x.end(), x.strand()], sorted(otherTrx.exons(), key=exonOrderKey_coords)))
        if selfExons == otherExons:
            identical = True
        return identical
    def add_gtfCount(self, count=1):
        if self.__gtf_count == '': self.__gtf_count = 1
        if not isinstance(self.__gtf_count, int): self.__gtf_count = int(self.__gtf_count)
        self.__gtf_count += count
    def reset_gtfCount(self): self.__gtf_count = 1
    def add_skip(self):
        self.__skip = 1

def trxOrderKey_num(trx):
    idPattern = '[G\.\_]\d+'
    geneID = trx.geneID()
    if geneID == '':
        raise Exception('Check transcript gene IDs')
    geneID_nums = re.findall(idPattern, geneID)
    if len(geneID_nums) == 1:
        geneID_num = int(geneID_nums[0].strip('._').lstrip('0'))
        geneID_prefix = geneID.split(geneID_nums[0])
    elif len(geneID_nums) == 0:
        geneID_num = 0
        geneID_prefix = geneID
    else:
        raise Exception('Check gene IDs')
    floatPattern = '\d+\.\d+'
    trxID = trx.trxID()
    trxID_num = re.findall(floatPattern, trxID)[0]
    isoform_num = int(trxID_num.split('.')[1])
    chr = trx.chr()
    if chr.lower() == 'chrx': chr = 'chr100'
    elif chr.lower() == 'chry': chr = 'chr101'
    elif chr.lower() == 'chrm': chr = 'chr102'
    chrOrder = int(chr.strip('chr'))
    if trx.strand() == '+': strand = 0
    elif trx.strand() == '-': strand = 1
    elif trx.strand() == '.': strand = 2
    else: strand = 3
    gtfCount = trx.gtfCount()
    if trx.gtfCount() == '': gtfCount = 0
    gtfCount = - gtfCount
    return chrOrder, geneID_prefix, geneID_num, isoform_num, trx.start(), trx.end(), strand, -trx.exonCount(), -trx.exonLenTotal(), gtfCount

def trxOrderKey_coords(trx):
    chr = trx.chr()
    if chr.lower() == 'chrx': chr = 'chr100'
    elif chr.lower() == 'chry': chr = 'chr101'
    elif chr.lower() == 'chrm': chr = 'chr102'
    chrOrder = int(chr.strip('chr'))
    if trx.strand() == '+': strand = 0
    elif trx.strand() == '-': strand = 1
    elif trx.strand() == '.': strand = 2
    else: strand = 3
    gtfCount = trx.gtfCount()
    if trx.gtfCount() == '': gtfCount = 0
    gtfCount = - gtfCount
    return chrOrder, trx.start(), trx.end(), strand, -trx.exonCount(), -trx.exonLenTotal(), gtfCount

class gene():
    def __init__(self, gene_id, chr, strand, trxs,
                 gene_type='', gene_name='', attribute='', start='.', end='.', source='.', score='.', frame='.'):
        self.__gene_id = gene_id
        if not chr.startswith('chr'): chr = 'chr'+chr
        self.__chr = chr
        self.__strand = strand
        self.__gene_type = gene_type
        self.__gene_name = gene_name
        self.__attribute = attribute
        self.__source = source
        self.__score = score
        self.__frame = frame
        if trxs == '.': self.__trxs = []
        if len(trxs) > 0:
            trxs = sorted(trxs, key=trxOrderKey_coords)
            trxNum = 0
            for trx in trxs:
                trxNum += 1
                trxID = self.__gene_id + '.' + str(trxNum)
                trx.update(gene_id=self.__gene_id, trx_id=trxID,
                           gene_type=self.geneType(), gene_name=self.geneName(), attribute=self.attribute())
            self.__start = min(map(lambda x: int(x.start()), trxs))
            self.__end = max(map(lambda x: int(x.end()), trxs))
        else:
            self.__start = start
            self.__end = end
        self.__trxs = trxs
    def trxs(self): return self.__trxs
    def geneID(self): return self.__gene_id
    def chr(self): return self.__chr
    def start(self):
        if self.__start != '.': return (int(self.__start))
        return self.__start
    def end(self):
        if self.__end != '.': return (int(self.__end))
        return self.__end
    def strand(self): return self.__strand
    def geneType(self): return self.__gene_type
    def geneName(self): return self.__gene_name
    def attribute(self): return self.__attribute
    def source(self): return self.__source
    def score(self): return self.__score
    def frame(self): return self.__frame
    def set_trxs(self, newTrxs):
        self.__trxs = newTrxs
    def set_geneID(self, newGene_id):
        self.__gene_id = newGene_id
        trxs = sorted(self.trxs(), key=trxOrderKey_coords)
        trxNum = 0
        for trx in trxs:
            trx.set_geneID(newGene_id)
            trxNum += 1
            trxID = newGene_id + '.' + str(trxNum)
            trx.set_trxID(trxID)
        self.set_trxs(trxs)
    def set_chr(self, newChr):
        self.__chr = newChr
        trxs = sorted(self.trxs(), key=trxOrderKey_num)
        for trx in trxs:
            trx.set_chr(newChr)
        self.set_trxs(trxs)
    def set_start(self, newStart):
        self.__start = newStart
    def set_end(self, newEnd):
        self.__end = newEnd
    def set_strand(self, newStrand):
        self.__strand = newStrand
        trxs = sorted(self.trxs(), key=trxOrderKey_num)
        for trx in trxs:
            trx.set_strand(newStrand)
        self.set_trxs(trxs)
    def set_geneType(self, newGene_type):
        self.__gene_type = newGene_type
        trxs = sorted(self.trxs(), key=trxOrderKey_num)
        for trx in trxs:
            trx.set_geneType(newGene_type)
        self.set_trxs(trxs)
    def set_geneName(self, newGeneName):
        self.__gene_name = newGeneName
        trxs = sorted(self.trxs(), key=trxOrderKey_num)
        for trx in trxs:
            trx.set_geneName(newGeneName)
        self.set_trxs(trxs)
    def set_attribute(self, newAttribute):
        self.__attribute = newAttribute
        trxs = sorted(self.trxs(), key=trxOrderKey_num)
        for trx in trxs:
            trx.set_attribute(newAttribute)
        self.set_trxs(trxs)
    def set_source(self, newSource):
        self.__source = newSource
        trxs = sorted(self.trxs(), key=trxOrderKey_num)
        for trx in trxs:
            trx.set_source(newSource)
        self.set_trxs(trxs)
    def set_score(self, newScore):
        self.__score = newScore
        trxs = sorted(self.trxs(), key=trxOrderKey_num)
        for trx in trxs:
            trx.set_score(newScore)
        self.set_trxs(trxs)
    def set_frame(self, newFrame):
        self.__frame = newFrame
        trxs = sorted(self.trxs(), key=trxOrderKey_num)
        for trx in trxs:
            trx.set_frame(newFrame)
        self.set_trxs(trxs)
    def update(self, gene_id='.', chr='.', strand='.', source='.', score='.', frame='.',
               gene_type='.', gene_name='.', attribute='.', trxs='.'):
        if gene_id != '.': self.set_geneID(gene_id)
        if chr != '.': self.set_chr(chr)
        if strand != '.': self.set_strand(strand)
        if gene_type != '.': self.set_geneType(gene_type)
        if gene_name != '.': self.set_geneName(gene_name)
        if attribute != '.': self.set_attribute(attribute)
        if source != '.': self.set_source(source)
        if score != '.': self.set_score(score)
        if frame != '.': self.set_frame(frame)
        if trxs != '.' and isinstance(trxs, list):
            if len(trxs) > 0:
                trxs = sorted(trxs, key=trxOrderKey_coords)
                trxNum = 0
                for trx in trxs:
                    tssPattern = 'TSS\.\d+'
                    tssMatches = re.findall(tssPattern, trx.trxID())
                    if len(tssMatches) >= 1:
                        trxTSS = '_' + tssMatches[0]
                    else:
                        trxTSS = ''
                    cpsPattern = 'CPS\.\d+'
                    cpsMatches = re.findall(cpsPattern, trx.trxID())
                    if len(cpsMatches) >= 1:
                        trxCPS = '_' + cpsMatches[0]
                    else:
                        trxCPS = ''
                    tag = '' + trxTSS + trxCPS
                    trxNum += 1
                    if not trx.geneID().startswith('Gene-') and not trx.trxID().startswith('Trx-'):
                        trxID = self.geneID() + '.' + str(trxNum) + tag
                    else:
                        trxID = trx.trxID()
                    trx.update(gene_id=self.geneID(),trx_id=trxID,
                                gene_type=self.geneType(),gene_name=self.geneName(),attribute=self.attribute())
                self.set_start(min(map(lambda x: int(x.start()), trxs)))
                self.set_end(max(map(lambda x: int(x.end()), trxs)))
                self.set_trxs(trxs)
            else:
                raise Exception("'trxs' given as an empty list")
    def len(self): return abs(int(self.__end) - int(self.__start))
    def trxCount(self): return len(self.__trxs)
    # def isidentical(self, otherGene):
    #     identical = False
    #     selfTrxs = list(map(lambda x: [x.chr(), x.start(), x.end(), x.strand(), list(map(lambda y: [y.chr(), y.start(), y.end(), y.strand()], sorted(x.exons(), key=exonOrderKey_coords))],
    #                         sorted(self.trxs(), key=trxOrderKey_num))))
    #     otherTrxs = list(map(lambda x: [x.chr(), x.start(), x.end(), x.strand(), list(map(lambda y: [y.chr(), y.start(), y.end(), y.strand()], sorted(x.exons(), key=exonOrderKey_coords))],
    #                          sorted(otherGene.trxs(), key=trxOrderKey_num))))
    #     if selfTrxs == otherTrxs:
    #         identical = True
    #     return identical

def geneOrderKey_num(gene):
    idPattern = '[G\.\_]\d+'
    geneID = gene.geneID()
    if geneID == '':
        raise Exception('Check transcript gene IDs')
    geneID_nums = re.findall(idPattern, geneID)
    if len(geneID_nums) == 1:
        geneID_num = int(geneID_nums[0].strip('._').lstrip('0'))
        geneID_prefix = geneID.split(geneID_nums[0])
    elif len(geneID_nums) == 0:
        geneID_num = 0
        geneID_prefix = geneID
    else:
        raise Exception('Check gene IDs')
    chr = gene.chr()
    if chr.lower() == 'chrx': chr = 'chr100'
    elif chr.lower() == 'chry': chr = 'chr101'
    elif chr.lower() == 'chrm': chr = 'chr102'
    chrOrder = int(chr.strip('chr'))
    if gene.strand() == '+': strand = 0
    elif gene.strand() == '-': strand = 1
    elif gene.strand() == '.': strand = 2
    else: strand = 3
    return chrOrder, geneID_prefix, geneID_num, gene.start(), gene.end(), -gene.trxCount(), strand

def geneOrderKey_coords(gene):
    chr = gene.chr()
    if chr.lower() == 'chrx': chr = 'chr100'
    elif chr.lower() == 'chry': chr = 'chr101'
    elif chr.lower() == 'chrm': chr = 'chr102'
    chrOrder = int(chr.strip('chr'))
    if gene.strand() == '+': strand = 0
    if gene.strand() == '-': strand = 1
    if gene.strand() == '.': strand = 2
    else: strand = 3
    return chrOrder, gene.start(), gene.end(), strand

def gtfOrderKey_gene(gtfLine):
    if gtfLine.count('\t') != 8:
        notHeader = 0
        chrOrder=0;geneID_prefix=0;geneID_num=0;isoform_num=0;feature=0;start=0;end=0;strand=0
    else:
        notHeader = 1
        lineList = gtfLine.split('\t')
        chr = lineList[0]; source = lineList[1]; feature = lineList[2]; start = lineList[3]; end = lineList[4]
        score = lineList[5]; strand = lineList[6]; frame = lineList[7]; attributes = lineList[8]
        geneID = attributes.split('gene_id')[1].split(';')[0].split('"')[1]
        if feature == 'transcript' or feature == 'exon':
            trxID = attributes.split('transcript_id')[1].split(';')[0].split('"')[1]
            floatPattern = '\d+\.\d+'
            trxID_num = re.findall(floatPattern, trxID)[0]
            isoform_num = int(trxID_num.split('.')[1])
            numStart = trxID.find(trxID_num)
            numPeriod = trxID_num.find('.')
            if geneID == '':
                if trxID[numStart + numPeriod] != '.':
                    raise Exception('Check transcript IDs')
                geneID = trxID[:numStart+numPeriod]
        else:
            isoform_num = 0
        idPattern = '[G\.\_]\d+'
        geneID_nums = re.findall(idPattern, geneID)
        if len(geneID_nums) == 1:
            geneID_num = int(geneID_nums[0].strip('._').lstrip('0'))

            geneID_prefix = geneID.split(geneID_nums[0])
        elif len(geneID_nums) == 0:
            geneID_num = 0
            geneID_prefix = geneID
        else:
            raise Exception('Check gene IDs')
        if chr.lower() == 'chrx': chr = 'chr100'
        elif chr.lower() == 'chry': chr = 'chr101'
        elif chr.lower() == 'chrm': chr = 'chr102'
        chrOrder = int(chr.strip('chr'))
        if feature == 'gene': feature = 0
        elif feature == 'transcript': feature = 1
        elif feature == 'exon': feature = 2
        else: feature = 3
    return notHeader, chrOrder, geneID_prefix, geneID_num, isoform_num, feature, start, end, strand

#***************  functions  ***************

def find_path(program):
    find = Popen(['which', program], stdout=PIPE, stderr=STDOUT)
    result, stderr = find.communicate()
    result = str(result)
    if ':' in str(result): print("couldn't find {} path: aborting".format(program)); sys.exit()
    else:
        result = result.strip('\n')
        print("{} path: '{}'".format(program, result))
        return result

def qsub_execute(job, command, nodes2avoid='.', log='.'):
    '''executes qsub'''
    if log == '.':
        log = job
    if log.strip('/').count('/') == 0:
        print('log directory not specified\n'+'generating log in current working directory')
        log = os.getcwd() + '/' + log
    if not log.endswith('.log'): log += '.log'
    node = str(random.choice([i for i in range(1,7) if str(i) not in nodes2avoid]))
    qsub = 'qsub -N {} -lvnode=biglab{} -j oe -o {} -- '.format(job,node,log)
    qsub += command
    print('**********{:^13}**********\n'.format('executing') + qsub+'\n' + '{}\n'.format('*'*33))
    # print('{}\n'.format('*'*29)+'********{:^13}********\n'.format('executing') + qsub + '{}\n{}\n'.format('*'*29,'*'*29))
    cmd = Popen(qsub, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    stdout, stderr = cmd.communicate()
    status = cmd.wait()
    if status != 0:
        print('There was an error')
        sys.exit()
    return stdout, stderr

def qsub_I_string(job, command, nodes2avoid='.', log='.'):
    '''creates qsub -I string to copy and run'''
    if log == '.':
        log = job
    if log.strip('/').count('/') == 0:
        print('log directory not specified\n'+'generating log in current working directory')
        log = os.getcwd() + '/' + log
    if not log.endswith('.log'): log += '.log'
    node = str(random.choice([i for i in range(1,7) if str(i) not in nodes2avoid]))
    qsub = 'qsub -I -N {} -lvnode=biglab{} -j oe -o {} -- '.format(job,node,log)
    qsub += command
    print('**********{:^20}**********\n'.format('string to execute') + qsub+'\n' + '{}\n'.format('*'*40))
    # return stdout, stderr

def qsub_wait(job, nJob=0, t=30, user='sjlee'):
    '''waits until number of certain job by certain user is under threshold (nJob, default=0)'''
    if len(job) > 15: job = job[:15]#; print("qsub 'job' too long: shortened to {} (15 char.)".format(job))
    spaces = ' '*(17-len(job))
    while True:
        cmd = Popen('qstat', shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        qstatOut, stderr = cmd.communicate()
        if qstatOut.count(job+spaces+user) > int(nJob): time.sleep(int(t))
        else: break

#------

def add_GTFattribute(gtfLine, attributeName, attributeValue='', order='.'):
    '''Adds an attribute to a single GTF line.
    If 'order' is given as an integer, 'attribute' is added as 'order'th attribute. Added to end if exceeds limit.
    If 'order' is given as another attribute, 'attribute' is added after that attribute. Added to end if no such attribute'''
    lineList = gtfLine.strip('\n').split('\t')
    Attributes = lineList[8].strip()
    AttributeList = Attributes.split(';')
    attribute = ' {} "{}"'.format(attributeName, attributeValue)
    if order == '.':
        if Attributes.find(attributeName) > -1:
            idx = next(i for i,s in enumerate(AttributeList) if s.split('"')[0].strip() == attributeName)
            AttributeList[idx] = attribute
        else:
            AttributeList.append(attribute)
    elif isinstance(order, int):
        if Attributes.find(attributeName) > -1:
            AttributeList = [s for s in AttributeList if s.split('"')[0].strip() != attributeName]
        AttributeList.insert(order-1, attribute)
    elif isinstance(order, str):
        if Attributes.find(attributeName) > -1:
            AttributeList = [s for s in AttributeList if s.split('"')[0].strip() != attributeName]
        if Attributes.find(order) > -1:
            idx = next(i for i,s in enumerate(AttributeList) if s.split('"')[0].strip() == order)
            AttributeList.insert(idx+1, attribute)
        else: #'order' not in Attributes
            AttributeList.append(attribute)
    lineList[8] = ';'.join(AttributeList).strip()
    newGtfLine = '\t'.join(lineList) + '\n'
    return newGtfLine

def filter_ncGenes(geneList, cpc_cutoff='.', cpat_cutoff='.'):
    if cpc_cutoff == '.': cpc_cutoff = -0.3
    else: cpc_cutoff = float(cpc_cutoff)
    if cpat_cutoff == '.': cpat_cutoff = 0.44
    else: cpat_cutoff = float(cpat_cutoff)
    ncGeneList = []
    for gene in geneList:
        trxs = gene.trxs()
        coding = False; cpc_blank = []; cpat_blank = []
        for trx in trxs:
            if coding: break
            if trx.cpcScore() == '': cpc_blank.append(trx)
            else:
                if trx.cpcScore() >= cpc_cutoff: coding = True
            if trx.cpatScore() == '': cpat_blank.append(trx)
            else:
                if trx.cpatScore() >= cpat_cutoff: coding = True
        if coding: continue
        if len(cpc_blank) == len(trxs) and len(cpat_blank) == len(trxs):
            print('gene with no cpc/cpat score: {}'.format(gene.geneID()))
            continue
        ncGeneList.append(gene)
    return ncGeneList

def get_commonTrxs(*GTFs):
    '''input type must be list'''
    num = len(GTFs)
    print('GTFs: {}'.format(GTFs))
    print('len(GTFs): {}'.format(len(GTFs)))
    allTrx = []
    commonTrxDic = {}
    i = 0
    for gtf in GTFs:
        i += 1
        iF = open(gtf, 'r')
        gtfLines = iF.readlines(); iF.close()
        # print('len(gtfLines): {}'.format(len(gtfLines)))
        trxLines = list(filter(lambda x: x.split('\t')[2] == 'transcript', gtfLines))
        # print('len(trxLines): {}'.format(len(trxLines)))
        trxList = list(map(lambda x: next(j for j in x.split('\t')[8].split(';') if 'transcript_id' in j).split('"')[1], trxLines))
        uniqTrxList = list(set(trxList))
        d = defaultdict(int)
        for trx in trxList:
            d[trx] += 1
        items = d.items()
        twoUpList = [x for x in items if x[1] >= 2]
        # lowerTrxList = list(map(lambda x: x.lower(), trxList))
        allTrx.extend(uniqTrxList)
        #stats
        print('gtf{}: {}'.format(i, gtf))
        print('number of total transcripts: {}'.format(len(trxList)))
        print('number of unique transcripts: {}'.format(len(uniqTrxList)))
        print('number of transcripts with 2 or more copies: {}'.format(len(twoUpList)))
        print('transcripts with 2 or more copies: {}\n'.format(twoUpList))
    for trx in allTrx:
        if allTrx.count(trx) == num:
            if trx not in commonTrxDic: commonTrxDic[trx] = ''
        else: continue
    commonTrxs = commonTrxDic.keys()
    print('number of common transcripts: {}\n'.format(len(commonTrxs)))
    return commonTrxs

def get_GTFgenes():
    return

def get_GTFtrxs(GTF):
    '''Returns a {key:[trxs]} dictionary'''
    chrList = list(map(lambda x: 'chr'+str(x), list(range(1, 23)) + ['X','Y','M']))
    iF = open(GTF, 'r')
    gtfLines = iF.readlines(); iF.close()
    chrTrxDic = {}
    for i in range(len(gtfLines)):
        line = gtfLines[i]
        if line.count('\t') != 8: continue
        lineList = line.strip('\n ').split('\t')
        chr = lineList[0]; feature = lineList[2]
        if not chr.startswith(('chr', 'Chr', 'CHR')): chr = 'chr'+chr
        if chr not in chrList: continue #only mouse/human chromosomes
        # if feature == 'gene': continue
        # elif feature == 'exon': continue
        if feature != 'transcript': continue
        else:
            trxSource = lineList[1]; trxStart = int(lineList[3]); trxEnd = int(lineList[4])
            trxScore = lineList[5]; trxStrand = lineList[6]; trxFrame = lineList[7]; trxAttributes = lineList[8]
            if trxAttributes.find('gene_id') > -1: trxGeneID = trxAttributes.split('gene_id')[1].split(';')[0].split('"')[1]
            else: print('trx with no gene_id: {}'.format(line)); continue
            if trxAttributes.find('transcript_id') > -1: trxTrxID = trxAttributes.split('transcript_id')[1].split(';')[0].split('"')[1]
            else: print('trx with no transcript_id: {}'.format(line)); continue
            if trxAttributes.find('gene_type') > -1: trxGeneType = trxAttributes.split('gene_type')[1].split(';')[0].split('"')[1]
            else: trxGeneType = ''
            if trxAttributes.find('gene_name') > -1: trxGeneName = trxAttributes.split('gene_name')[1].split(';')[0].split('"')[1]
            else: trxGeneName = ''
            if trxAttributes.find('attribute') > -1: trxAttribute = trxAttributes.split('attribute')[1].split(';')[0].split('"')[1]
            else: trxAttribute = ''
            if trxAttributes.find('cpc_score') > -1: trxCpcScore = trxAttributes.split('cpc_score')[1].split(';')[0].split('"')[1]
            else: trxCpcScore = ''
            if trxAttributes.find('cpat_score') > -1: trxCpatScore = trxAttributes.split('cpat_score')[1].split(';')[0].split('"')[1]
            else: trxCpatScore = ''
            if trxAttributes.find('gtf_count') > -1: trxGtfCount = int(trxAttributes.split('gtf_count')[1].split(';')[0].split('"')[1])
            else: trxGtfCount = ''
            tssPattern = 'TSS\.\d+'
            tssMatches = re.findall(tssPattern, trxAttributes)
            if len(tssMatches) >= 1: trxTSS = tssMatches[0]
            else: trxTSS = ''
            cpsPattern = 'CPS\.\d+'
            cpsMatches = re.findall(cpsPattern, trxAttributes)
            if len(cpsMatches) >= 1: trxCPS = cpsMatches[0]
            else: trxCPS = ''
            exons = []
            try:
                #transcript with most exons in 16 B cell samples merge had 146 exons -> check next 1000 lines to be safe
                if ((len(gtfLines)-1) - (i+1)) > 1000:
                    nextTrx = next(i+1+idx for idx, s in enumerate(gtfLines[i+1:i+1+1000]) if s.split('\t')[2] =='transcript')
                else: #elif (len(gtfLines)-1) <= 1000:
                    nextTrx = next(i+1+idx for idx, s in enumerate(gtfLines[i+1:]) if s.split('\t')[2] == 'transcript')
                trxExonLines = gtfLines[i+1:nextTrx]
                #in case transcripts with identical trx_ids messed up line order
                if len(trxExonLines) == 0:
                    print('{}\n{}'.format(gtfLines[i],gtfLines[i+1]))
                    try:
                        if ((len(gtfLines)-1) - (nextTrx)) > 1000:
                            nextnextTrx = next(nextTrx+idx for idx, s in enumerate(gtfLines[nextTrx:nextTrx+1000]) if s.split('\t')[2] =='transcript')
                        else:
                            nextnextTrx = next(nextTrx + idx for idx, s in enumerate(gtfLines[nextTrx:]) if s.split('\t')[2] == 'transcript')
                        trxExonLines = list(set(gtfLines[nextTrx:nextnextTrx]))
                    except StopIteration:
                        # last transcript
                        trxExonLines = gtfLines[nextTrx:]
            except StopIteration:
                # last transcript
                trxExonLines = gtfLines[i+1:]
            exNum = 0
            for j in range(len(trxExonLines)):
                exLine = trxExonLines[j]
                if exLine.count('\t') != 8: continue
                exLineList = exLine.split('\t')
                exChr = exLineList[0]; exFeature = exLineList[2]
                if exFeature != 'exon': continue
                exSource = exLineList[1]; exStart = int(exLineList[3]); exEnd = int(exLineList[4]); exScore = exLineList[5]
                exStrand = exLineList[6]; exFrame = exLineList[7]; exAttributes = exLineList[8]
                if exAttributes.find('gene_id') > -1: exGeneID = exAttributes.split('gene_id')[1].split(';')[0].split('"')[1]
                else: print('exon with no gene_id: {}'.format(line)); continue
                if exGeneID != trxGeneID: continue
                if exAttributes.find('transcript_id') > -1: exTrxID = exAttributes.split('transcript_id')[1].split(';')[0].split('"')[1]
                else: print('exon with no transcript_id: {}'.format(line)); continue
                if exTrxID != trxTrxID: continue
                if exAttributes.find('gene_type') > -1: exGeneType = exAttributes.split('gene_type')[1].split(';')[0].split('"')[1]
                else: exGeneType = ''
                if exAttributes.find('gene_name') > -1: exGeneName = exAttributes.split('gene_name')[1].split(';')[0].split('"')[1]
                else: exGeneName = ''
                if exAttributes.find('attribute') > -1: exAttribute = exAttributes.split('attribute')[1].split(';')[0].split('"')[1]
                else: exAttribute = ''
                exNum += 1
                nExon = exon(gene_id=exGeneID,trx_id=exTrxID,exon_number=exNum,chr=exChr,start=exStart,end=exEnd,
                             strand=exStrand,gene_type=exGeneType,gene_name=exGeneName,attribute=exAttribute,
                             source=exSource,score=exScore,frame=exFrame)
                exons.append(nExon)
            nTrx = transcript(gene_id=trxGeneID,trx_id=trxTrxID,chr=chr,strand=trxStrand,exons=exons,
                              gene_type=trxGeneType,gene_name=trxGeneName,attribute=trxAttribute,
                              cpc_score=trxCpcScore,cpat_score=trxCpatScore,gtf_count=trxGtfCount,
                              source=trxSource,score=trxScore,frame=trxFrame,tss=trxTSS,cps=trxCPS)
            if chr not in chrTrxDic: chrTrxDic[chr] = []
            chrTrxDic[chr].append(nTrx)
    return chrTrxDic

def get_headerLines(file):
    headerLines = []
    oF = open(file, 'r')
    nonHeader = 0
    for line in oF:
        if nonHeader >= 2: break
        if line.startswith('#'):
            headerLines.append(line)
            if nonHeader > 0:
                nonHeader -= 1
        else:
            nonHeader += 1
    oF.close()
    if len(headerLines) == 0: print('no header lines')
    return headerLines

def get_uniqTrxs(trxList, count='O'):
    '''Max. 5 replicates'''
    trxList = sorted(trxList, key=trxOrderKey_coords)
    nTrxList = []
    for i in range(len(trxList)-4):
        trx = trxList[i]
        next1Trx = trxList[i + 1]; next2Trx = trxList[i + 2]; next3Trx = trxList[i + 3]; next4Trx = trxList[i + 4]
        if i != len(trxList)-5:
            if trx.skip() == 1: continue
            if trx.isidentical(next1Trx):
                if count == 'O': trx.add_gtfCount()
                next1Trx.add_skip()
                if trx.isidentical(next2Trx):
                    if count == 'O': trx.add_gtfCount()
                    next2Trx.add_skip()
                    if trx.isidentical(next3Trx):
                        if count == 'O': trx.add_gtfCount()
                        next3Trx.add_skip()
                        if trx.isidentical(next4Trx):
                            if count == 'O': trx.add_gtfCount()
                            next4Trx.add_skip()
            nTrxList.append(trx)
        if i == len(trxList)-5:
            if trx.skip() != 1:
                if trx.isidentical(next1Trx):
                    if count == 'O': trx.add_gtfCount()
                    next1Trx.add_skip()
                    if trx.isidentical(next2Trx):
                        if count == 'O': trx.add_gtfCount()
                        next2Trx.add_skip()
                        if trx.isidentical(next3Trx):
                            if count == 'O': trx.add_gtfCount()
                            next3Trx.add_skip()
                            if trx.isidentical(next4Trx):
                                if count == 'O': trx.add_gtfCount()
                                next4Trx.add_skip()
                nTrxList.append(trx)
            if next1Trx.skip() != 1:
                if next1Trx.isidentical(next2Trx):
                    if count == 'O': next1Trx.add_gtfCount()
                    next2Trx.add_skip()
                    if next1Trx.isidentical(next3Trx):
                        if count == 'O': next1Trx.add_gtfCount()
                        next3Trx.add_skip()
                        if next1Trx.isidentical(next4Trx):
                            if count == 'O': next1Trx.add_gtfCount()
                            next4Trx.add_skip()
                nTrxList.append(next1Trx)
            if next2Trx.skip() != 1:
                if next2Trx.isidentical(next3Trx):
                    if count == 'O': next2Trx.add_gtfCount()
                    next3Trx.add_skip()
                    if next2Trx.isidentical(next4Trx):
                        if count == 'O': next2Trx.add_gtfCount()
                        next4Trx.add_skip()
                nTrxList.append(next2Trx)
            if next3Trx.skip() != 1:
                if next3Trx.isidentical(next4Trx):
                    if count == 'O': next3Trx.add_gtfCount()
                    next4Trx.add_skip()
                nTrxList.append(next3Trx)
            if next4Trx.skip() != 1:
                nTrxList.append(next4Trx)
    return nTrxList

# def get_uniqTrxs(trxList, count='O'):
#     '''Max. 5 replicates'''
#     trxList = sorted(trxList, key=trxOrderKey_coords)
#     nTrxList = []
#     for i in range(len(trxList)-4):
#         # if i == 150: sys.exit()
#         trx = trxList[i]
#         next1Trx = trxList[i + 1]; next2Trx = trxList[i + 2]; next3Trx = trxList[i + 3]; next4Trx = trxList[i + 4]
#         trxExons = list(map(lambda x: [x.chr(), x.start(), x.end(), x.strand()], sorted(trx.exons(), key=exonOrderKey_coords)))
#         print('trx{} - {}: start - {}, end, - {}, exon count - {}'.format(i+1, trx.trxID(), trx.start(), trx.end(), trx.exonCount()))
#         print('exons - {}'.format(trxExons))
#         print('skip: {}'.format(trx.skip()))
#         print('before gtfCount: {}'.format(trx.gtfCount()))
#         if i != len(trxList)-5:
#             if trx.skip() == 1: continue
#             if trx.isidentical(next1Trx):
#                 if count == 'O': trx.add_gtfCount(); print('gtfCount +1')
#                 next1Trx.add_skip()
#                 if trx.isidentical(next2Trx):
#                     if count == 'O': trx.add_gtfCount(); print('gtfCount +1')
#                     next2Trx.add_skip()
#                     if trx.isidentical(next3Trx):
#                         if count == 'O': trx.add_gtfCount(); print('gtfCount +1')
#                         next3Trx.add_skip()
#                         if trx.isidentical(next4Trx):
#                             if count == 'O': trx.add_gtfCount(); print('gtfCount +1')
#                             next4Trx.add_skip()
#             print('after gtfCount: {}'.format(trx.gtfCount()))
#             nTrxList.append(trx)
#             print('added to nTrxList\n')
#         if i == len(trxList)-5:
#             if trx.skip() != 1:
#                 if trx.isidentical(next1Trx):
#                     if count == 'O': trx.add_gtfCount(); print('gtfCount +1')
#                     next1Trx.add_skip()
#                     if trx.isidentical(next2Trx):
#                         if count == 'O': trx.add_gtfCount(); print('gtfCount +1')
#                         next2Trx.add_skip()
#                         if trx.isidentical(next3Trx):
#                             if count == 'O': trx.add_gtfCount(); print('gtfCount +1')
#                             next3Trx.add_skip()
#                             if trx.isidentical(next4Trx):
#                                 if count == 'O': trx.add_gtfCount(); print('gtfCount +1')
#                                 next4Trx.add_skip()
#                 print('after gtfCount: {}'.format(trx.gtfCount()))
#                 nTrxList.append(trx)
#                 print('added to nTrxList\n')
#             if next1Trx.skip() != 1:
#                 trxExons = list(map(lambda x: [x.chr(), x.start(), x.end(), x.strand()],
#                                     sorted(next1Trx.exons(), key=exonOrderKey_coords)))
#                 print('trx{} - {}: start - {}, end, - {}, exon count - {}'.format(i + 2, next1Trx.trxID(), next1Trx.start(), next1Trx.end(), next1Trx.exonCount()))
#                 print('exons - {}'.format(trxExons))
#                 print('skip: {}'.format(trx.skip()))
#                 print('before gtfCount: {}'.format(trx.gtfCount()))
#                 if next1Trx.isidentical(next2Trx):
#                     if count == 'O': next1Trx.add_gtfCount(); print('gtfCount +1')
#                     next2Trx.add_skip()
#                     if next1Trx.isidentical(next3Trx):
#                         if count == 'O': next1Trx.add_gtfCount(); print('gtfCount +1')
#                         next3Trx.add_skip()
#                         if next1Trx.isidentical(next4Trx):
#                             if count == 'O': next1Trx.add_gtfCount(); print('gtfCount +1')
#                             next4Trx.add_skip()
#                 print('after gtfCount: {}'.format(trx.gtfCount()))
#                 nTrxList.append(next1Trx)
#                 print('added to nTrxList\n')
#             if next2Trx.skip() != 1:
#                 trxExons = list(map(lambda x: [x.chr(), x.start(), x.end(), x.strand()],
#                                     sorted(next2Trx.exons(), key=exonOrderKey_coords)))
#                 print('trx{} - {}: start - {}, end, - {}, exon count - {}'.format(i + 3, next2Trx.trxID(), next2Trx.start(), next2Trx.end(), next2Trx.exonCount()))
#                 print('exons - {}'.format(trxExons))
#                 print('skip: {}'.format(trx.skip()))
#                 print('before gtfCount: {}'.format(trx.gtfCount()))
#                 if next2Trx.isidentical(next3Trx):
#                     if count == 'O': next2Trx.add_gtfCount(); print('gtfCount +1')
#                     next3Trx.add_skip()
#                     if next2Trx.isidentical(next4Trx):
#                         if count == 'O': next2Trx.add_gtfCount(); print('gtfCount +1')
#                         next4Trx.add_skip()
#                 print('after gtfCount: {}'.format(trx.gtfCount()))
#                 nTrxList.append(next2Trx)
#                 print('added to nTrxList\n')
#             if next3Trx.skip() != 1:
#                 trxExons = list(map(lambda x: [x.chr(), x.start(), x.end(), x.strand()],
#                                     sorted(next3Trx.exons(), key=exonOrderKey_coords)))
#                 print('trx{} - {}: start - {}, end, - {}, exon count - {}'.format(i + 4, next3Trx.trxID(), next3Trx.start(), next3Trx.end(), next3Trx.exonCount()))
#                 print('exons - {}'.format(trxExons))
#                 print('skip: {}'.format(trx.skip()))
#                 print('before gtfCount: {}'.format(trx.gtfCount()))
#                 if next3Trx.isidentical(next4Trx):
#                     if count == 'O': next3Trx.add_gtfCount(); print('gtfCount +1')
#                     next4Trx.add_skip()
#                 print('after gtfCount: {}'.format(trx.gtfCount()))
#                 nTrxList.append(next3Trx)
#                 print('added to nTrxList\n')
#             if next4Trx.skip() != 1:
#                 trxExons = list(map(lambda x: [x.chr(), x.start(), x.end(), x.strand()],
#                                     sorted(next4Trx.exons(), key=exonOrderKey_coords)))
#                 print('last transcript')
#                 print('trx{} - {}: start - {}, end, - {}, exon count - {}'.format(i + 5, next4Trx.trxID(), next4Trx.start(), next4Trx.end(), next4Trx.exonCount()))
#                 print('exons - {}'.format(trxExons))
#                 print('skip: {}'.format(trx.skip()))
#                 print('gtfCount: {}'.format(trx.gtfCount()))
#                 nTrxList.append(next4Trx)
#                 print('added to nTrxList\n')
#     return nTrxList

def group_overlapTrxs(trxList, genePrefix='.', addTssCps='O'):
    '''Groups overlapping input transcripts as new genes.
    If new name isn't given, 'Gene-ChrX_', 'Trx-ChrX_' prefixes are used.'''
    trxList = sorted(trxList, key=trxOrderKey_coords)
    lociDic = {}; lociList = []; nGeneList = []
    if genePrefix == '.': genePrefix = 'Gene'
    k = 2
    for i in range(len(trxList)):
        nTrx = trxList[i]
        chr = nTrx.chr(); strand = nTrx.strand()
        if i == 0:
            newLoci = '{}-Chr{}_1'.format(genePrefix, chr.lstrip('cChHrR'))
            lociDic[newLoci] = []
            lociDic[newLoci].append(nTrx)
            lociList.append(newLoci)
        else:
            #check last 10 loci
            if len(lociList) < 10:
                for prevLoci in lociList:
                    prevTrxs = lociDic[prevLoci]
                    overlap = False
                    for prevTrx in prevTrxs:
                        if nTrx.strandedOverlaps(prevTrx):
                            lociDic[prevLoci].append(nTrx)
                            overlap = True
                            break
                    if overlap: break
                else:
                    newLoci = '{}-Chr{}_{}'.format(genePrefix, chr.lstrip('cChHrR'), k)
                    if newLoci in lociDic:
                        raise Exception('new loci is in code: check for error')
                    lociDic[newLoci] = []
                    k += 1
                    lociDic[newLoci].append(nTrx)
                    lociList.append(newLoci)
            else:
                for prevLoci in lociList[len(lociList)-10:]:
                    prevTrxs = lociDic[prevLoci]
                    overlap = False
                    for prevTrx in prevTrxs:
                        if nTrx.strandedOverlaps(prevTrx):
                            lociDic[prevLoci].append(nTrx)
                            overlap = True
                            break
                    if overlap: break
                else:
                    newLoci = '{}-Chr{}_{}'.format(genePrefix, chr.lstrip('cChHrR'), k)
                    if newLoci in lociDic:
                        raise Exception('new loci is in code: check for error')
                    lociDic[newLoci] = []
                    k += 1
                    lociDic[newLoci].append(nTrx)
                    lociList.append(newLoci)
    def lociOrderKey(loci):
        numPattern = '\_\d+'
        locis = re.findall(numPattern, loci)
        if len(loci) == 0:
            raise Exception('weird loci name: {}'.format(loci))
        else:
            lociOrder = int(locis[0].strip('_'))
        return lociOrder
    k = 0
    for loci in sorted(lociDic.keys(), key=lociOrderKey):
        k += 1
        trxList = lociDic[loci]
        nGene_ID = loci
        nChr = trxList[0].chr(); nStrand = trxList[0].strand(); nGene_type = trxList[0].geneType(); nGene_name = trxList[0].geneName()
        nAttribute = trxList[0].attribute(); nSource = trxList[0].source(); nScore = trxList[0].score(); nFrame = trxList[0].frame()
        nGene = gene(gene_id=nGene_ID, chr=nChr, strand=nStrand, gene_type=nGene_type, gene_name=nGene_name,
                     attribute=nAttribute, source=nSource, score=nScore, frame=nFrame, trxs=trxList)
        if nGene_ID.startswith('Gene-'):
            nTrxs = nGene.trxs()
            for nTrx in nTrxs:
                nTrx_ID = nTrx.trxID().replace('Gene-', 'Trx-')
                nTrx.set_trxID(nTrx_ID)
            nGene.update(trxs=nTrxs)
        # print('gene{} - {}: start - {}, end - {}, trx count - {}'.format(k, nGene.geneID(), nGene.start(), nGene.end(), nGene.trxCount()))
        # trxs = list(map(lambda x: [x.chr(), x.start(), x.end(), x.strand(), x.exonCount(), x.exonLenTotal()], sorted(nGene.trxs(), key=trxOrderKey_coords)))
        # print('trxs - {}'.format(trxs))
        if addTssCps == 'O':
            nTrxs = nGene.trxs()
            for nTrx in nTrxs:
                if nTrx.tss() != '':
                    nTrxID = nTrx.trxID() + '_' + nTrx.tss()
                    if nTrx.cps() != '':
                        nTrxID += '_' + nTrx.cps()
                    nTrx.set_trxID(nTrxID)
                elif nTrx.cps() != '':
                    nTrxID = nTrx.trxID() + '_' + nTrx.cps()
                    nTrx.set_trxID(nTrxID)
            nGene.update(trxs=nTrxs)
        nGeneList.append(nGene)
    return nGeneList


def rem_geneLines(inGTF):
    outGTF = inGTF.replace('.gtf', '.no_genes.gtf')
    iF = open(inGTF, 'r')
    gtfLines = iF.readlines(); iF.close()
    headerLines = []; geneLines = []; nonGeneLines = []
    for line in gtfLines:
        if line.startswith('#') or line.count('_') != 8:
            headerLines.append(line+'\n')
            continue
        elif line.split('\t')[2] == 'gene':
            geneLines.append(line+'\n')
        else:
            nonGeneLines.append(line+'\n')
    with open(outGTF, 'w') as oF:
        for header in headerLines:
            oF.write(header)
        for line in nonGeneLines:
            oF.write(line)
    return headerLines + nonGeneLines

def rem_GTFattribute(inGTF, attributeName, feature='.'):
    '''Removes attribute from GTF file, can limit to specified feature (e.g. transcript)'''
    outGTF = inGTF.replace('.gtf', '.{}X.gtf'.format(attributeName))
    iF = open(inGTF, 'r')
    gtfLines = iF.readlines(); iF.close()
    with open(outGTF, 'w') as oF:
        for line in gtfLines:
            if line.count('\t') != 8:
                oF.write(line)
                print('weird line: {}'.format(line)); continue
            lineList = line.strip('\n').split('\t')
            Attributes = lineList[8].rstrip(';')
            if Attributes.find(attributeName) < 0: continue
            lineFeature = lineList[2]
            if feature != '.':
                if lineFeature != feature: continue
            AttributeList = Attributes.split(';')
            newAttributeList = [x.strip() for x in AttributeList if x.split('"')[0].strip() != attributeName]
            lineList[8] = '; '.join(newAttributeList) + ';'
            newLine = '\t'.join(lineList) + '\n'
            oF.write(newLine)

# def write_combinedGTF(outGTF, *inGTFs):
#     '''Change so input can be genes OR transcripts!!!!!!!!!@@'''
#     allLines = []
#     k = 0
#     for gtf in inGTFs:
#         k += 1
#         iF = open(gtf, 'r')
#         gtfLines = iF.readlines(); iF.close()
#         trxLines = filter(lambda x: x.split('\t')[2] == 'transcript', gtfLines)
#         print('gtf{}: {}'.format(k, gtf))
#         print('gtf{} lines: {}'.format(k, len(gtfLines)))
#         print('gtf{} transcripts: {}\n'.format(k, len(trxLines)))
#         allLines.extend(gtfLines)
#     outGtfLines = list(set(allLines))
#     outGtfLines.sort(key=gtfOrderKey_gene)
#     outTrxLines = list(filter(lambda x: x.split('\t')[2] == 'transcript', outGtfLines))
#     with open(outGTF, 'w') as oF:
#         for line in outGtfLines:
#             oF.write(line)
#     # stats
#     print('combined GTF: {}'.format(outGTF))
#     print('combined GTF lines: {}'.format(len(outGtfLines)))
#     print('combined GTF transcripts: {}'.format(len(outTrxLines)))
#     return outGtfLines

def write_GTF(chrDic, outGTF, feature, mode='.'):
    '''a'''
    start = time.time()
    print('writing: {}'.format(outGTF))
    outLines = []
    if mode == '.': mode = 'w'
    with open(outGTF, mode) as oF:
        chrList = sorted(chrDic.keys(), key=chrOrderKey)
        for chr in chrList:
            if feature == 'gene':
                geneList = sorted(chrDic[chr], key=geneOrderKey_coords)
                for gene in geneList:
                    source = gene.source(); geneFeature = 'gene'
                    geneStart = gene.start(); geneEnd = gene.end(); geneScore = gene.score(); geneStrand = gene.strand(); geneFrame = gene.frame()
                    geneGene_id = gene.geneID(); geneGene_type = gene.geneType(); geneGene_name = gene.geneName(); geneAttribute = gene.attribute()
                    geneLine = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(chr, source, geneFeature, geneStart, geneEnd, geneScore, geneStrand, geneFrame)\
                              + 'gene_id "{}";'.format(geneGene_id)
                    if geneGene_type != '': geneLine += ' gene_type "{}";'.format(geneGene_type)
                    if geneGene_name != '': geneLine += ' gene_name "{}";'.format(geneGene_name)
                    if geneAttribute != '': geneLine += ' attribute "{}";'.format(geneAttribute)
                    outLines.append(geneLine + '\n')
                    oF.write(geneLine + '\n')
                    trxList = sorted(gene.trxs(), key=trxOrderKey_coords)
                    for trx in trxList:
                        trxFeature = 'transcript'
                        trxStart = trx.start(); trxEnd = trx.end(); trxScore = '.'; trxStrand = trx.strand(); trxFrame = '.'
                        trxGene_id = trx.geneID(); trxTrx_id = trx.trxID(); trxGene_type = trx.geneType(); trxGene_name = trx.geneName()
                        trxAttribute = trx.attribute(); cpc_score = trx.cpcScore(); cpat_score = trx.cpatScore(); gtf_count = trx.gtfCount()
                        tssPattern = 'TSS\.\d+'
                        tssMatches = re.findall(tssPattern, trxTrx_id)
                        if len(tssMatches) >= 1:
                            trxTSS = tssMatches[0]
                        else:
                            trxTSS = trx.tss()
                        cpsPattern = 'CPS\.\d+'
                        cpsMatches = re.findall(cpsPattern, trxTrx_id)
                        if len(cpsMatches) >= 1:
                            trxCPS = cpsMatches[0]
                        else:
                            trxCPS = trx.cps()
                        trxLine = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(chr, source, trxFeature, trxStart, trxEnd, trxScore, trxStrand, trxFrame) \
                                  + 'gene_id "{}"; transcript_id "{}";'.format(trxGene_id, trxTrx_id)
                        if trxTSS != '': trxLine += ' TSS "{}";'.format(trxTSS)
                        if trxCPS != '': trxLine += ' CPS "{}";'.format(trxCPS)
                        if trxGene_type != '': trxLine += ' gene_type "{}";'.format(trxGene_type)
                        if trxGene_name != '': trxLine += ' gene_name "{}";'.format(trxGene_name)
                        if trxAttribute != '': trxLine += ' attribute "{}";'.format(trxAttribute)
                        if cpc_score != '': trxLine += ' cpc_score "{}";'.format(cpc_score)
                        if cpat_score != '': trxLine += ' cpat_score "{}";'.format(cpat_score)
                        if gtf_count != '': trxLine += ' gtf_count "{}";'.format(gtf_count)
                        outLines.append(trxLine + '\n')
                        oF.write(trxLine + '\n')
                        exons = sorted(trx.exons(), key=exonOrderKey_coords)
                        exNum = 0
                        for exon in exons:
                            exFeature = 'exon'
                            exStart = exon.start(); exEnd = exon.end(); exScore = '.'; exStrand = exon.strand(); exFrame = '.'
                            exGene_id = exon.geneID(); exTrx_id = exon.trxID(); exGene_type = exon.geneType(); exGene_name = exon.geneName()
                            exAttribute = exon.attribute()
                            if exGene_type == '': exGene_type = trxGene_type
                            if exGene_name == '': exGene_name = trxGene_name
                            if exAttribute == '': exAttribute = trxAttribute
                            # exNum = exon.exonNum()
                            exNum += 1
                            exLine = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(chr, source, exFeature, exStart, exEnd, exScore, exStrand, exFrame) \
                                     + 'gene_id "{}"; transcript_id "{}"; exon_number "{}";'.format(exGene_id, exTrx_id, exNum)
                            if exGene_type != '': exLine += ' gene_type "{}";'.format(exGene_type)
                            if exGene_name != '': exLine += ' gene_name "{}";'.format(exGene_name)
                            if exAttribute != '': exLine += ' attribute "{}";'.format(exAttribute)
                            oF.write(exLine + '\n')
                            outLines.append(exLine + '\n')
            elif feature == 'transcript':
                trxList = sorted(chrDic[chr], key=trxOrderKey_coords)
                for trx in trxList:
                    source = trx.source(); trxFeature = 'transcript'
                    trxStart = trx.start(); trxEnd = trx.end(); trxScore = trx.score(); trxStrand = trx.strand(); trxFrame = trx.frame()
                    trxGene_id = trx.geneID(); trxTrx_id = trx.trxID(); trxGene_type = trx.geneType(); trxGene_name = trx.geneName()
                    trxAttribute = trx.attribute(); cpc_score = trx.cpcScore(); cpat_score = trx.cpatScore(); gtf_count = trx.gtfCount()
                    tssPattern = 'TSS\.\d+'
                    tssMatches = re.findall(tssPattern, trxTrx_id)
                    if len(tssMatches) >= 1:
                        trxTSS = tssMatches[0]
                    else:
                        trxTSS = trx.tss()
                    cpsPattern = 'CPS\.\d+'
                    cpsMatches = re.findall(cpsPattern, trxTrx_id)
                    if len(cpsMatches) >= 1:
                        trxCPS = cpsMatches[0]
                    else:
                        trxCPS = trx.cps()
                    trxLine = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(chr, source, trxFeature, trxStart, trxEnd, trxScore, trxStrand, trxFrame)\
                              +'gene_id "{}"; transcript_id "{}";'.format(trxGene_id, trxTrx_id)
                    if trxTSS != '': trxLine += ' TSS "{}";'.format(trxTSS)
                    if trxCPS != '': trxLine += ' CPS "{}";'.format(trxTSS)
                    if trxGene_type != '': trxLine += ' gene_type "{}";'.format(trxGene_type)
                    if trxGene_name != '': trxLine += ' gene_name "{}";'.format(trxGene_name)
                    if trxAttribute != '': trxLine += ' attribute "{}";'.format(trxAttribute)
                    if cpc_score != '': trxLine += ' cpc_score "{}";'.format(cpc_score)
                    if cpat_score != '': trxLine += ' cpat_score "{}";'.format(cpat_score)
                    if gtf_count != '': trxLine += ' gtf_count "{}";'.format(gtf_count)
                    outLines.append(trxLine + '\n')
                    oF.write(trxLine + '\n')
                    exons = sorted(trx.exons(), key=exonOrderKey_coords)
                    exNum = 0
                    for exon in exons:
                        exFeature = 'exon'
                        exStart = exon.start(); exEnd = exon.end(); exScore = '.'; exStrand = exon.strand(); exFrame = '.'
                        exGene_id = exon.geneID(); exTrx_id = exon.trxID(); exGene_type = exon.geneType(); exGene_name = exon.geneName()
                        exAttribute = exon.attribute()
                        if exGene_type == '': exGene_type = trxGene_type
                        if exGene_name == '': exGene_name = trxGene_name
                        if exAttribute == '': exAttribute = trxAttribute
                        # exNum = exon.exonNum()
                        exNum += 1
                        exLine = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(chr, source, exFeature, exStart, exEnd, exScore, exStrand, exFrame)\
                                  + 'gene_id "{}"; transcript_id "{}"; exon_number "{}";'.format(exGene_id, exTrx_id, exNum)
                        if exGene_type != '': exLine += ' gene_type "{}";'.format(exGene_type)
                        if exGene_name != '': exLine += ' gene_name "{}";'.format(exGene_name)
                        if exAttribute != '': exLine += ' attribute "{}";'.format(exAttribute)
                        oF.write(exLine + '\n')
                        outLines.append(exLine + '\n')
    print('done: {} sec.'.format(time.time() - start))
    return outLines

def write_trxFilteredGTF(inGTF, trxList, outGTF):
    if outGTF.strip('/').count('/') == 0:
        outGTF = os.getcwd() + '/' + outGTF
        print('no directory specified for output GTF: writing as {}'.format(outGTF))
    trxList = list(map(lambda x: x.lower(), trxList))
    iF = open(inGTF, 'r')
    inGtfLines = iF.readlines(); iF.close()
    outGtfLines = list(set(filter(lambda x: next(j for j in x.split('\t')[8].split(';') if 'transcript_id' in j).split('"')[1].lower() in trxList, inGtfLines)))
    outGtfLines.sort(key=gtfOrderKey_gene)
    if len(outGtfLines) == 0:
        print('no shared transcripts in input GTF: aborting'); sys.exit()
    outGtfTrxLines = filter(lambda x: x.split('\t')[2] == 'transcript', outGtfLines)
    trxList = list(map(lambda x: next(j for j in x.split('\t')[8].split(';') if 'transcript_id' in j).split('"')[1], outGtfTrxLines))
    uniqTrxList = list(set(trxList))
    d = defaultdict(int)
    for trx in trxList:
        d[trx] += 1
    items = d.items()
    twoUpList = [x for x in items if x[1] >= 2]
    with open(outGTF, 'w') as outGTF:
        for line in outGtfLines:
            outGTF.write(line)
    #stats
    print('number of total transcripts in output GTF: {}'.format(len(outGtfTrxLines)))
    print('number of unique transcripts in output GTF: {}'.format(len(uniqTrxList)))
    print('number of transcripts with 2 or more copies: {}'.format(len(twoUpList)))
    print('transcripts with 2 or more copies: {}'.format(twoUpList))
    return outGtfLines
