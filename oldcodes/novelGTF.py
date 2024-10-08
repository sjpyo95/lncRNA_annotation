##### novelGTF.py
##### ""
##### Made by Sangjin Lee
##### Last modified: 2020.09.05
import os, sys
import SJ_module as sj
import re

mode = sys.argv[1] # novel/update or 1/2
if mode == '1' or mode == 'novel':
    gene_prefix = sys.argv[2]
    outGTF = sys.argv[3]
    inGTFs = sys.argv[4:]
    if outGTF == '.':
        outGTF = '/home/sjlee/data/Bcell/RPDs/Meta-Assembly/novelGTF/novelGenes.FILL.gtf'
        # /home/sjlee/data/Bcell/RPDs/Meta-Assembly/novelGTF/novelLncGenes.b1b2.gtf
    if '.noncoding.' in outGTF:
        print("Can't have 'noncoding' in output GTF"); sys.exit()
    if len(inGTFs) > 5:
        print("Max. 5 GTFs at a time"); sys.exit()
    if os.path.exists(outGTF):
        print('{} already exists: aborting'.format(outGTF)); sys.exit()

    allDic = {}; ncDic = {}; geneCount = 0; ncGeneCount = 0
    for i in range(len(inGTFs)):
        gtf = inGTFs[i]
        if not gtf.startswith('/home/'):
            idx = gtf.find('/home/')
            gtf = gtf[idx:]
            inGTFs[i] = gtf
        gtfDic = sj.get_GTFtrxs(gtf)
        for chr in gtfDic.keys():
            if chr not in allDic: allDic[chr] =[]
            gtfTrxs = gtfDic[chr]
            allDic[chr].extend(gtfTrxs)
    for chr in sorted(allDic.keys(), key=sj.chrOrderKey):
        # if chr.lower() != 'chr1': continue
        trxs = sj.get_uniqTrxs(allDic[chr])
        genes = sj.group_overlapTrxs(trxs, gene_prefix, addTssCps='X')
        allDic[chr] = genes; geneCount += len(genes)
        ncGenes = sj.filter_ncGenes(genes, cpc_cutoff=-0.3, cpat_cutoff=0.44)
        ncDic[chr] = ncGenes; ncGeneCount += len(ncGenes)

    # sys.exit()
    #
    # if os.path.exists(outGTF):
    #     print('{} already exists: aborting'.format(outGTF)); sys.exit()

    # write comprehensive GTF (coding + noncoding genes)
    with open(outGTF, 'w') as oF:
        oF.write('#number of GTFs combined: {}\n'.format(len(inGTFs)))
        num = 0
        for gtf in inGTFs:
            num += 1
            oF.write('#GTF{}: {}\n'.format(num, gtf))
    sj.write_GTF(allDic, outGTF, 'gene', 'a')

    # write noncoding GTF (only noncoding genes)
    with open(outGTF.replace('.gtf', '.noncoding.gtf'), 'w') as oF:
        oF.write('#number of GTFs combined: {}\n'.format(len(inGTFs)))
        num = 0
        for gtf in inGTFs:
            num += 1
            oF.write('#GTF{}: {}\n'.format(num, gtf))
    sj.write_GTF(ncDic, outGTF.replace('.gtf', '.noncoding.gtf'), 'gene', 'a')

    # write stats
    with open(outGTF.replace('.gtf', '.stats'), 'w') as oF:
        oF.write('#number of GTFs combined: {}\n'.format(len(inGTFs)))
        num = 0
        for gtf in inGTFs:
            num += 1
            oF.write('#GTF{}: {}\n'.format(num, gtf))
        oF.write('\n')
        headerLine = '#novel_GTF\tgenes\tnon-coding_genes\tcoding_genes\n'
        oF.write(headerLine)
        statLine = '{}\t{}\t{}\t{}\n'.format(os.path.basename(outGTF), geneCount, ncGeneCount, geneCount-ncGeneCount)
        oF.write(statLine)


elif mode == '2' or mode == 'update':
    gene_prefix = sys.argv[2]
    gtf2update = sys.argv[3]
    inGTFs = sys.argv[4:]
#/home/sjlee/data/Bcell/RPDs/Meta-Assembly/all/lncRNA_anno/allG_overlap_filtered/no_overlap/transcripts.no_overlap.combined.cp.gtf
#/home/sjlee/data/Bcell/RPDs/Meta-Assembly/b1/lncRNA_anno/allG_overlap_filtered/no_overlap/transcripts.no_overlap.combined.cp.gtf
#/home/sjlee/data/Bcell/RPDs/Meta-Assembly/b2/lncRNA_anno/allG_overlap_filtered/no_overlap/transcripts.no_overlap.combined.cp.gtf

    versPattern = '\.v\d+.'
    gtfVersions = re.findall(versPattern, gtf2update)
    if len(gtfVersions) == 0:
        outGTF = gtf2update.replace('.gtf', '.v2.gtf')
    elif len(gtfVersions) == 1:
        version = int(gtfVersions[0].split('v')[1].split('.')[0])
        outGTF = gtf2update.split(gtfVersions[0])[0] + '.v{}.gtf'.format(version+1)
    else:
        raise Exception('check GTF version')
    if os.path.exists(outGTF):
        print('{} already exists: aborting'.format(outGTF)); sys.exit()
    if len(inGTFs) > 4:
        print("Max. 5 GTFs at a time (including GTF to update)"); sys.exit()
    if '.noncoding.' in gtf2update:
        print("Can't update 'noncoding' gtf: update comprehensive gtf"); sys.exit()
    headerLines = sj.get_headerLines(gtf2update)
    prevGTFlines = [x for x in headerLines if x.startswith('#GTF') and x.endswith('.gtf\n')]
    prevGTFs = list(map(lambda x:x.split(':')[1].strip(' \n'), prevGTFlines))
    num = int(prevGTFlines[-1].split(':')[0].strip('#GTF'))
    newHeaderLines = prevGTFlines

    allDic = {}; ncDic = {}; geneCount = 0; ncGeneCount = 0
    gtfDic = sj.get_GTFtrxs(gtf2update)
    for chr in gtfDic.keys():
        if chr not in allDic: allDic[chr] = []
        gtfTrxs = gtfDic[chr]
        allDic[chr].extend(gtfTrxs)
    noNew = 0
    for gtf in inGTFs:
        if not gtf.startswith('/home/'):
            idx = gtf.find('/home/')
            gtf = gtf[idx:]
        if gtf in prevGTFs:
            print('{} already in novel GTF'.format(gtf))
            noNew += 1; continue
        num += 1
        newHeaderLines.append('#GTF{}: {}\n'.format(num, gtf))
        gtfDic = sj.get_GTFtrxs(gtf)
        for chr in gtfDic.keys():
            if chr not in allDic: allDic[chr] =[]
            gtfTrxs = gtfDic[chr]
            allDic[chr].extend(gtfTrxs)
    if noNew == len(inGTFs):
        print('all given GTFs are already in novel GTF'); sys.exit()
    newHeaderLines.insert(0, '#number of GTFs combined: {}\n'.format(num))
    for chr in sorted(allDic.keys(), key=sj.chrOrderKey):
        # if chr.lower() != 'chr1': continue
        trxs = sj.get_uniqTrxs(allDic[chr])
        genes = sj.group_overlapTrxs(trxs, gene_prefix, addTssCps='X')
        allDic[chr] = genes; geneCount += len(genes)
        ncGenes = sj.filter_ncGenes(genes, cpc_cutoff=-0.3, cpat_cutoff=0.44)
        ncDic[chr] = ncGenes; ncGeneCount += len(ncGenes)

    # sys.exit()
    #
    # if os.path.exists(outGTF):
    #     print('{} already exists: aborting'.format(outGTF)); sys.exit()

    # write comprehensive GTF (coding + noncoding genes)
    with open(outGTF, 'w') as oF:
        for line in newHeaderLines:
            oF.write(line)
    sj.write_GTF(allDic, outGTF, 'gene', 'a')

    # write comprehensive GTF (coding + noncoding genes)
    with open(outGTF.replace('.gtf', '.noncoding.gtf'), 'w') as oF:
        for line in newHeaderLines:
            oF.write(line)
    sj.write_GTF(ncDic, outGTF.replace('.gtf', '.noncoding.gtf'), 'gene', 'a')

    # write stats
    with open(outGTF.replace('.gtf', '.stats'), 'w') as oF:
        for line in newHeaderLines:
            oF.write(line)
        oF.write('\n')
        headerLine = '#novel_GTF\tgenes\tnon-coding_genes\tcoding_genes\n'
        oF.write(headerLine)
        statLine = '{}\t{}\t{}\t{}\n'.format(os.path.basename(outGTF), geneCount, ncGeneCount, geneCount-ncGeneCount)
        oF.write(statLine)
