#!/usr/bin/env python
import sys, os
import glob

def readREF(iFP):
    idD = {}
    for line in open(iFP):
        gene_id = line.split('gene_id "')[-1].split(';')[0].rstrip().replace('"','')
        gene_name = line.split('gene_name "')[-1].split(';')[0].rstrip().replace('"','')
        if gene_id in idD and gene_name != idD[gene_id]:
            print 'ERROR, ', gene_id, gene_name
            exit()
        else:
            idD[gene_id] = gene_name
    return idD

def getAssigned(iFP):
    if not os.path.exists(iFP):
        print ('summary file not found')
        print (iFP)
        exit()
    for line in open(iFP):
        if line.startswith('Assigned'):
            nAssigned = int(line.strip().split()[-1])
            return nAssigned
    else:
        print('ERROR 1: Assigned line not found')
        exit()

def readSAF(iFP, nTotalRead):
    expD = {}
    for line in open(iFP):
        if line.startswith('#'): continue
        if line.startswith('Geneid'): continue
        items = line.rstrip('\n').split('\t')
        geneID = items[0]
        length = int(items[-2])
        reads = int(items[-1])
        FPKM = 1000000000.0 * reads / length / nTotalRead
        if geneID in expD:
            print('ERROR 2: key exists\n', geneID, expD)
            exit()
        expD[geneID] = FPKM
    return expD

def printTXT(expDD, oFP, idD):
    oF = open(oFP, 'w')
    oF.write( '\t'.join(['NAME', 'DESCRIPTION'] + sorted(expDD.keys())) + '\n')
    geneList = expDD.values()[0].keys()
    for gene in geneList:
        gene_name = idD[gene]
        oF.write ( '\t'.join([gene_name, '-'] + \
                [str(expDD[sample][gene]) for sample in sorted(expDD.keys())]) + '\n')
    oF.close()

def printCLS(expDD, oFP):
    oF = open(oFP, 'w')
    totalSample = str(len(expDD.keys()))
    totalCase = str(2) # tumor and normal
    oF.write( totalSample + ' ' + totalCase + ' 1\n')
    oF.write( '# '+group1+' '+group2+' \n' ) 
    oF.write( ' '.join(['0' if sample.startswith(group1) else '1' for sample in sorted(expDD.keys())]))
    oF.close()

if __name__=='__main__':
#    if len(sys.argv) != 3:
#        print('usage: python ~.py input_folder_path')
#        print('')
#        print('this code will automatically find all files ending with ".saf" extender')
#        print('from given directory. and convert to GCT format with FPKM value')
#        exit()
	iFolP = sys.argv[1]
	group = sys.argv[2]	#group1.group2
	group1 = group.split('.')[0]; group2 = group.split('.')[1]
	iFPL = glob.glob(iFolP + '/*.txt')
	oFP_txt = iFolP + '/expressions.txt'
	oFP_cls = iFolP + '/samples.cls'

	expDD = {}
	for iFP in sorted(iFPL):
		sample = iFP.split('/')[-1]
		iFP_summary = iFP + '.summary'
		nTotalRead = getAssigned(iFP_summary)
		expD = readSAF(iFP, nTotalRead)
		expDD[sample] = expD
	ref = '/home/sjpyo/new/data/anno/mouse/gencode.vM24.annotation.gtf'
	idD = readREF(ref)
	printTXT(expDD, oFP_txt, idD)
	printCLS(expDD, oFP_cls)





