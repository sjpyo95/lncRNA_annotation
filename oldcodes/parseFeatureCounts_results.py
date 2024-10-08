#!/usr/bin/env python
import sys, os
import glob

def getAssigned(iFP):
    for line in open(iFP):
        if line.startswith('Assigned'):
            nAssigned = int(line.strip().split()[-1])
            return nAssigned
    else:
        print 'ERROR 1: Assigned line not found'

def getExpDD(iFPL):
    expDD = dict()
    for iFP in iFPL:
        nTotalRead = getAssigned(iFP + '.summary')
        iFP_ = os.path.basename(iFP)
        expDD[iFP_] = dict()
        for line in open(iFP):
            if line.startswith('#'): continue
            if line.startswith('Geneid'): continue
            itemL = line.strip().split()
            gid = itemL[0]
            nRead = int(itemL[-1])
            expDD[iFP_][gid] = nRead
    return expDD

def printReadCountTable(expDD, iFPL, oFP):
    oF = open(oFP, 'w')
    ## aquire list of gene
    iFP_random = expDD.keys()[0]
    gidL = expDD[iFP_random].keys()
    gidL.sort(key=lambda x: int(x.lstrip('ENSG')[-1]))

    ## print tag
    oF.write('genename\t' + '\t'.join(map(os.path.basename, iFPL)) + '\n')
    for gid in gidL:
        line = gid + '\t'
        for iFP in iFPL:
            iFP = os.path.basename(iFP)
            exp = expDD[iFP][gid]
            line += str(exp) + '\t'
        oF.write(line.strip() + '\n')
    oF.close()

def printLabelTable(iFPL, oFP):
    oF = open(oFP, 'w')
    oF.write('Sample\tgroup\n')
    for iFP in iFPL:
        iFP = os.path.basename(iFP)
        group = 'tumor' if iFP.startswith('tumor') else 'normal'
        oF.write(iFP.strip() + '\t' + group + '\n')
    oF.close()

if __name__=='__main__':
    if len(sys.argv) != 2:
        print 'usage: python ~.py input_folder'
        print ''
        print 'this code will find all input_folder/*.saf and create'
        print 'readCountTable.txt, labelTable.txt in given directory'
        exit()
    iFolP = sys.argv[1]
    oFP_rct = iFolP + 'readCountTable.txt'
    oFP_lt = iFolP + 'labelTable.txt'
    iFPL = sorted(glob.glob(iFolP + '/*.txt'))
    expDD = getExpDD(iFPL)
    printReadCountTable(expDD, iFPL, oFP_rct)
    printLabelTable(iFPL, oFP_lt)
