import sys
import module as mdl

infile = sys.argv[1]
source = sys.argv[2]
genefix = sys.argv[3]
attri = sys.argv[4]
outfile = open(sys.argv[5], 'w')
stat = sys.argv[6]
gtfD = mdl.getGtf(infile)
geneDic = dict()

def getGeneloc(trxs):
	minstart = min([int(x.start()) for x in trxs])
	maxend = max([int(x.end()) for x in trxs])
	return minstart, maxend

def gettsscps(trxid):
	sup = []
	if 'TSS' in trxid: sup.append('TSS')
	if 'CPS' in trxid: sup.append('CPS')
	else: sup = ['None']
	sup = ','.join(sup)
	return sup

def getGeneattri(trxs):
	end = ','.join(list(set(sum(map(lambda x: x.attri().split('|')[0].split(','), trxs), []))))
	celltypes = ','.join(list(set(sum(map(lambda x: x.attri().split('|')[-1].split(','), trxs),[]))))
	geneattri = end+'|'+celltypes
	return geneattri
	

for chrom in gtfD.keys():
	trxs = gtfD[chrom]
	for trx in trxs:
		geneid = trx.geneid()
		if not geneDic.has_key(geneid):
			geneDic[geneid] = []
		geneDic[geneid].append(trx)

i = 0
for geneid in geneDic.keys():
#	newGeneid = genefix + '-' + geneid.split('.')[0]
	if 'novel' in stat:
		newGeneid = geneid
		i += 1
		j = '0'*(4-len(str(i)))
		gidn = j+str(i)
		newGeneid = genefix + '_' + gidn
		geneattri = attri
#		newGeneid = geneid.split('lncRNA-')[-1]
	else:
		newGeneid = geneid
	
	trxs = geneDic[geneid]
	geneStart, geneEnd = getGeneloc(trxs)
	genetype = list(set(map(lambda x: x.genetype(), trxs)))
#	if len(genetype) > 1:
#		genetype = '+'.join(genetype)
#	else:
	genetype = genetype[0]
#	attri = trxs[0].attri()
#	attri = geneid
#		print geneid
#	if attri == '.':
#		attri = trxs[0].attri()
	if attri == 'stay':
		geneattri = getGeneattri(trxs)
	geneline = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; gene_type \"%s\"; gene_name \"%s\"; attribute \"%s\";' % (trxs[0].chr(), source, 'gene', geneStart, geneEnd, '.', trxs[0].sense(), '.', newGeneid, genetype, trxs[0].symb(), geneattri)
	outfile.write(geneline + '\n')
	x = 0
	for trx in trxs:
		if attri == 'stay':
			trxattri = trx.attri()
		else:
			endsupstat = gettsscps(trx.trxid())
			trxattri = endsupstat + '|' + attri
#		trxid = trx.trxid()
		if 'novel' in stat:
			x += 1
#			y = '0'*(3-len(str(x)))
#			tidn = y+str(x)
			tidn = str(x)
			newTrxid = newGeneid + '.' + tidn
#			newTrxid = trxid
#		elif 'known' in trxid:
#			newTrxid = trxid.split('lncRNA-')[-1]
		else:
			newTrxid = trx.trxid()
		trxline = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; gene_type \"%s\"; gene_name \"%s\"; attribute \"%s\";' % (trx.chr(), source, 'transcript', trx.start(), trx.end(), '.', trx.sense(), '.', newGeneid, newTrxid, trx.genetype(), trx.symb(), trxattri)
		exons = trx.exons()
		outfile.write(trxline + '\n')
		for ex in exons:
			exonline = '%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; gene_type \"%s\"; gene_name \"%s\"; attribute \"%s\";' % (ex.chr(), source, 'exon', ex.start(), ex.end(), '.', trx.sense(), '.', newGeneid, newTrxid, trx.genetype(), trx.symb(), trxattri)
			outfile.write(exonline+'\n')
outfile.close()
