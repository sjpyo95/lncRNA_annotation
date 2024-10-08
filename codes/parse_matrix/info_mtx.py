import sys, os
import module as mdl

ig = ['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene', 'IG_V_gene', 'TR_C_gene', 'TR_J_gene', 'TR_V_gene', 'TR_D_gene', 'IG_pseudogene', 'IG_C_pseudogene', 'IG_J_pseudogene', 'IG_V_pseudogene', 'TR_V_pseudogene', 'TR_J_pseudogene']

ncRNA = ['Mt_rRNA', 'Mt_tRNA', 'miRNA', 'misc_RNA', 'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'ribozyme', 'sRNA', 'scaRNA', 'Mt_tRNA_pseudogene', 'tRNA_pseudogene', 'snoRNA_pseudogene', 'snRNA_pseudogene', 'scRNA_pseudogene', 'rRNA_pseudogene', 'misc_RNA_pseudogene', 'miRNA_pseudogene', 'vaultRNA']

def getGenetype(rnas):
	genetypedic = dict()
	for chr in rnas.keys():
		genes = rnas[chr]
		for gene in genes:
			id = gene.geneid(); genetype = gene.genetype()
			genetypedic[id] = genetype
	return genetypedic

mtx_file = sys.argv[1]
ref_file = sys.argv[2]
known_file = sys.argv[3]
novel_file = sys.argv[4]
outfile = sys.argv[5]

mtx, samples = mdl.getMatrix(mtx_file)

ref = getGenetype(mdl.getGtf2(ref_file))
knownlncRNAs = getGenetype(mdl.getGtf2(known_file))
novellncRNAs = getGenetype(mdl.getGtf2(novel_file))
lncs = knownlncRNAs.keys() + novellncRNAs.keys()
elseGlist = [x for x in ref.keys() if x not in lncs]
elseGenes = dict()
for gene in ref.keys():
	if gene in elseGlist:
		elseGenes[gene] = ref[gene]

newMtx = dict()
pcgN = 0; igN = 0; ncN = 0; kl = 0; nl = 0; pl = 0
for idsymb in mtx.keys():
	id = idsymb[0]; symb = idsymb[1]
	expdic = mtx[idsymb]
	if id in knownlncRNAs.keys():
		type = 'Known_lncRNA'; genetype = knownlncRNAs[id]
		kl += 1
	elif id in novellncRNAs.keys():
		type = 'Novel_lncRNA'; genetype = novellncRNAs[id]
		nl += 1
	elif id in elseGenes.keys():
		genetype = elseGenes[id]
		if genetype == 'protein_coding':
			type = 'PCG'
			pcgN += 1
		elif genetype in ig:
			type = 'Ig'
			igN += 1
		elif genetype in ncRNA:
			type = 'ncRNA'
			ncN += 1
		elif genetype == 'lncRNA':
			type = 'Known_lncRNA'
			kl+=1
		elif genetype == 'lincRNA':
			type = 'Known_lncRNA'
			kl += 1
		elif genetype == 'sense_intronic':
			type = 'Known_lncRNA'
			genetype = 'intervening'
			kl += 1
		else:
			type = 'pseudo'
			pl += 1
	else:
		print idsymb
#		type = genetype
		continue
	
	newMtx[(id,symb,type,genetype)] = expdic

print 'PCG:', pcgN
print 'Ig:', igN
print 'ncRNA:', ncN
print 'Known lncRNA:', kl
print 'Novel lncRNA:', nl
print 'Pseudogene:', pl

mdl.writeinfoMatrix(newMtx, samples, outfile)


