import sys, os
import module as mdl

matrix_file = sys.argv[1]
ref_file = sys.argv[2]
novel_file = sys.argv[3]
out_file = sys.argv[4]

mtx, samples = mdl.getMatrix(matrix_file)
ref = mdl.getGtf2(ref_file)
novel = mdl.getGtf2(novel_file)

refgenes = sum(ref.values(), [])
novelgenes = sum(novel.values(), [])
genes = refgenes + novelgenes

pcgs = ['protein_coding','IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene', 'IG_V_gene', 'TR_C_gene', 'TR_J_gene', 'TR_V_gene', 'TR_D_gene']
known_lncRNAs = ['lncRNA','3prime_overlapping_ncRNA', 'antisense', 'bidirectional_promoter_lncRNA', 'lincRNA, macro_lncRNA', 'non_coding', 'processed_transcript', 'sense_intronic' ,'sense_overlapping']

new_mtx = dict()
for idsymb in mtx.keys():
	id = idsymb[0]; symb = idsymb[1]
	exps = mtx[idsymb]
	gene = filter(lambda x: x.geneid() == id, genes)[0]
	genetype = gene.genetype()
	if gene in refgenes:
		if genetype in pcgs:
			typ = 'PCG'
		elif genetype in known_lncRNAs:
			typ = 'Known_lncRNA'
		else:
#			print genetype
			continue
	elif gene in novelgenes:
		typ = 'Novel_lncRNA'
	else:
#		print genetype
		continue
	new_mtx[(id,symb,typ,genetype)] = exps

outfile = open(out_file, 'w')
outfile.write('ID\tsymbol\ttype\tgenetype\t'+'\t'.join(samples))
lines = ''
for gene in new_mtx.keys():
	id = gene[0]; symbol = gene[1]; type = gene[2]; genetype = gene[3]
	line = id+'\t'+symbol+'\t'+type+'\t'+genetype
	express = new_mtx[gene]
	for sample in samples:
		line += '\t'+ str(express[sample])
	lines += '\n'+line
outfile.write(lines)
outfile.close()

