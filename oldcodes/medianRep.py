import sys
def getMedian(a):
	a_len = len(a)
	if (a_len == 0): return None
	a_center = int(a_len / 2)
	if (a_len % 2 == 1):
		return a[a_center]
	else:
		return (float(a[a_center - 1]) + float(a[a_center])) / 2.0

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')
#metaData = open(sys.argv[3], 'w')
lines = infile.readlines(); infile.close()
#samples = map(lambda x: x.split('_rep')[0], lines[0].strip().split('\t')[1:])
samples = map(lambda x: x, lines[0].strip().split('\t')[1:])
#sams = map(lambda x: x.split('_rep')[0], samples)
####MONACO#####
#celltypes =  {'CD4T_cells': ['Terminal_effector_CD4_T_cells', 'T_regulatory_cells', 'Th1_cells', 'Th2_cells', 'Th1Th17_cells', 'Th17_cells', 'Naive_CD4_T_cells', 'Follicular_helper_T_cells'], 'CD8T_cells': ['Central_memory_CD8_T_cell','Effector_memory_CD8_T_cells','Naive_CD8_T_cells','Terminal_effector_CD8_T_cells'], 'VdT_cells':['Vd2_gd_T_cells','Non-Vd2_gd_T_cells'], 'B_cells':['Exhausted_B_cells', 'Naive_B_cells','Non-switched_memory_B_cells','Plasmablasts', 'Switched_memory_B_cells'], 'NK_clles': ['Natural_killer_cells'],'Monocytes': ['Classical_monocytes', 'Intermediate_monocytes','Non_classical_monocytes'], 'Neutrophils': ['Low-density_neutrophils'], 'Basophils': ['Low-density_basophils'], 'Dendritic_cells':['Plasmacytoid_dendritic_cells', 'Myeloid_dendritic_cells'], 'Progenitor': ['Progenitor_cells'], 'MAIT':['MAIT_cells']}
#celltypes = {'B_cell': ['Exhausted_B_cells', 'Naive_B_cells', 'Non-switched_memory_B_cells', 'Plasmablasts', 'Switched_memory_B_cells'], 'T_CD4':['Terminal_effector_CD4_T_cells', 'T_regulatory_cells', 'Th1_cells', 'Th2_cells', 'Th1Th17_cells', 'Th17_cells', 'Naive_CD4_T_cells', 'Follicular_helper_T_cells'], 'T_CD8':['Central_memory_CD8_T_cell','Effector_memory_CD8_T_cells','Naive_CD8_T_cells','Terminal_effector_CD8_T_cells'], 'NK_cell':['Natural_killer_cells'],'Monocyte':['Classical_monocytes', 'Intermediate_monocytes','Non_classical_monocytes'], 'Neutrophil': ['Low-density_neutrophils']}
celltypes = {'B_cell': ['Naive_B_cells'], 'Th1':['Th1_cells'], 'Th2': ['Th2_cells'], 'Naive-CD4T':['Naive_CD4_T_cells'],'CD8T-TEM':['Effector_memory_CD8_T_cells'],'Naive-CD8T':['Naive_CD8_T_cells'], 'NK_cell':['Natural_killer_cells'],'Monocyte':['Classical_monocytes'], 'Neutrophil': ['Low-density_neutrophils']}

####SONG####
#celltypes = {'B_cell': ['Non-switched_memory_B_cells', 'Plasmablasts', 'Switched_memory_B_cells'], 'T_CD4':['Terminal_effector_CD4_T_cells', 'T_regulatory_cells', 'Th1_cells', 'Th2_cells', 'Th1Th17_cells', 'Th17_cells', 'Naive_CD4_T_cells', 'Follicular_helper_T_cells'], 'T_CD8':['Central_memory_CD8_T_cell','Effector_memory_CD8_T_cells','Naive_CD8_T_cells','Terminal_effector_CD8_T_cells'], 'NK_cell':['Natural_killer_cells'],'Monocyte':['Classical_monocytes', 'Intermediate_monocytes','Non_classical_monocytes'], 'Neutrophil': ['Low-density_neutrophils']}
#celltypes = {'B_cell': ['0604_Bcell'], 'Naive-CD4T':['G2_CD4-T-Naive'],'Th1':['G3_Th1'],'Th2':['G3_Th2'], 'CD8T-TEM':['G6_CD8-TEM'],'Naive-CD8T':['G6_CD8-Tnaive'], 'NK_cell':['0604_Nkcell'],'Monocyte':['0604_monocyte'], 'Neutrophil': ['0604_Neutrophil']}


#metaData.write('sampels\tcelltypes\n')
#for cell in celltypes.keys():
#	sams = celltypes[cell]
#	for sample in sams:
#		reps = list(set(filter(lambda x: sample == x.split('_rep')[0], lines[0].strip().split('\t')[1:])))
#		for rep in reps:
#			metaData.write(rep+'\t'+cell+'\n')
#metaData.close()
matdic = dict()
for i in xrange(1,len(lines)):
	line = lines[i].strip().split('\t')
	geneid = line[0]; exp = line[1:]
	if not matdic.has_key(geneid):
		matdic[geneid] = dict()
	for j in xrange(len(samples)):
		sam = samples[j]; express = exp[j]
		for type in celltypes.keys():
			cells = celltypes[type]
#			if sam in cells:
			if type in sam:
#				print type, sam
#				if not matdic[geneid].has_key(type):
				if not matdic[geneid].has_key(sam):
#					matdic[geneid][type] = []
					matdic[geneid][sam] = 0
#				matdic[geneid][type].append(express)
				matdic[geneid][sam] = express
#print matdic
#exit()
#		print sam, express
#		if not matdic[geneid].has_key(samples[j]):
#			matdic[geneid][samples[j]] = []
#		matdic[geneid][samples[j]].append(exp[j])
#	print matdic
#	exit()
sortSam = map(lambda x: x, sorted(list(set(samples))))
#sortSam = sorted(celltypes.keys())
print sortSam
#exit()
outfile.write('ID\t' + '\t'.join(sortSam)+'\n')
for geneid in matdic.keys():
	samdic = matdic[geneid]
	outfile.write(geneid)
	for name in sorted(samdic.keys()):
#		print name
#		print name, sorted(samdic[name])
#		median = getMedian(sorted(samdic[name]))
		median = samdic[name]
#		print median
		outfile.write('\t'+str(median))
	outfile.write('\n')
outfile.close()

