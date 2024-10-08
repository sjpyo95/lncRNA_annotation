import sys, os
import module as mdl


gtf_file = sys.argv[1]
filename = gtf_file.split('/')[-1].split('.gtf')[0]
outputdir = sys.argv[2]

gtf = mdl.getGtf2(gtf_file)

lowGtf = dict()
hiGtf = dict()
lowN = 0
hiN = 0
grades = dict()
for chr in gtf.keys():
	genes = gtf[chr]
	if not lowGtf.has_key(chr):
		lowGtf[chr] = []
	if not hiGtf.has_key(chr):
		hiGtf[chr] = []
	for gene in genes:
		grade = gene.elseinfo()['grade']
		if not grades.has_key(grade):
			grades[grade] = 0
		grades[grade] += 1
		if grade == 'NA': continue
		if float(grade) < 6:
			lowN += 1
			lowGtf[chr].append(gene)
		elif float(grade) >= 6:
			hiN += 1
			hiGtf[chr].append(gene)
		else:
			print grade

lowoutfile = outputdir + filename + '.low-conf.gtf'
hioutfile = outputdir + filename + '.high-conf.gtf'
#gradestat = open(outputdir + 'grade_score.stat', 'w')

#gradestat.write('GRADE\tGENES\n')

#for gr in grades.keys():
#	gradestat.write(gr+'\t'+str(grades[gr])+'\n')
#gradestat.close()
print 'low-confidence genes No. :',lowN
print 'high-confidence genes No. :',hiN
mdl.writeGtf3(lowGtf, lowoutfile)
mdl.writeGtf3(hiGtf, hioutfile)

