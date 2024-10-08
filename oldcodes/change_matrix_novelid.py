import sys
import module as md

matrixfile = open(sys.argv[1], 'r')
gtffile = sys.argv[2]
outfile = open(sys.argv[3], 'w')

gtf = sum(md.getGene(gtffile).values(), [])
ids = dict()
for gene in gtf:
	oldid = gene.attri()
	newid = gene.geneid()
	ids[oldid] = newid

lines = matrixfile.readlines(); matrixfile.close()

for i in xrange(len(lines)):
	line = lines[i].strip()
	if line.startswith('ID'):
		col = line
		outfile.write(col+'\n')
	else:
		id = line.split('\t')[0]
		if id.split('TCELL-')[-1] in ids.keys():
			newid = ids[id.split('TCELL-')[-1]]
			newline = newid+'\t'+newid+'\t'+'\t'.join(line.split('\t')[2:])
			outfile.write(newline+'\n')
		else:
			outfile.write(line+'\n')
outfile.close()
#matrix = md.getMatrix(matrixfile)
#
#newmat = dict()
#for idsym in matrix:
#	id = idsym[0]; symb = idsym[1]
#	express = matrix[idsym]
#	if id.split('TCELL-')[-1] in ids.keys():
#		newid = ids[id.split('TCELL-')[-1]]
#		newmat[(newid,newid)] = express
#	else:
#		newmat[idsym] = express
#md.writeMatrix(newmat, outfile)
