import sys
import module as mdl

matrixfile = sys.argv[1]
pcgfile = sys.argv[2]
lncfile = sys.argv[3]
novelfile = sys.argv[4]

matrix = mdl.getMatrix(matrixfile)
pcg = mdl.getGene(pcgfile)
lnc = mdl.getGene(lncfile)
novel = mdl.getGene(novelfile)

for idsymb in matrix:

