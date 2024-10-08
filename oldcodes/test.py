import sys
import module as mdl
tablefile = sys.argv[1]
geneLen = 10

def getfpkm(count, totalcount, length):
	return float((10**9)*count)/(totalcount*length)

def getTotalfpkm(countDic):
	totalcount = sum([x[-1] for x in countDic.values()])
	totalfpkm = 0
	for geneid in countDic.keys():
		sample = countDic[geneid]
		length = 10
		count = countDic[geneid][2]
		totalfpkm += getfpkm(count, totalcount, length)
	return totalfpkm

def gettpm(count, totalfpkm, length):
	fpkm = getfpkm(count, totalcount, length)
	return float(fpkm*10**6)/(totalfpkm)

A = {'gene1' : 1000, 'gene2' : 2000, 'gene3': 3000, 'gene4':4000, 'gene5':5000}
totalcount = sum(A.values())
totalfpkm = 0
for gene in A.keys():
	count = A[gene]
	fpkm = getfpkm(count, totalcount, 10)
	totalfpkm += fpkm
for gene in A.keys():
	count = A[gene]
	tpm = gettpm(count, totalfpkm, 10)
	print gene + '\t' + str(tpm)
