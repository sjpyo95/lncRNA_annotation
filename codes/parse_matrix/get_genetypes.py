import sys, os
import module as mdl

mtx_file = sys.argv[1]
utype = sys.argv[2]
out_file = sys.argv[3]

mtx, samples = mdl.getinfoMatrix(mtx_file)
allgenetypes = [x[2] for x in mtx.keys()]
print list(set(allgenetypes))
nmtx = dict()
for key in mtx.keys():
	type = key[2]
	if utype == type:
		nmtx[key] = mtx[key]

	elif utype == 'lncRNA':
		if type == 'Known_lncRNA':
			nmtx[key] = mtx[key]
		elif type == 'Novel_lncRNA':
			nmtx[key] = mtx[key]

	elif utype == 'PCG_lncRNA':
		if type == 'PCG':
			nmtx[key] = mtx[key]
		elif type == 'Known_lncRNA':
			nmtx[key] = mtx[key]
		elif type == 'Novel_lncRNA':
			nmtx[key] = mtx[key]
	else:
		print 'No matching genetype in this matrix. Check the genetype name.'
		exit()
mdl.writeinfoMatrix(nmtx, samples, out_file)
