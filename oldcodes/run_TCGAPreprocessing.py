import sys, os
import config as cf
import commands
import time
inputdir = sys.argv[1]
infofile = sys.argv[2]
outputdir = sys.argv[3]
def parseInfo(infile):
	filein = open(infile, 'r')
	lines = filein.readlines(); filein.close()
	infodic = dict()
	for i in xrange(1,len(lines)):
		line = lines[i].strip().split('\t')
		fileid = line[1]; caseid = line[2]; condi = line[28]; 
		newname = caseid + '_' + ''.join(condi.split(' '))
#		print fileid
#		print caseid
#		print condi
#		exit()
		infodic[fileid] = newname
	return infodic

infodic = parseInfo(infofile)
samples = filter(lambda x: x, os.listdir(inputdir))
for i in xrange(len(samples)):
	sample = samples[i]
	if sample in infodic.keys():
		newname = infodic[sample]
		soutputdir = outputdir + newname + '/'
		if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	sinputdir = inputdir + sample + '/'
	targz = sinputdir + filter(lambda x: '.tar.gz' in x, os.listdir(sinputdir))[0]
	comm = 'tar -xvf ' + targz + ' -C '+ soutputdir
	print comm + '\n'
	commands.getoutput(comm)
	
