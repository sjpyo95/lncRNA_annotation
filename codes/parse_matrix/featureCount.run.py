##conda activate express
import argparse, os
import config as cf

def main(thread, inputdir, refFasta, refAnno, refAnnoFormat, paired, strand, feature, metaFeature, multiread, outputdir, multi):
	if not os.path.exists(outputdir): os.makedirs(outputdir)
	if multi:
		sinputdirs = filter(lambda x: '.txt' not in x and '.out' not in x and '.xls' not in x, os.listdir(inputdir))
		for i in range(len(sinputdirs)):
			sinputdir = inputdir + '/' + sinputdirs[i] + '/'
			bamdir = sinputdir + '/1.alignment/'
			bamfile = bamdir + filter(lambda x: '.bam' in x and '.bai' not in x, os.listdir(bamdir))[0]
			outfile = outputdir + sinputdirs[i] + '.count.txt'
			featureCounts(bamfile, outfile, thread, refFasta, refAnno, refAnnoFormat, paired, strand, feature, metaFeature, multiread, outputdir)
	else:
		bamfiles = [inputdir + x for x in filter(lambda x: '.bam' in x and '.bai' not in x, os.listdir(inputdir))]
		for bamfile in bamfiles:
			outfile = outputdir + bamfile.split('/')[-1].split('_Aligned')[0] + '.count.txt'
			featureCounts(bamfile, outfile, tread, refFasta, refAnno, refAnnoFormat, paired, strand, feature, metaFeature, multiread, outputdir)

			
def featureCounts(bamfile, outfile, thread, refFasta, refAnno, refAnnoFormat, paired, strand, feature, metaFeature, multiread, outputdir):	
	fcCmd = 'featureCounts -T ' + thread + ' -t ' + feature + ' -s ' + strand + ' -g ' + metaFeature + ' -o ' + outfile
	if paired:
		fcCmd += ' -p '
		if multiread:
			fcCmd += ' -M '
		if refFasta:
			fcCmd += ' -G ' + refFasta
		if refAnno:
			fcCmd += ' -a ' + refAnno + ' -F ' + refAnnoFormat
		fcCmd += ' ' + bamfile
#		print fcCmd
		cf.run(fcCmd, True)


if __name__ == '__main__':
	ps = argparse.ArgumentParser(description = 'Run FeatureCount to quantify expression')
	ps.add_argument('-d', help = "A directory that contains BAM and indexed file", dest = 'inputdir', required = True)
	ps.add_argument('-T', help = "Number of threads", dest = 'thread')
	ps.add_argument('-p', help = "Indicate if reads are paired-end", action = 'store_true', default = False, dest = 'paired')
	ps.add_argument('-s', help = "Indicate if reads are strand-specific. 0 (unstranded), 1 (stranded) and 2 (reversely stranded)", default = 0, dest = 'strand')
	ps.add_argument('-t', help = "Specify the feature type(s)", dest = 'feature')
	ps.add_argument('-g', help = "Specify the attribute type used to group features into meta-features (eg. gene_id)", dest = 'metaFeature')
	ps.add_argument('-G', help = "The name of a FASTA-format file that contains the reference sequences used in read mapping", default = False, dest = 'refFasta')
	ps.add_argument('-M', help = "Multimapping reads will be counted", action = 'store_true', default = False, dest = 'multiread')
	ps.add_argument('-a', help = "Annotation file", default = False, dest = 'refAnno')
	ps.add_argument('-F', help = "Format of annotation file (eg. GTF)", default  = 'GTF',  dest = 'refAnnoFormat')
	ps.add_argument('-o', help = "Output directory", dest = 'outputdir', required = True)
	ps.add_argument('-multirun', help = "Run in multiple directories", action = 'store_true', default = False, dest = 'multi')
	
	ag = ps.parse_args()
	main(ag.thread, ag.inputdir, ag.refFasta, ag.refAnno, ag.refAnnoFormat, ag.paired, ag.strand, ag.feature, ag.metaFeature, ag.multiread, ag.outputdir, ag.multi)
