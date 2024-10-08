import sys
sys.path.append('/home/jihoon/')
import settings as st

ref_gtf='/home/sjpyo/new/data/anno/human/GENCODE/lift37/gencode.v34lift37.annotation.gtf'

input_gtf='/home/sjpyo/new/immune/human/human_immune_analysis/ref/monaco_ysu_M0.novel_lncRNAs.gene.newID.gtf'

ref_dic=st.get_gtf_gene(ref_gtf,dict())
novel_dic=st.get_gtf_gene(input_gtf,dict())

n=0

fileout=open('/home/sjpyo/new/immune/human/human_immune_analysis/ref/check_novel.txt','w')

for chrom in st.ms_chrs:
	if not chrom in ref_dic.keys() or not chrom in novel_dic.keys():
		continue
	for ref_gene in ref_dic[chrom]:
		check=False
		k=0
		for ref_trx in ref_gene.transcripts():
			for ref_exon in ref_trx.exons():
				ref_start=ref_exon.start()
				ref_end=ref_exon.end()
				ref_strand=ref_exon.strand()
				for novel_gene in novel_dic[chrom]:
					for novel_trx in novel_gene.transcripts():
						for novel_exon in novel_trx.exons():
							novel_start=novel_exon.start()
							novel_end=novel_exon.end()
							novel_strand=novel_exon.strand()
							if ref_start==novel_start and ref_end==novel_end:
								if ref_strand!=novel_strand:
									k+=1
									if k<2: continue
									else:
										n+=1
										fileout.write(ref_gene.genename() + '\t' + novel_gene.genename() + '\n')
										check=True
							if check: break
						if check: break
					if check: break
				if check: break
			if check: break

print n
			
