#!/usr/bin/python

# CAFE
sampling_number = 1000000
position_number = 25000000

MaxEntScan = '/home/bhyou/cafe/src/build/MaxEntScan/'

# CASOL
length_threshold = int(200)

hm_cpc_cutoff = float(-0.3); hm_te_cutoff = float(0.71)
ms_cpc_cutoff = float(-0.2); ms_te_cutoff = float(0.08)

hm_rnaseq = '/home/bhyou/project/lincRNA/data_set/mRNAseq/HeLa/GSM546921_trimmed_filtered/accepted_hits.bam'
hm_rnaseq_readtype = 'single-end'
hm_rnaseq_strand = 'strandSpecific'
hm_rnaseq_library = 'fr-unstranded'

hm_rprofile = '/home/bhyou/project/lincRNA/data_set/RIBOseq/HeLa/GSM546920_trimmed_filtered/accepted_hits.bam'
hm_rpro_readtype = 'single-end'
hm_rpro_strand = 'strandSpecific'
hm_rpro_library = 'fr-unstranded'

ms_rnaseq = '/home/bhyou/project/lincRNA/data_set/mRNAseq/Neutrophil/GSM546988_trimmed_filtered/accepted_hits.bam'
ms_rnaseq_readtype = 'single-end'
ms_rnaseq_strand = 'strandSpecific'
ms_rnaseq_library = 'fr-unstranded'

ms_rprofile = '/home/bhyou/project/lincRNA/data_set/RIBOseq/Neutrophil/GSM546987_trimmed_filtered/accepted_hits.bam'
ms_rpro_readtype = 'single-end'
ms_rpro_strand = 'strandSpecific'
ms_rpro_library = 'fr-unstranded'

# GENCODE
hm_gencode_others = '/home/hyeonlee17/BIG.Projects/lncRNAs/Goldstandard/updated.outputs/hg19/others/gencode.v19.annotation.UpdateCAGEFilter.UpdatePolyAFilter.others.gtf'
hm_gencode_lincRNAs = '/home/hyeonlee17/BIG.Projects/lncRNAs/Goldstandard/updated.outputs/hg19/lncRNAs/gencode.v19.annotation.UpdateCAGEFilter.UpdatePolyAFilter.lncRNAs.gtf'
ms_gencode_others = '/home/hyeonlee17/BIG.Projects/lncRNAs/Goldstandard/updated.outputs/mm9/others/gencode.vM1.annotation.UpdateCAGEFilter.UpdatePolyAFilter.others.gtf'
ms_gencode_lincRNAs = '/home/hyeonlee17/BIG.Projects/lncRNAs/Goldstandard/updated.outputs/mm9/lncRNAs/gencode.vM1.annotation.UpdateCAGEFilter.UpdatePolyAFilter.lncRNAs.gtf'

# RefFlat
hm_refFlat_others = '/home/hyeonlee17/BIG.Projects/lncRNAs/Goldstandard/updated.outputs/hg19/refFlat/others/refFlat_hg19_20130909.UpdateCAGEFilter.UpdatePolyAFilter.others.gtf'
hm_refFlat_lincRNAs = '/home/hyeonlee17/BIG.Projects/lncRNAs/Goldstandard/updated.outputs/hg19/refFlat/ncRNAs/refFlat_hg19_20130909.UpdateCAGEFilter.UpdatePolyAFilter.lncRNAs.gtf'
ms_refFlat_others = '/home/hyeonlee17/BIG.Projects/lncRNAs/Goldstandard/updated.outputs/mm9/refFlat/others/refFlat_mm9_20130909.UpdateCAGEFilter.UpdatePolyAFilter.others.gtf'
ms_refFlat_lincRNAs = '/home/hyeonlee17/BIG.Projects/lncRNAs/Goldstandard/updated.outputs/mm9/refFlat/ncRNAs/refFlat_mm9_20130909.UpdateCAGEFilter.UpdatePolyAFilter.lncRNAs.gtf'

#MiTranscriptome
hm_mitrans_others = '/home/bhyou/cafe/dataset/annotation/mitranscriptome/human/hg19/mitranscriptome.v2.others.gtf'
hm_mitrans_lincRNAs = '/home/bhyou/cafe/dataset/annotation/mitranscriptome/human/hg19/mitranscriptome.v2.lincRNAs.gtf'

#Ensemble
ms_ensemble_others = '/home/bhyou/cafe/dataset/annotation/ensemble/mouse/mm9/Mus_musculus.NCBIM37.67.others.gtf'
ms_ensemble_lincRNAs = '/home/bhyou/cafe/dataset/annotation/ensemble/mouse/mm9/Mus_musculus.NCBIM37.67.lincRNAs.gtf'

# Human
hm_chrs = map(lambda x: 'chr' + x, map(str, range(1, 23)) + ['X', 'Y'])
#hm_chrs = map(lambda x: 'chr' + x, map(str, range(1, 23)) + ['X', 'Y', 'M'])
hm_fasta = '/export/home/sjpyo/lncRNA_anno_profiling_exercise/hepg2/data/genome/ucsc.hg19.fasta'
hm_nibDir = '/Data_Set/Genome/human/hg19/blat/'
hm_genome = '/home/bhyou/cafe/dataset/info/hg19.genome'
hm_intronmin = 61; hm_intronmax = 25099
#hm_intronmin = 61; hm_intronmax = 265006
#hm_intronmin = 85; hm_intronmax = 25099

hm_gencode_protein = '/home/bhyou/project/Prototype/dataset/cuffcompare/human/gencode.v19.annotation_N_protein.gtf'
hm_gencode_lincRNA = '/home/bhyou/project/Prototype/dataset/cuffcompare/human/gencode.v19.annotation_N_lincRNA.gtf'

# CAGEseq
hm_fpseqAnnotation = '/home/bhyou/cafe/dataset/annotation/cage/human/hg19/fantom/hg19.cage_peak_coord_robust.bed'
hm_fpseqThreshold = float(0.01); hm_fpseqInterval = int(3000)

# PolyAseq
hm_tpseqAnnotation = '/home/bhyou/cafe/dataset/annotation/polya/human/hg19/hg19_all.15.sumCM.bed'
hm_tpseqHeLa = '/home/bhyou/cafe/dataset/annotation/polya/human/hg19/GSM1268942_hela.15.sumCM.bed'
hm_tpseqHEK293 = '/home/bhyou/cafe/dataset/annotation/polya/human/hg19/GSM1268943_hek293.15.sumCM.bed'
hm_tpseqHuh7 = '/home/bhyou/cafe/dataset/annotation/polya/human/hg19/GSM1268944_huh7.15.sumCM.bed'
hm_tpseqIMR90 = '/home/bhyou/cafe/dataset/annotation/polya/human/hg19/GSM1268945_imr90.15.sumCM.bed'
hm_tpseqThreshold = float(0.01); hm_tpseqInterval = int(5000)

# HeLa
hela_np_gtfFile = '/home/bhyou/cafe/dataset/cufflinks/human/HeLa/unstranded/parameter/transcripts.gtf'
hela_sp_gtfFile = '/home/bhyou/cafe/dataset/cufflinks/human/HeLa/stranded/parameter/transcripts.gtf'
hela_np_bamFile = '/home/bhyou/cafe/dataset/tophat/mRNAseq/human/HeLa/unstranded/HeLa_Np.sorted.bam'
hela_sp_bamFile = '/home/bhyou/cafe/dataset/tophat/mRNAseq/human/HeLa/stranded/HeLa_Sp.sorted.bam'
hela_all_bamFile = '/home/bhyou/cafe/dataset/tophat/mRNAseq/human/HeLa/HeLa_all.sorted.bam'
hela_outputdir = '/home/bhyou/cafe/results/HeLa/'

# HeLaV2
helaV2_np_gtfFile = '/home/bhyou/cafe/dataset/cufflinks/human/HeLaV2/unstranded/parameter/transcripts.gtf'
helaV2_sp_gtfFile = '/home/bhyou/cafe/dataset/cufflinks/human/HeLaV2/stranded/parameter/transcripts.gtf'
helaV2_np_bamFile = '/home/bhyou/cafe/dataset/tophat/mRNAseq/human/HeLaV2/unstranded/HeLa_Np.sorted.bam'
helaV2_sp_bamFile = '/home/bhyou/cafe/dataset/tophat/mRNAseq/human/HeLaV2/stranded/HeLa_Sp.sorted.bam'
helaV2_all_bamFile = '/home/bhyou/cafe/dataset/tophat/mRNAseq/human/HeLaV2/HeLa_all.sorted.bam'
helaV2_outputdir = '/home/bhyou/cafe/results/HeLaV2/'
helaV2_cnp_gtfFile = '/home/bhyou/cafe/results/HeLaV2/combine/default/transcripts.gtf'
helaV2_mnp_gtfFile = '/home/bhyou/cafe/results/HeLaV2/cuffmerge/transcripts.gtf'

hela_strand_threshold = float(0.95); hela_read_threshold = int(5)

hela_protein_rpkm_1 = '/home/bhyou/cafe/dataset/annotation/cuffcompare/human/HeLa/gencode.v19_protein_rpkm_1.gtf'
hela_protein_rpkm_01 = '/home/bhyou/cafe/dataset/annotation/cuffcompare/human/HeLa/gencode.v19_protein_rpkm_0.1.gtf'
hela_lncRNA_rpkm_1 = '/home/bhyou/cafe/dataset/annotation/cuffcompare/human/HeLa/gencode.v19_lncRNA_rpkm_1.gtf'
hela_lncRNA_rpkm_01 = '/home/bhyou/cafe/dataset/annotation/cuffcompare/human/HeLa/gencode.v19_lncRNA_rpkm_0.1.gtf'

hela_raw_inputdir = '/home/bhyou/cafe/results/HeLaV2/updating_cps/npsp/all/'
hela_length_inputdir = '/home/bhyou/cafe/results/HeLaV2/casol/length/'
hela_gencode_inputdir = '/home/bhyou/cafe/results/HeLaV2/casol/overlap/others/gencode/'
hela_refFlat_inputdir = '/home/bhyou/cafe/results/HeLaV2/casol/overlap/others/refFlat/'
hela_mitrans_inputdir = '/home/bhyou/cafe/results/HeLaV2/casol/overlap/others/mitranscriptome/'

hela_cpc_file = '/home/bhyou/cafe/results/HeLaV2/casol/cpc/results.txt'
hela_te_file = '/home/bhyou/cafe/results/HeLaV2/casol/te/results.txt'

hela_filter_inputdir = '/home/bhyou/cafe/results/HeLaV2/casol/filter/'
hela_lncRNA_inputdir = '/home/bhyou/cafe/results/HeLaV2/casol/overlap/lincRNAs/gencode/'
hela_lncRNA_inputdir2 = '/home/bhyou/cafe/results/HeLaV2/casol/overlap/lincRNAs/refFlat/'

# Esophageal_carcinoma
ec_raw_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/dataset/cufflinks/'
ec_bam_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/dataset/tophat/novel/'
fn_all_bamFile = '/home/bhyou/cancer/esophageal_carcinoma/dataset/tophat/novel/fn.bam'
fc_all_bamFile = '/home/bhyou/cancer/esophageal_carcinoma/dataset/tophat/novel/fc.bam'
ec_outputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/'

ec_exj_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/updating_exj/'
ec_cuf_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/cufflinks/'
ec_cuf_inputgtf = '/home/bhyou/cancer/esophageal_carcinoma/results/cuffmerge/transcripts.gtf'
ec_tss_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/updating_tss/'

ec_cps_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/updating_cps/'
ec_length_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/length/'
ec_gencode_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/overlap/others/gencode/'
ec_refFlat_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/overlap/others/refFlat/'
ec_mitrans_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/overlap/others/mitranscriptome/'

ec_cpc_file = '/home/bhyou/cancer/esophageal_carcinoma/results/cpc/results.txt'
ec_te_file = '/home/bhyou/cancer/esophageal_carcinoma/results/te/results.txt'

ec_filter_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/filter/'
ec_lincRNA_inputdir = '/home/bhyou/cancer/esophageal_carcinoma/results/overlap/lincRNAs/gencode/'
ec_lincRNA_inputdir2 = '/home/bhyou/cancer/esophageal_carcinoma/results/overlap/lincRNAs/refFlat/'

# ENCODE
en_inputdir = '/home/bhyou/cafe/dataset/encode/'
en_cpc_file = '/home/bhyou/outputdir/cocoa/casol/cpc/tmp/result.txt'
en_te_file = '/home/bhyou/outputdir/cocoa/casol/te/tmp/result.txt'

# CATHOLIC
ca_inputdir = '/home/bhyou/cancer/hepatocellular_carcinoma/'
ca_cpc_file = '/home/bhyou/cancer/hepatocellular_carcinoma/results/cuffmerge/sample/casol/cpc/tmp/result.txt'
ca_te_file = '/home/bhyou/cancer/hepatocellular_carcinoma/results/cuffmerge/sample/casol/te/tmp/result.txt'

# Mouse
ms_chrs = map(lambda x: 'chr' + x, map(str, range(1, 20)) + ['X', 'Y'])
ms_fasta = '/Data_Set/Genome/mouse/mm9/mm9_fasta/mm9.fa'
ms_nibDir = '/Data_Set/Genome/mouse/mm9/blat/'
ms_genome = '/home/bhyou/cafe/dataset/info/mm9.genome'
ms_intronmin = 52; ms_intronmax = 18914
#ms_intronmin = 52; ms_intronmax = 240764
#ms_intronmin = 95; ms_intronmax = 18914

# CAGEseq
ms_fpseqAnnotation = '/home/bhyou/cafe/dataset/annotation/cage/mouse/mm9/fantom/mm9.cage_peak_coord_robust.bed'
ms_fpseqThreshold = float(0.01); ms_fpseqInterval = int(2000)

# PolyAseq
ms_tpseqAnnotation = '/home/bhyou/cafe/dataset/annotation/polya/mouse/mm9/mm9_all.15.sumCM.bed'
ms_tpseqHeart = '/home/bhyou/cafe/dataset/annotation/polya/mouse/mm9/GSM1268946_heart_wt_pA.15.sumCM.bed'
ms_tpseqMuscle = '/home/bhyou/cafe/dataset/annotation/polya/mouse/mm9/GSM1268947_muscle_wt_pA.15.sumCM.bed'
ms_tpseqLiver = '/home/bhyou/cafe/dataset/annotation/polya/mouse/mm9/GSM1268948_liver_wt_pA.15.sumCM.bed'
ms_tpseqLung = '/home/bhyou/cafe/dataset/annotation/polya/mouse/mm9/GSM1268949_lung_wt_pA.15.sumCM.bed'
ms_tpseqWat = '/home/bhyou/cafe/dataset/annotation/polya/mouse/mm9/GSM1268950_wat_wt_pA.15.sumCM.bed'
ms_tpseqKidney = '/home/bhyou/cafe/dataset/annotation/polya/mouse/mm9/GSM1268951_kidney_wt_pA.15.sumCM.bed'
ms_tpseqES = '/home/bhyou/cafe/dataset/annotation/polya/mouse/mm9/GSM1268958_mESC_pA.15.sumCM.bed'
ms_tpseq3T3 = '/home/bhyou/cafe/dataset/annotation/polya/mouse/mm9/GSM1268959_3T3_pA.15.sumCM.bed'
ms_tpseqThreshold = float(0.01); ms_tpseqInterval = int(3000)

# mES
mes_np_gtfFile = '/home/bhyou/cafe/dataset/cufflinks/mouse/mES/unstranded/parameter/transcripts.gtf'
mes_sp_gtfFile = '/home/bhyou/cafe/dataset/cufflinks/mouse/mES/stranded/pe/J1_x/transcripts.gtf'
mes_np_bamFile = '/home/bhyou/cafe/dataset/tophat/mRNAseq/mouse/mES/unstranded/J1_x/accepted_hits.bam'
mes_sp_bamFile = '/home/bhyou/cafe/dataset/tophat/mRNAseq/mouse/mES/stranded/pe/J1_x/accepted_hits.bam'
mes_all_bamFile = '/home/bhyou/cafe/dataset/tophat/mRNAseq/mouse/mES/mESC_all.sorted.bam'
mes_outputdir = '/home/bhyou/cafe/results/mES/'
mes_cnp_gtfFile = '/home/bhyou/cafe/results/mES/combine/default/transcripts.gtf'
mes_mnp_gtfFile = '/home/bhyou/cafe/results/mES/cuffmerge/transcripts.gtf'

mes_strand_threshold = float(0.95); mes_read_threshold = int(5)

mes_protein_rpkm_1 = '/home/bhyou/cafe/dataset/annotation/cuffcompare/mouse/mES/gencode.vM1_protein_rpkm_1.gtf'
mes_protein_rpkm_01 = '/home/bhyou/cafe/dataset/annotation/cuffcompare/mouse/mES/gencode.vM1_protein_rpkm_0.1.gtf'
mes_lncRNA_rpkm_1 = '/home/bhyou/cafe/dataset/annotation/cuffcompare/mouse/mES/gencode.vM1_lncRNA_rpkm_1.gtf'
mes_lncRNA_rpkm_01 = '/home/bhyou/cafe/dataset/annotation/cuffcompare/mouse/mES/gencode.vM1_lncRNA_rpkm_0.1.gtf'

mes_raw_inputdir = '/home/bhyou/cafe/results/mES/updating_cps/npsp/all/'
mes_length_inputdir = '/home/bhyou/cafe/results/mES/casol/length/'
mes_gencode_inputdir = '/home/bhyou/cafe/results/mES/casol/overlap/others/gencode/'
mes_refFlat_inputdir = '/home/bhyou/cafe/results/mES/casol/overlap/others/refFlat/'
mes_ensemble_inputdir = '/home/bhyou/cafe/results/mES/casol/overlap/others/ensemble/'

mes_cpc_file = '/home/bhyou/cafe/results/mES/casol/cpc/results.txt'
mes_te_file = '/home/bhyou/cafe/results/mES/casol/te/results.txt'

mes_filter_inputdir = '/home/bhyou/cafe/results/mES/casol/filter/'
mes_lncRNA_inputdir = '/home/bhyou/cafe/results/mES/casol/overlap/lincRNAs/gencode/'
mes_lncRNA_inputdir2 = '/home/bhyou/cafe/results/mES/casol/overlap/lincRNAs/refFlat/'
