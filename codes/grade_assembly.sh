#!/usr/bin/bash
set -euo pipefail

# old_PYTHONPATH=$PYTHONPATH
# PYTHONPATH="/home/sjlee/programs/Python-3.9.4/lib/python3.9/site-packages"

### Confidence scores
# 256 <= FPK       --> 8+ : High confidence
# 128 <= FPK < 256 --> 7  : Intermediate confidence
# 64  <= FPK < 128 --> 6  : Intermediate confidence
#        FPK < 64  --> 5- : Low confidence

### Step 1. Receives a StringTie assembly/meta-assembly and BAM(s) used to create the assembly as inputs.
## Multiple BAM inputs should be separated by commas.
usage() {
        echo "Usage: ${0} [-h] -g <in_gtf> [-s <0|1|2>] [-b <bam1,...,bamN>]
	-g: (GTF) StringTie assembly/meta-assembly to grade (graded GTF with only -g given simply outputs grade stats).
	-s: (strand) strand parameter for featureCounts (0: unstranded / 1: stranded / 2: reversely stranded)
	-b: (BAM) BAM file(s) used for StringTie assembly/meta-assembly." >&2;	
exit 1;
}

while getopts ":g:s:b:h" opt; do
        case "${opt}" in
                g)
                        in_gtf=${OPTARG}
                        ;;
                s)
                        strand=${OPTARG}
                        ;;
                b)
                        bam_files=${OPTARG}
                        ;;
                h)
			usage
			;;
		:)
                        echo "ERROR: Option -${OPTARG} requires an argument"
                        usage
                        ;;
                \?)
			echo "Invalid option -${OPTARG}" >&2
                        usage
                        ;;
        esac
done

### Check required switches exist
#if [ -z "${in_gtf}" ] || [ -z "${strand}" ] || [ -z "${bam_files}" ]; then
#        usage
#fi
if [[ -z "${in_gtf}" ]]; then
	usage
fi

echo
echo "Input arguments"
echo "-- in_gtf: ${in_gtf}"
in_gtf_dir=$(dirname ${in_gtf})


### If -s and -b are given, add gene lines and grades to ${in_gtf}
if [ ${OPTIND} -eq 7 ]
then
	echo "-- strand: ${strand}"
	# echo "-- bam_files: $bam_files"
	bam_files_arr=("${bam_files//,/ }")
	echo "-- bam_files: ${bam_files_arr}"	

	### Step 2. Run a python3 script that takes the assembly GTF as input and outputs a new GTF with gene lines added.
	echo "$(which python3) /home/sjlee/projects/002_grade_assembly/src/python/GA_add_gene_lines.py ${in_gtf}"
	/share/apps/python3.4/bin/python3.4 /home/sjlee/projects/002_grade_assembly/src/python/GA_add_gene_lines.py ${in_gtf}
	in_gtf="${in_gtf/.gtf/.graded.gtf}"
	in_gtf_dir=$(dirname ${in_gtf})

	### Step 3. For each input BAM, run featureCounts with the new GTF given as reference.
	for bam in $bam_files_arr
	do
        	out_bam_dir=$(dirname ${bam})
        	out_fCounts="${out_bam_dir}/featureCounts.txt"
		#echo "/home/sjlee/programs/subread-2.0.2-Linux-x86_64/bin/featureCounts $bam -a $in_gtf -s $strand -o $out_fCounts -p -B -C"
		$(which featureCounts) ${bam} -a ${in_gtf} -s ${strand} -o ${out_fCounts} -p -B -C
	done

	### Step 4. Run a python3 script that takes the new GTF and featureCounts results as inputs and adds grades (confidence scores) to transcripts (calculated at the gene level).
	echo
	echo "Outputs"
	echo "-- out_gtf: ${in_gtf}"
	/share/apps/python3.4/bin/python3.4 /home/sjlee/projects/002_grade_assembly/src/python/GA_grade_assembly.py ${in_gtf} ${bam_files_arr}
fi

### Step 5. Write and plot grade stats of genes and transcripts within input GTF.
/share/apps/python3.4/bin/python3.4 /home/sjlee/projects/002_grade_assembly/src/python/GA_results.py ${in_gtf}
echo
echo "Outputs"
echo "-- out_num_stats: ${in_gtf_dir}/grade_assembly_results.num.stats" 
echo "-- out_num_stats_plot: ${in_gtf_dir}/grade_assembly_results.num.png"
echo "-- out_class_stats: ${in_gtf_dir}/grade_assembly_results.class.stats"


# ### Return PYTHONPATH back to normal
# export PYTHONPATH=$old_PYTHONPATH

