
# This script will submit the jobs for processing all of the tnseq data

cd /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/Talia
####################################
#First directory to manage is recent tnseq run
####################################
all_reads=/uufs/chpc.utah.edu/common/home/karasov-group1/SRA/weigel_tnseq/TnSeq_2020_08_28/

#  
for dir in `ls $all_reads | grep p7_g9`; \
	    do sample=`echo $dir | sed 's/illumina_ST-J00101_flowcellA_SampleId//g' | sed 's/_RunId0183_LaneId.//g'`; 

# Need to specify index for each run. If uninfected do against DC3000, else choose appropriate
	    if [[ "$sample" == *"uninfected"* ]]; then
	    index="DC3000"; fi
	    if [[ "$sample" == *"DC3000"* ]]; then
	    index="DC3000";
	    fi
	    if [[ "$sample" == *"p7_g9"* ]]; then
	    index="p7_g9";
	    fi
	    if [[ "$sample" == *"p25_c2"* ]]; then
	    index="p25_c2";
	    fi
	    echo $sample; echo $index; 

#now for submitting the script
	    sbatch --export=index=$index,read_direc=$all_reads/$dir,sample=$sample /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/scripts/tnseq_pseudomonas/scripts/generic_adapter_trimming.sh

	    done

####################################
# Peter's tnseq run
####################################

	    all_reads=/uufs/chpc.utah.edu/common/home/karasov-group1/SRA/weigel_tnseq/Tnseq_2019_10_24_Peter 
	    index="p25_c2"

	    for dir in `ls $all_reads`; \
		    do sample=`echo $dir | sed 's/illumina_ST-J00101_flowcellA_SampleId//g' | sed 's/_RunId0159_LaneId.//g'`; 
		    echo $sample; echo $index; 
		    sbatch --export=index=$index,read_direc=$all_reads/$dir,sample=$sample /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/scripts/tnseq_pseudomonas/scripts/generic_adapter_trimming.sh
		    done
