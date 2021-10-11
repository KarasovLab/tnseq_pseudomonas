#!/bin/bash

#SBATCH -t 1:00:00 # assign time limit

#SBATCH -J multicov # assign job name

#SBATCH -N 1 # node number

#SBATCH -o out%tnseq # assign standard output file name

#SBATCH -e err%tnseq # assign standard error file name

#SBATCH -n 1 # assign CPU number (the owner cluster should have 128 threads for 64 cores)

#SBATCH -A mpag-np  # group name (the owner node name)

#SBATCH -p mpag-shared-np  # microbe, plant and animal genetics (the partition name). Note: mpag-np (to have the entire node) or mpag-shared-np (to request a portion of the node)

#SBATCH --mail-type=begin,end # email you when begin, end, fail, requeue

#SBATCH --mail-user=aubrey.hawks@utah.edu # your email address

#SBATCH --mem=5000 # request memory in Gb  (the owner node have maximum 500 GB memory, assign it if you know the exact memory you want to use)

# The following lines are commented out as they work with the notch-peak freecycle
## SBATCH -A karasov-group1 # my group name

## SBATCH --partition notchpeak-freecycle # request notchpeak-freecycle for script testing





#edited now for p laurin's data this script takes the bam files and overlaps them with gene annotations. This can then go into edgeR

#go to the directory where we want the output
cd /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/Aubrey #should this go in my folder? 
#test on one file to get this working

#this loop iterates over the files in the directory with the tnseq outputs, looks for bam files that are not sorted, sorts them and appends 'sorted' to the name
#for sample in `ls /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/Talia | grep bam | grep -v sorted`;
#for sample in `ls /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/Talia/ | grep bam | grep -v sorted`;  
#	do
#	#need to sort and index the bam files
#		samtools sort /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/Talia/$sample > sorted.$sample
#		samtools index sorted.$sample
#done
#
##This block may not be necessary, depending on if the chromosomes are named consistently - there should be an error if there is an issue, but until then, I'll comment this out
##with the original gff  the chromosome name is wrong. We need to relabel the gff file
##sed -e 's/utg000001l_p25.C2/utg000001c:1.0-5963307.0_pilon/' \
##/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon_no_fasta.gff > plate25.c2_fin.gff
#
##some of the files are empty and will cause  index errors for bedtools, so we need to move them out of the directory:
#
for file in `find /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/Aubrey/ -type f -empty`;
        do mv $file /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/Talia/empty_bams/;
done
#
#and call bedtools (on all bams) with overlap with full genome
#We have multiple genomes, so we need to match each file to the correct reference
#I think that since I specified data/Aubrey as the directory, all outputs to this point should be here and I should only have sorted bam files

#echo "Starting p25"
#for p25
#bedtools multicov -bams ./sorted.*p25*.bam -bed /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/p25_c2/plate25.C2.annotation_no_fasta_rename.gff  > \
#p25c2_counts.gff


echo "Starting p7"
#for p7
bedtools multicov -bams ./sorted.TnSeq_p7_g9_*.bam -bed /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/p7_g9/annotation/plate7.G9.annotation.wOutFasta_rename.gff  > \
p7g9_counts.gff


#
#echo "Starting DC3000"
##for DC3000
#bedtools multicov -bams ./sorted.*DC3000*.bam -bed /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/DC3000_GCF_000007805.1/genomic.gff> \
#DC3000_counts.gff
#
#echo "Starting MgSO4"
##for MgSO4
#bedtools multicov -bams ./sorted.MgSO4*.bam -bed /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/DC3000_GCF_000007805.1/genomic.gff> \
#MgSO4_counts.gff
#
#
#echo "Starting NA"
##for NA
#bedtools multicov -bams ./sorted.NA*.bam -bed /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/DC3000_GCF_000007805.1/genomic.gff> \
#NA_counts.gff
#
#echo "Starting plate"
##for plate
#bedtools multicov -bams ./sorted.plate*.bam -bed /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/p25_c2/plate25.C2.annotation_no_fasta_rename.gff> \
#plate_counts.gff

#ls ./sorted.plate*.bam > plaurin_order.txt
#
#echo "Starting uninfected"
##for uninfected
#bedtools multicov -bams ./sorted.*uninfected*.bam -bed /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/D3000_GCF_000007805.1/genomic.gff> \
#uninfected_counts.gff

#map_ref = ()
#for sample in `ls /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/Aubrey | grep -v p25 | grep -v p7`;
#	do
#		map_ref += (sample)
#done

#for s in ${map_ref[@]};
#	do

#		echo $s
#done
#		bedtools multicov -bams $s  -bed /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/D3000_GCF_000007805.1/genomic.gff


#head p25c2_counts.gff
#sed 's/^.*locus_tag=//' counts_plaurin.gff > gene_counts.tab

