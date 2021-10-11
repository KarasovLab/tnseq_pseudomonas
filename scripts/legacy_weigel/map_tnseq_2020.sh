#!/bin/bash

# Modified from Manuela's script in Mar. 2020 (based off of her Jan. 2020 script)
# This is a shell script for mapping the Tn-insertion mutants to the P25.C2 Pseudomonas genome, in order to identify the location of insertions
# HiSeq Run 0159, lane 6, October 2019
# Script adapted from s1_general_adapter_trimming.sh from Talia
# A copy of this can be found in /ebio/abt6_projects7/small_projects/mhoelscher/Tn_mapping/s1_general_adapter_trimming_Manu.sh

# It is necessary to remove the transposon (pIT2) sequence from the samples, as well as the transposon (pIT2) sequence, as they will not map to the Pseudomonas genome
# https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapter
# we have 5 different ones: P5_Nex_Tn_amend_S502, P5_Nex_Tn_amend_S503, P5_Nex_Tn_amend_S504, P5_Nex_Tn_amend_S505, P5_Nex_Tn_amend_S506
# example:
# P5_Nex_Tn_amend_S502 AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCAGGACGCTACTTGTGTATAAGA
# The transposon part of those primers are the same for all, so not necessary to change the sequence in the script:
# pIT2
# CAGGACGCTACTTGTGTATAAGA


# # Establish read and processed directories, and suffix and prefix
read_direc=/ebio/abt6_projects8/Pseudomonas_mapping/data/raw_reads/hiseq_0159
processed_direc=/ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process #mneumann_tnseq/TnSeq_Oct19
suffix=_RunId0159_LaneId6
prefix=/illumina_ST-J0010_flowcellA_SampleId

cd $processed_direc

 # for file in `ls $read_direc`;
 # do
 #    foo=${file#$prefix};
 #    foo=${foo%$suffix};
 #    sample=$foo;
 #    export sample;    #exports the values of $sample so that they can be used as arguments in the script


# Here is where the mapping script starts.
# To run it, define the read_direc, processed_direc, suffix and prefix,
# and then copy the for loop above and uncomment it, and run it in the terminal.

# Print current sample ID
echo $sample

# Super-Deduper to remove PCR duplicates
echo "Executing super_deduper"
/ebio/abt6_projects9/Pseudomonas_diversity/Programs/Super-Deduper/super_deduper -1 $read_direc/*$sample*/*R1* -2 $read_direc/*$sample*/*R2* -p $sample.dedup -g
rm $sample.dedup_nodup_PE2.fastq.gz

# cutadapt to remove pIT2 sequence from the reads
# Discard sequences that are not untrimmed (do not start by the transposon sequence, i.e. do not correspond to transposon insertions)
# Trimmed sequences are stores in a new file = sample.R1_trimmed
echo "Executing cutadapt"
cutadapt -g ^CAGGACGCTACTTGTGTATAAG --discard-untrimmed $sample.dedup_nodup_PE1.fastq.gz  | grep -A2 -B1 ^AGTCA  | sed '/^--$/d' > temp.$sample.fastq
cutadapt -g ^AGTCA --discard-untrimmed  temp.$sample.fastq > $sample.R1_trimmed.fastq

# Then, map to P25.C2 genome
# bwa maps reads against the reference genome
# but first the ref genome has to be indexed (this has to be done only once, then
# it can be commented out, and the index is defined as below)
# bwa index /ebio/abt6_projects8/Pseudomonas_mapping/data/ale_tn_mapping/p25.c2_genome/plate25.C2.pilon.contigs_renamed.fasta
#index=/ebio/abt6_projects8/Pseudomonas_mapping/data/ale_tn_mapping/p25.c2_genome/plate25.C2.pilon.contigs_renamed.fasta
index=/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon.fasta

# map reads to index, generates a .sam file as output
echo "Executing mapping"
bwa mem -t 4 $index  $sample.R1_trimmed.fastq > $sample.sam

# samtools
samtools view -bT $index $sample.sam >$sample.bam

# to get the stats for mapping
# generate BAM files
echo "Generating stats"
echo "flagstat"
samtools flagstat $sample.bam
echo "sort"
samtools sort $sample.bam -o $sample.bam
echo "index"
samtools index $sample.bam $sample.bai

# #output the start and stop positions of each read that maps
# echo "Generating tab file"
# samtools view $sample.bam|awk '{print $3 "\t" $4 "\t" $4+length($10)-1}' > $sample.tab

# #output bed file and overlap with P25.C2
# p25c2=/ebio/abt6_projects8/Pseudomonas_mapping/data/ale_tn_mapping/p25.c2_RAST/p25c2_rast.gff

# # This bed file can be loaded into the IGV as a track of P25.C2
# echo "Generating bed file for IGV"
# samtools view $sample.bam |awk '{print $3 "\t" $4 "\t" $4+length($10)-1}' > $processed_direc/$sample.bed

# # The intersected bed file can be used for further analysis of insertions
# echo "Generating intersectBed"
# intersectBed -a $p25c2 -b $sample.bam > $processed_direc/$sample.p25c2.bed
# done
