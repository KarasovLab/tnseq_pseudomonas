#!/bin/sh
# 
#  Request 20G of RAM
#$ -l h_vmem=20G
#$ -o $HOME/tmp/stdout_of_job
#$ -j y
#$ -S /bin/bash

#need to remove multiple sequences from the samples
#https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapter
#P5_Nex_Tn_amend_S502	
#AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCAGGACGCTACTTGTGTATAAGA
#pIT2
#CAGGACGCTACTTGTGTATAAGA
#adapter 1
#CTGTCTCTTATACACATCTGACGCTGCCGACGA

#for file in `ls $read_direc`; do foo=${file#$prefix}; foo=${foo%$suffix}; sample=$foo; qsub  -v read_direc=$read_direc,sample=$sample /ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/code/s1_general_adapter_trimming.sh; done

read_direc=/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/raw_reads/hiseq_0141 #$read_direc
sample=$sample
processed_direc=/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/processed_reads/hiseq_0141 #$processed_direc
suffix=_RunId0091_LaneId1 #$suffix #suffix=_LaneId1
prefix=illumina_ST-J00101_flowcellA_SampleId

cd $processed_direc

/ebio/abt6_projects9/Pseudomonas_diversity/Programs/Super-Deduper/super_deduper -1 $read_direc/*$sample*/*R1* -2 $read_direc/*$sample*/*R2* -p $sample.dedup -g
rm $sample.dedup_nodup_PE2.fastq.gz
/ebio/abt6/tkarasov/.local/bin/cutadapt -g ^CAGGACGCTACTTGTGTATAAG --discard-untrimmed $sample.dedup_nodup_PE1.fastq.gz  | grep -A2 -B1 ^AGTCA  | sed '/^--$/d' > temp.$sample.fastq
/ebio/abt6/tkarasov/.local/bin/cutadapt -g ^AGTCA --discard-untrimmed  temp.$sample.fastq > $sample.R1_trimmed.fastq
                                           
#now to map to the NP29.1a genome which has already been indexed with following command
#bwa index /ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/NP29_index/NP29.1a_11252011.fasta
index=/ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/NP29_index/NP29.1a_11252011.fasta

bwa mem -t 4 $index  $sample.R1_trimmed.fastq > $sample.sam
samtools view -bT $index $sample.sam >$sample.bam

#to get the stats for mapping
samtools flagstat $sample.bam
samtools sort $sample.bam $sample
samtools index $sample.bam $sample.bai

#output the start and stop positions of each read that maps
samtools view $sample.bam|awk '{print $3 "\t" $4 "\t" $4+length($10)-1}' > $sample.tab

#output bed file and overlap with NP29
NP29=/ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/NP29_index/317.111_no_dup.gff
samtools view $sample.bam |awk '{print $3 "\t" $4 "\t" $4}' > $sample.bed
intersectBed -a $NP29 -b $sample.bed  > $sample.NP29.bed; done


