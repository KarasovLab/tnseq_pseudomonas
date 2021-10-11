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

#for file in `ls $read_direc`; do foo=${file#$prefix}; foo=${foo%$suffix}; sample=$foo; qsub  -v sample=$sample /ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/code/adapter_trimming.sh; done

cd /ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/hiseq0081
read_direc=/ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/raw_reads/hiseq_run0081
prefix=illumina_ST-J00101_flowcellA_SampleId
suffix=_LaneId2

sample=$sample
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










#zcat $read_direc/*$sample*/*R1* | awk '{print (NR%4 == 1) ? "@NP29_2_" ++i : $0}'|grep TATAAGAGTCA -A2 -B1 | sed '/^--$/d'  > $sample.R1.fastq
#zcat $read_direc/*$sample*/*R1* | grep TATAAGAGTCA -A2 -B1 | sed '/^--$/d'  > $sample.R1.fastq
#zcat $read_direc/*$sample*/*R1*  > $sample.R1.fastq
#agrep -2   CAGGACGCTACTTGTGTATAAGAGTCA   $sample.R1.fastq >  $sample.R1.txt
#grep -f $sample.R1.txt  -A2 -B1  $sample.R1.fastq > $sample.R1.tn.fastq
#grep -f $sample.R1.txt  -A2 -B1  $sample.R1.fastq | sed '/^--$/d'   $sample.R1.fastq.gz> $sample.R1.fastq
#rm -r  $read_direc/*$sample*/
#cat $read_direc/$read1 | grep CAGGACGCTACTTGTGTATAAGAGTCA -A2 -B1 | sed '/^--$/d' >  $read1.tn_only.fastq
#trying to use fast clipper to remove tn insertion site but this does not work well
#fastx_clipper -a CAGGACGCTACTTGTGTATAAGAGTCA -c -i $sample.R1.fastq -Q33 > $sample.R1_trimmed.fastq
#if using second read (we are not)
#replace="2:N"
#cat $sample.R1.fastq | grep "^@" | sed -e "s/1:N/$replace/g"> $sample.temp
#zcat $read_direc/*$sample*/*R2* | grep -Ff $sample.temp  -A3 > $sample.R2.fastq
#select for only reads that start with transposons
#cat $read_direc/$read1 | grep CAGGACGCTACTTGTGTATAAGAGTCA -A2 -B1 | sed '/^--$/d' >  $read1.tn_only.fastq
#cat $read1.tn_only.fastq | grep "^@"| cut -f 1 -d ' ' > temp
#grep -Ff temp $read2 -A3 > $read2.tn_only.fastq



#samtools rmdup -S $sample.bam $sample.no_dup.bam
#samtools sort $sample.no_dup.bam $sample.no_dup
#samtools index $sample.no_dup.bam $sample.no_dup.bai
#samtools mpileup $sample.tn_removed_sort.bam > $sample.pileup
#samtools rmdup -s both_reads_tn_removed_sort.bam both_reads_tn_removed_sort_no_dup.bam #but this removes poorly
#bedtools merge -i both_reads_tn_removed_sort_no_dup.bam -o both_reads_tn_removed_sort_no_dup_merge.bam
#bedtools merge -i both_reads_tn_removed_sort.bam -o both_reads_tn_removed_sort_merge.bam
#bedtools coverage -abam both_reads_tn_removed_sort_no_dup.bam -b cnv_regions.bed
#samtools index both_reads_tn_removed_sort_no_dup.bam both_reads_tn_removed_sort_no_dup.bai
#convert to bed and extract start and stop positions for each read then limit to only those reads with uniqueness
#bamToBed -i both_reads_tn_removed_sort_no_dup.bam > both_reads_tn_removed_sort_no_dup.bed
#cat  both_reads_tn_removed_sort_no_dup.bed | awk '$5 == "37"' > both_reads_tn_removed_sort_no_dup_unique.bed
#gff2bed < 317.111.gff > 317.111.gff.bed

