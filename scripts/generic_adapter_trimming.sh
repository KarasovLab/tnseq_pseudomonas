#!/bin/bash

#SBATCH -t 14-00:00:00 # assign time limit

#SBATCH -J assemble # assign job name

#SBATCH -N 1 # node number

#SBATCH -o out%tnseq # assign standard output file name

#SBATCH -e err%tnseq # assign standard error file name

#SBATCH -n 8 # assign CPU number (the owner cluster should have 128 threads for 68 cores)

#SBATCH -A mpag-np  # group name (the owner node name)

#SBATCH -p mpag-shared-np  # microbe, plant and animal genetics (the partition name). Note: mpag-np (to have the entire node) or mpag-shared-np (to request a portion of the node)

#SBATCH --mail-type=begin,end # email you when begin, end, fail, requeue

#SBATCH --time=01:00:00

## SBATCH --mail-user=t.karasov@utah.edu # your email address

## SBATCH --mem=500000 # request memory (the owner node have maximum 500 GB memory, assign it if you know the exact memory you want to use)

## SBATCH -A karasov-group1 # my group name

## SBATCH --partition notchpeak-freecycle # request notchpeak-freecycle for script testing



sleep 15

##need to remove multiple sequences from the samples
##https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapter
##P5_Nex_Tn_amend_S502	
##AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCAGGACGCTACTTGTGTATAAGA
##pIT2
##CAGGACGCTACTTGTGTATAAGA
##adapter 1
##CTGTCTCTTATACACATCTGACGCTGCCGACGA
#
###for file in `ls $read_direc`; do foo=${file#$prefix}; foo=${foo%$suffix}; sample=$foo; qsub  -v sample=$sample /ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/code/adapter_trimming.sh; done
#
#for file in `ls $read_direc`; do foo=${file#$prefix}; foo=${foo%$suffix}; sample=$foo; qsub  -v sample=$sample /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/scripts/tnseq_pseudomonas/scripts/adapter_trimming.sh; done
#
#
cd /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/Talia/
#read_direc=$read_direc
full_name=$read_direc/*R1*
#sample=$new_name
#index=$index
## sample is read1 from illumina_ST.....
#
##/ebio/abt6_projects9/Pseudomonas_diversity/Programs/Super-Deduper/super_deduper -1 $read_direc/*R1* -2 $read_direc/*R2* -p $sample.dedup -g
#
/uufs/chpc.utah.edu/common/home/u0572090/anaconda3/bin/hts_SuperDeduper -1 $read_direc/*R1* -2 $read_direc/*R2* -f $sample.dedup -F
#
##rm $sample.dedup_nodup_PE2.fastq.gz
#
##/ebio/abt6/tkarasov/.local/bin/cutadapt -g ^CAGGACGCTACTTGTGTATAAG --discard-untrimmed $sample.dedup_nodup_PE1.fastq.gz  | grep -A2 -B1 ^AGTCA  | sed '/^--$/d' > temp.$sample.fastq
##/ebio/abt6/tkarasov/.local/bin/cutadapt -g ^AGTCA --discard-untrimmed  temp.$sample.fastq > $sample.R1_trimmed.fastq
/uufs/chpc.utah.edu/common/home/karasov-group1/bin/cutadapt -g ^CAGGACGCTACTTGTGTATAAG --discard-untrimmed $sample.dedup_R1.fastq.gz  | sed '/^--$/d' > temp.$sample.fastq
/uufs/chpc.utah.edu/common/home/karasov-group1/bin/cutadapt -g ^AGTCA --discard-untrimmed  temp.$sample.fastq > $sample.R1_trimmed.fastq
#                                          
##now to map to the NP29.1a genome which has already been indexed with following command
#bwa index  /uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/DC3000_GCF_000007805.1/GCF_000007805.1_ASM780v1_genomic.fna
##index=/ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/NP29_index/NP29.1a_11252011.fasta
#index=/uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/p7_g9/plate7.G9.pilon.contigs_renamed.fasta
index=/uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/p25_c2/plate25.C2.pilon.contigs_renamed.fasta

if [[ "$index" == "DC3000" ]]; then
index=/uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/DC3000_GCF_000007805.1/GCF_000007805.1_ASM780v1_genomic.fna; fi

if [[ "$index" == "p7_g9" ]]; then
index=/uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/p7_g9/plate7.G9.pilon.contigs_renamed.fasta; fi

if [[ "$index" == "p25_c2" ]]; then
index=/uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/p25_c2/plate25.C2.pilon.contigs_renamed.fasta; fi
# index=/uufs/chpc.utah.edu/common/home/karasov-group1/tnseq/data/ref_database/references_sequences/p25_c2/p25.C2.contigs.second_polished.pilon.fasta ; fi
## SBATCH -p notchpeak-freecycle # request notchpeak-freecycle for script testing





#
/uufs/chpc.utah.edu/common/home/u0572090/anaconda3/bin/bwa mem -t 4 $index  $sample.R1_trimmed.fastq > $sample.sam
samtools view -bT $index $sample.sam >$sample.bam
#
##to get the stats for mapping
samtools flagstat $sample.bam
samtools sort -o $sample.bam $sample.bam
samtools index $sample.bam $sample.bai
#
##output the start and stop positions of each read that maps
samtools view $sample.bam|awk '{print $3 "\t" $4 "\t" $4+length($10)-1}' > $sample.tab
#









