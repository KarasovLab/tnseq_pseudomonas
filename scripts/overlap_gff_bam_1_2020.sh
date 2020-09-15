#edited now for
cd /ebio/abt6_projects8/Pseudomonas_mapping/poo
for sample in `ls /ebio/abt6_projects8/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/ | grep bam`;
do

#need to sort and index the bam files
samtools sort /ebio/abt6_projects8/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/$sample > sorted.$sample
samtools index sorted.$sample

#and call bedtools (on all bams)
bedtools multicov -bams ./*.bam -bed plate25.C2.annotation_no_fasta.gff > counts.gff
bedtools multicov -bams ./*.bam -bed /ebio/abt6_projects8/Pseudomonas_mapping/poo/plate25.c2_fin.gff > /ebio/abt6_projects8/Pseudomonas_mapping/poo/counts2.gff

head gene_counts.gff
sed 's/^.*locus_tag=//' counts2.gff > gene_counts.tab





#with the original gff  the chromosome name is wrong. We need to relabel the gff file
# sed -e 's/gnl|Prokka|DAKCFMEO_//' plate25.C2.annotation_no_fasta.gff > plate25.c2_fin.gff

bedtools multicov -bams ./*.bam -bed /ebio/abt6_projects8/Pseudomonas_mapping/poo/plate25.c2_fin.gff > counts2.gff