#edited now for p laurin's data this script takes the bam files and overlaps them with gene annotations. This can then go into edgeR
cd /ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process
for sample in `ls /ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process/ | grep bam | grep -v sorted`;
	do
	#need to sort and index the bam files
		samtools sort /ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process/$sample > sorted.$sample
		samtools index sorted.$sample
done


#with the original gff  the chromosome name is wrong. We need to relabel the gff file
sed -e 's/utg000001l_p25.C2/utg000001c:1.0-5963307.0_pilon/' \
/ebio/abt6_projects8/Pseudomonas_mapping/data/mapping/references_sequences/p25.C2.contigs.second_polished.pilon_no_fasta.gff > plate25.c2_fin.gff


#and call bedtools (on all bams) with overlap with full genome
bedtools multicov -bams ./sorted.*.bam -bed plate25.c2_fin.gff  > \
/ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process/counts_plaurin.gff

head gene_counts.gff
sed 's/^.*locus_tag=//' counts_plaurin.gff > gene_counts.tab





