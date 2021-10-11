w#merge bam files
cd /ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/hiseq0081
#for file in `ls | grep day8 | grep bam | grep -v all`; do samtools view $file  |awk '{print $3 "\t" $4 "\t" $4}' > $file.bed; done
#for rec in `ls | grep bed`; do NP29=/ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/NP29_index/317.111_no_dup.gff; intersectBed -a $NP29 -b $rec  > $rec.NP29.bed; done

samtools merge all_Canola_8d.bam Canola_8d*.bam
samtools merge all_NP29.bam NP29*.bam
samtools merge all_Peas_8d.bam Peas_8d*.bam
samtools merge all_Barley_8d.bam Barley_8d*.bam
samtools merge all_Col-0_8d.bam Col-0_8d*.bam


#index bam files
samtools index all_Barley_8d.bam
samtools index all_Col-0_8d.bam
samtools index all_Peas_8d.bam
samtools index all_Canola_8d.bam
samtools index all_NP29.bam

#output start positions
#output the start and stop positions of each read that maps
samtools view all_Barley_8d.bam |awk '{print $3 "\t" $4 "\t" $4}' > all_Barley_8d.bed
# $4+length($10)-1}' > all_Barley_8d.tab
samtools view all_Canola_8d.bam |awk '{print $3 "\t" $4 "\t" $4}' > all_Canola_8d.bed
samtools view all_Col-0_8d.bam |awk '{print $3 "\t" $4 "\t" $4}' > all_Col-0_8d.bed
samtools view all_NP29.bam |awk '{print $3 "\t" $4 "\t" $4}' > all_NP29.bed
samtools view all_Peas_8d.bam |awk '{print $3 "\t" $4 "\t" $4}' > all_Peas_8d.bed

NP29=/ebio/abt6_projects9/Pseudomonas_diversity/Tnseq/processed_reads/NP29_index/317.111_no_dup.gff
#now i want to intersect bed with gff
array=( "all_Barley_8d.bed" "all_Canola_8d.bed" "all_NP29.bed" "all_Peas_8d.bed" "all_Col-0_8d.bed" "317.111_no_dup.gff" )
for rec in "${array[@]}"; do
intersectBed -a $NP29 -b $rec  > $rec.NP29.bed; done
