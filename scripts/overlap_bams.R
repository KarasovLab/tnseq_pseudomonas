library("GenomicFeatures")
library("DESeq2")
#BiocManager::install("DESeq2")
library( "Rsamtools" )
library("GenomicAlignments")
library('tidyr')
library('dplyr')

#Tutorial on how to use deseq
#before you start connect to server smb://reo.eb.local/abt6_projects8/Pseudomonas_mapping
#Load the gene annotation file
hse <- makeTxDbFromGFF( "/Volumes/Pseudomonas_mapping/data/ale_tn_mapping/p25.c2_genome/annotation/plate25.C2.annotation.gff", format="auto" )
exonsByGene <- exonsBy( hse, by="gene" )

#Load the bam files
fls <- list.files( "/Volumes/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19", pattern="bam$", full=TRUE )



#Specify lines to read in at a time (Not clear why but we are doing it anyways)
bamLst <- BamFileList( fls, yieldSize=100000 )

#Run overlap analysis
se <- summarizeOverlaps( exonsByGene, bamLst,
                         mode="Union",
                         singleEnd=FALSE,
                         ignore.strand=TRUE,
                         fragments=TRUE )
#! needs time to run! check data:colnames 
se
#class: RangedSummarizedExperiment 
#dim: 6 162 
#metadata(0):
# assays(1): counts
#rownames(6): Hgd_1 Hgd_2 ... aaeB_1 aaeB_2
#rowData names(0):
#  colnames(162): plate1_A1_33.bam plate1_A10_310.bam ... plate2_H8_310t0.bam plate2_H9_310t0.bam
#colData names(0):
head( assay(se) )
colSums( assay(se) )
colData(se)
rowData(se)


#using the data for analysis

# need fls (all bam files) as csv file:
write.csv(fls, file.path("/Volumes/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19", "bamfiles.csv"), row.names=FALSE )

# Read in sample info file
# for separate function you need "library('tidyr')"
# sampleInfo <- read.csv("/path/to/file.CSV" )

#sampleInfo <- read.csv( "/Volumes/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/name2.csv" )
#sampleInfo_sep = separate(sampleInfo, x, into = c("sample_name", "treatment"), sep = "-")
#sampleInfo_sep = separate(sampleInfo_sep, treatment, into = c("treatment", "timepoint"), sep = "       ")
#sampleInfo_sep$sample_name = rename(sampleInfo_sep$sample_name, /Volumes/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/ into 


sampleInfo2 <- read.csv( "/Volumes/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/bamfiles.csv")
sampleInfo_sep2 = separate(sampleInfo2, x, into = c("x1","x2", "x3", "x4", "X5", "X6"), sep = "/")
head(sampleInfo_sep2)
real_thing = sampleInfo_sep2$X6
real_thing = as.data.frame(real_thing)
sampleInfo_sep2 = separate(real_thing, 1, into = c("sample_name","position", "treatment"), sep = "_")
head(sampleInfo_sep2)
sampleInfo_sep2 = separate(sampleInfo_sep2, treatment, into = c("batch", "timepoint"), sep = "_")

# data resulting as factor not character

# create index
hm <- separate(sampleInfo2, x, into = c("x1","x2", "x3", "x4", "X5", "X6"), sep = "/")$X6
seIdx <- match(colnames(se), hm)
str(seIdx)
colData(se) <- cbind( colData(se), sampleInfo_sep2[ seIdx, ] )
ddsFull <- DESeqDataSet( se, design = ~ batch + treatment )
> se


#playing
#sample_all <- data.frame(sampleInfo2, sampleInfo_sep)
#cbind (sampleInfo_sep,sampleInfo2)
#clean data:
#remove: empty, negcon, control needed: > '%!in%' <- function(x,y)!('%in%'(x,y))
'%!in%' <- function(x,y)!('%in%'(x,y))
#new = sampleInfo_sep[sampleInfo_sep$treatment %!in% c("negcon", "control", "water", "empty"),]
#new


#rename treatment from 33 and 310 > 33t2 and 310t2 #doesn't work yet: "unexpected symbol in "temp_[temp_=="33"] = 33t2"":
temp_ =sampleInfo_sep2$treatment
temp_[which(temp_=="33.bam")] = "33_t2"
temp_[which(temp_=="310.bam")] = "310_t2"
temp_[which(temp_=="310t0.bam")] = "310_t0"
temp_[which(temp_=="33t0.bam")] = "33_t0"
temp_
sampleInfo_sep2$treatment=temp_
sampleInfo_sep2

#write.csv(sampleInfo_sep, file = "sampleInfo_sep.csv"
# csv file ends up in "home" here: mhoelscher 