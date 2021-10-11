library("GenomicFeatures") #this provides tools for extracting genomic features based on locations after specifying which genome browser to use 
#library("DESeq2")
#BiocManager::install("DESeq2")
library( "Rsamtools" ) # samtools interface for R
library("GenomicAlignments") # containers for storing gene alignments - read counting, computing coverage, nuc content etc.
library('tidyr') # package for 'tidying' data with limited set of functions for reshaping
library('dplyr') # set of functions for data manipulation
library(edgeR) # differential expression analysis
library(ggplot2) # makes fancy plots

#https://wikis.utexas.edu/display/bioiteam/Differential+gene+expression+analysisetp
#and more https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

#sample_info = read.csv("/ebio/abt6_projects8/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/sampleInfo_TnSeq_Oct2019.csv", header = T, row.names = 1)

#the order in the counts2.gff is the order in *.bam (assigns 'sample_order' to an order table from a bam_order text file that can be used to order the count files?)
# should there be a separate bam_order file? Is it the same for all count files or do we need to write this in a way that is inclusive of all the count files?
sample_order = read.table("/ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process/bam_order.txt")


# this reads in the count files from a delimited file, specifies that the first line is not column headers and assigns the data to the variable 'counts'
counts = read.delim("/ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process/counts_plaurin.gff", header=F)

# count_table is a subset of counts including all rows and all the columns from 10 on? (I'm understanding dim(counts)[2] giving the columns count from the table dimensions) Is this right? Why are the first 10 excluded?

count_table = counts[,c(10:dim(counts)[2])]
rownames(count_table) = counts$V9
colnames(count_table) = sample_order[,1]
head(counts)

#Setting up the factors
group_type <- factor(separate(sample_order, V1, into = c("plate", "pos", "time"), sep ="_"  )$time)
group_batch <- as.character(group_type)
group_batch[which(group_batch == "33")] <- "Treated"
group_batch[which(group_batch == "310")] <- "Treated"
group_batch[which(group_batch == "33t0")] <- "T0"
group_batch[which(group_batch == "310t0")] <- "T0"
group_type[which(group_type == "33t0")] <- "33"
group_type[which(group_type == "310t0")] <- "310"

group_batch = as.factor(group_batch)
samp_info <- data.frame(batch = group_type, class = group_batch) %>% filter(class %in% c("T0", "Treated"))
samp_info$class <- as.factor(as.character((samp_info$class)))
samp_info$batch <- as.factor(as.character((samp_info$batch)))
count_table = count_table[,which(group_batch %in% c("T0", "Treated"))]

dge = DGEList(counts=count_table)
# countsPerMillion <- cpm(dge)
# summary(countsPerMillion)
# countCheck <- countsPerMillion > 1
# keep <- which(rowSums(countCheck) >= 2)
dge.trimmed <- dge #[keep,]

#Normalize library size
dge <- calcNormFactors(dge.trimmed, method="none") # Normalize library sizes using TMM

#Perfomr MDS on dge.trimmed
mds <- t(dge.trimmed$counts) %>% dist() %>% cmdscale() %>% data.frame()
colnames(mds) <- c("MDS1", "MDS2")

pdf("/ebio/abt6_projects8/Pseudomonas_mapping/data/fig_misc/MDS_plaurin.pdf", useDingbats = FALSE, fonts = "ArialMT")

ggplot(data = mds, aes(x = MDS1, y = MDS2)) +
  geom_point( aes(color = samp_info$class, shape = samp_info$batch), cex = 3) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

###########################
# THIS IS WHERE THE REGRESSION COEFFICIENT ESIMATION BEGINS
###########################

# Make the model matrix
design = model.matrix(~samp_info$class + samp_info$batch) # Create design matrix for glm

# Estimate the GLM
dge = estimateGLMCommonDisp(dge, design)
dge = estimateGLMTagwiseDisp(dge, design)
fit = glmFit(dge,design)

#lrt2vs1 <- glmLRT(fit, coef=2)
#lrt3vs1 <- glmLRT(fit, coef=3)

# Perform the LRT
???LINES MISSING
# sampleInfo_sep2 = separate(real_thing, 1, into = c("sample_name","position", "treatment"), sep = "_")
# head(sampleInfo_sep2)
# sampleInfo_sep2 = separate(sampleInfo_sep2, treatment, into = c("batch", "timepoint"), sep = "_")
# 
# # data resulting as factor not character
# 
# # create index
# hm <- separate(sampleInfo2, x, into = c("x1","x2", "x3", "x4", "X5", "X6"), sep = "/")$X6
# seIdx <- match(colnames(se), hm)
# str(seIdx)
# colData(se) <- cbind( colData(se), sampleInfo_sep2[ seIdx, ] )
# ddsFull <- DESeqDataSet( se, design = ~ batch + treatment )
# > se
# 
# 
# #playing
# #sample_all <- data.frame(sampleInfo2, sampleInfo_sep)
# #cbind (sampleInfo_sep,sampleInfo2)
# #clean data:
# #remove: empty, negcon, control needed: > '%!in%' <- function(x,y)!('%in%'(x,y))
# '%!in%' <- function(x,y)!('%in%'(x,y))
# #new = sampleInfo_sep[sampleInfo_sep$treatment %!in% c("negcon", "control", "water", "empty"),]
# #new
# 
# 
# #rename treatment from 33 and 310 > 33t2 and 310t2 #doesn't work yet: "unexpected symbol in "temp_[temp_=="33"] = 33t2"":
# temp_ =sampleInfo_sep2$treatment
# temp_[which(temp_=="33.bam")] = "33_t2"
# temp_[which(temp_=="310.bam")] = "310_t2"
# temp_[which(temp_=="310t0.bam")] = "310_t0"
# temp_[which(temp_=="33t0.bam")] = "33_t0"
# temp_
# sampleInfo_sep2$treatment=temp_
# sampleInfo_sep2
# 
# #write.csv(sampleInfo_sep, file = "sampleInfo_sep.csv"
# # csv file ends up in "home" here: mhoelscher 
