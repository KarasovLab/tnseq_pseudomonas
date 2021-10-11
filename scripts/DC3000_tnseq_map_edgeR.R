library("GenomicFeatures") #this provides tools for extracting genomic features based on locations after specifying which genome browser to use 
#library("DESeq2")
#BiocManager::install("DESeq2")
library( "Rsamtools" ) # samtools interface for R
library("GenomicAlignments") # containers for storing gene alignments - read counting, computing coverage, nuc content etc.
library(tidyr) # package for 'tidying' data with limited set of functions for reshaping
library('dplyr') # set of functions for data manipulation
library(edgeR) # differential expression analysis
library(ggplot2) # makes fancy plots

#https://wikis.utexas.edu/display/bioiteam/Differential+gene+expression+analysisetp
#and more https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

#sample_info = read.csv("/ebio/abt6_projects8/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/sampleInfo_TnSeq_Oct2019.csv", header = T, row.names = 1)

#the order in the counts2.gff is the order in *.bam (make an order table from a bam_order text file that can be used to order the count files?)
# should there be a separate bam_order file? Is it the same for all count files or do we need to write this in a way that is inclusive of all the count files?
sample_order = read.table("tnseq/DC3000_order_list.txt")

# makes a new table, "counts", by reading in the delimited file (the gffs) - the first row does not contain column headers
counts = read.delim("tnseq/DC3000_counts.gff", header=F)
counts = counts %>% filter(counts$V2 != "Protein Homology")
# count_table takes all rows and columns 10 on from 'counts' (I understand dim(counts)[2] to give the number of columns in 'counts'). Why is this starting from 10? <- 10 is where the counts start. Before that is information about the region being counted. Is each column after for one of the samples/bams?
count_table = counts[,c(10:dim(counts)[2])]
# this uses the list named V9 from 'counts' as the row names - what is V9? Will this be the same for the new gffs?
colnames(count_table) = sample_order[,1]


total_counts <- count_table %>% slice(1)
above = which(total_counts > 1000)
count_table <- count_table[above]
sample_order <- sample_order[above,]
sample_order <- data.frame(sample_order)

# this uses column 1 from the bam_order file as the column names for 'count_table'

rownames(count_table) = counts$V9
head(counts)

#Setting up the factors
#I'm not sure everything this one is doing... It seems like it is separating column V1 from sample_order and storing the new names 'plate', 'pos' and 'time', with columns separated by '_', then it calls the vector you just named 'time', so that factor() is being done on 'time', but I'm not clear what 'factor' does or what is actually stored in 'time'. I'm guessing the time points since 'group-type' is called later and the values 33t0 and 310t0 are replaced? 
group_type <- factor(separate(sample_order, sample_order, into = c("tnseq","strain", "time", "plate","pos"), sep ="_"  )$time)
#Is this taking the factors from group type and converting them to strings? I'm not clear what is meant by 'character'?
group_batch <- as.character(group_type)
# 'which' gives the index of the specified string in the variable you name, so in this case, it finds the strings and uses the index values to replace it with the new string. In group_batch, does the exclusion of t0 refer to the batch at t3? 
group_batch[which(group_batch == "T3")] <- "Treated"
#group_batch[which(group_batch == "310")] <- "Treated"
#group_batch[which(group_batch == "33t0")] <- "T0"
group_batch[which(group_batch == "T0")] <- "T0"
#group_type[which(group_type == "33t0")] <- "33"
#group_type[which(group_type == "310t0")] <- "310"
# If I'm understanding this right, these are all from the list of bams and are based on the file naming conventions used for the earlier experiment. Do we want to use the same conventions here? Will we need to specify strains somewhere as well as time and treated/not?
# Group types will be DC300, p7g9 and p25c2 and each of the controls?
# Group batch will be T0 and T3?

# converts the batch back to a factor
group_batch = as.factor(group_batch)

# makes a data frame with group_type and group_batch and passes it through a filter so that only the samples with the class specified as TO or Treated will be output to samp_info
samp_info <- data.frame(class = group_batch) %>% filter(class %in% c("T0", "Treated"))
# batch and class converted back to characters
samp_info$class <- as.factor(as.character((samp_info$class)))
#samp_info$batch <- as.factor(as.character((samp_info$batch)))
# limits count table to all rows, but only the columns with TO or Treated (based on the index in group_batch)

#count_table = count_table[,which(group_batch %in% c("T0", "Treated"))]

# dge will be an object with rows that specify features, columns that specify samples and each cell will be the counts for the feature for each sample
dge = DGEList(counts=count_table)
# countsPerMillion <- cpm(dge)
# summary(countsPerMillion)
# countCheck <- countsPerMillion > 1
# keep <- which(rowSums(countCheck) >= 2)

#not sure what is being trimmed or how
dge.trimmed <- dge #[keep,]

#Normalize library size
dge <- calcNormFactors(dge.trimmed, method="none") # Normalize library sizes using TMM

#Perform MDS on dge.trimmed
# is this the volcano plot?
mds <- t(dge.trimmed$counts) %>% dist() %>% cmdscale() %>% data.frame()
colnames(mds) <- c("MDS1", "MDS2")

pdf("tnseq/MDS_DC3000test.pdf", useDingbats = FALSE, fonts = "ArialMT")

ggplot(data = mds, aes(x = MDS1, y = MDS2)) +
  geom_point( aes(color = samp_info$class), cex = 3) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

###########################
# THIS IS WHERE THE REGRESSION COEFFICIENT ESIMATION BEGINS
###########################

# Make the model matrix
design = model.matrix(~samp_info$class) # Create design matrix for glm

# Estimate the GLM
dge = estimateGLMCommonDisp(dge, design)
dge = estimateGLMTagwiseDisp(dge, design)
fit = glmFit(dge,design)

#lrt2vs1 <- glmLRT(fit, coef=2)
#lrt3vs1 <- glmLRT(fit, coef=3)

# Perform the LRT
lrt3vs2 <- glmLRT(fit, coef = 2)
etp = topTags(lrt3vs2, n=100000)
deGenes <- decideTestsDGE(lrt3vs2, p=0.001)

# dge <- estimateCommonDisp(dgList)
# dge <- estimateTagwiseDisp(dge)
# et <- exactTest(dge, pair = c("33t0", "33")) # This output the comparsion of B - A so genes with positive log-fold change are uregulated in group B
# etp <- topTags(et, n=100000)
etp$table$logFC = -etp$table$logFC


fin_data = etp$table


pdf("tnseq/tnseq_diff_DC3000_full_genome.pdf", useDingbats = FALSE, fonts = "ArialMT")

ggplot(data = fin_data, aes(x = logFC, y = -log10(FDR))) +
  geom_point(data = subset(fin_data, FDR < 0.01), col = "RED") +
  geom_point(data = subset(fin_data, FDR > 0.01), col = "GREY") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_linedraw() + 
  theme(panel.grid = element_blank()) +
  xlab("Fold change (log2)") +
  ylab("FDR (log10)") 

dev.off()


write.csv(etp$table, "tnseq/edgeR-DC3000.csv")

#underrepresented in plant
sig <- etp[which(etp$table$FDR < 0.01),]
under <- sig[which(sig$table$logFC < 0),]
over <- sig[which(sig$table$logFC > 0),]

# #Tutorial on how to use deseq
# #before you start connect to server smb://reo.eb.local/abt6_projects8/Pseudomonas_mapping
# #Load the gene annotation file
# hse <- makeTxDbFromGFF("~/work_main/abt6_projects8/Pseudomonas_mapping/poo/counts2.gff", format="auto" )
# 
# # exonsbygene
# exonsByGene <- exonsBy( hse, by="gene" )
# 
# #Load the bam files
# fls <- list.files( "~/work_main/abt6_projects8/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19", pattern="bam$", full=TRUE )
# 
# 
# 
# #Specify lines to read in at a time (Not clear why but we are doing it anyways)
# bamLst <- BamFileList( fls, yieldSize=100000 )
# 
# 
# 
# #Run overlap analysis
# se <- summarizeOverlaps( exonsByGene, bamLst,
#                          mode="Union",
#                          singleEnd=FALSE,
#                          ignore.strand=TRUE,
#                          fragments=TRUE )
# #! needs time to run! check data:colnames 
# se
# #class: RangedSummarizedExperiment 
# #dim: 6 162 
# #metadata(0):
# # assays(1): counts
# #rownames(6): Hgd_1 Hgd_2 ... aaeB_1 aaeB_2
# #rowData names(0):
# #  colnames(162): plate1_A1_33.bam plate1_A10_310.bam ... plate2_H8_310t0.bam plate2_H9_310t0.bam
# #colData names(0):
# head( assay(se) )
# colSums( assay(se) )
# colData(se)
# rowData(se)
# 
# 
# #using the data for analysis
# 
# # need fls (all bam files) as csv file:
# write.csv(fls, file.path("/Volumes/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19", "bamfiles.csv"), row.names=FALSE )
# 
# # Read in sample info file
# # for separate function you need "library('tidyr')"
# # sampleInfo <- read.csv("/path/to/file.CSV" )
# 
# #sampleInfo <- read.csv( "/Volumes/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/name2.csv" )
# #sampleInfo_sep = separate(sampleInfo, x, into = c("sample_name", "treatment"), sep = "-")
# #sampleInfo_sep = separate(sampleInfo_sep, treatment, into = c("treatment", "timepoint"), sep = "       ")
# #sampleInfo_sep$sample_name = rename(sampleInfo_sep$sample_name, /Volumes/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/ into 
# 
# 
# sampleInfo2 <- read.csv( "/Volumes/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/bamfiles.csv")
# sampleInfo_sep2 = separate(sampleInfo2, x, into = c("x1","x2", "x3", "x4", "X5", "X6"), sep = "/")
# head(sampleInfo_sep2)
# real_thing = sampleInfo_sep2$X6
# real_thing = as.data.frame(real_thing)
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
