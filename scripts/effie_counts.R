library("GenomicFeatures") #this provides tools for extracting genomic features based on locations after specifying which genome browser to use 
#library("DESeq2")
#BiocManager::install("DESeq2")
library( "Rsamtools" ) # samtools interface for R
library("GenomicAlignments") # containers for storing gene alignments - read counting, computing coverage, nuc content etc.
library(tidyr) # package for 'tidying' data with limited set of functions for reshaping
library('dplyr') # set of functions for data manipulation
library(edgeR) # differential expression analysis
library(ggplot2) # makes fancy plots

setwd("/Users/talia/Google Drive/My Drive/Utah_Professorship/projects/Tnseq/scripts/tnseq_pseudomonas/data")
#https://wikis.utexas.edu/display/bioiteam/Differential+gene+expression+analysisetp
#and more https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

#sample_info = read.csv("/ebio/abt6_projects8/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/sampleInfo_TnSeq_Oct2019.csv", header = T, row.names = 1)

#the order in the counts2.gff is the order in *.bam (make an order table from a bam_order text file that can be used to order the count files?)
# should there be a separate bam_order file? Is it the same for all count files or do we need to write this in a way that is inclusive of all the count files?
sample_order = read.table("./sample_info_210726_10_11_2021.txt", header = TRUE)
rownames(sample_order) <- sample_order$file_name


# makes a new table, "counts", by reading in the delimited file (the gffs) - the first row does not contain column headers
counts = read.delim("./corrected_effie_count_10_11_2021.txt", header = TRUE)
counts = counts %>% filter(counts$V2 != "Protein Homology")
# count_table takes all rows and columns 10 on from 'counts' (I understand dim(counts)[2] to give the number of columns in 'counts'). Why is this starting from 10? <- 10 is where the counts start. Before that is information about the region being counted. Is each column after for one of the samples/bams?
count_list <- names(which(colSums(counts[,c(10:dim(counts)[2])])>20000))
count_table = counts[,count_list]

#reorder sample_order
sample_order = sample_order[count_list,]
rownames(count_table) = counts$gene_id



# dge will be an object with rows that specify features, columns that specify samples and each cell will be the counts for the feature for each sample
dge = DGEList(counts=count_table, samples = sample_order)
# countsPerMillion <- cpm(dge)
# summary(countsPerMillion)
# countCheck <- countsPerMillion > 1
# keep <- which(rowSums(countCheck) >= 2)

#not sure what is being trimmed or how
dge.trimmed <- dge 

#Normalize library size
dge <- calcNormFactors(dge.trimmed, method="none") # Normalize library sizes using TMM

#Perform MDS on dge.trimmed
# is this the volcano plot?
make_mds<- function(dge.thing) {
  mds <- sqrt(t(dge.thing$counts)) %>% dist() %>% cmdscale(k=10) %>% data.frame()
 mds_eig <- sqrt(t(dge.thing$counts)) %>% dist() %>% cmdscale(eig=TRUE, k=10)
 colnames(mds) <- c("MDS1", "MDS2", "MDS3", "MDS4", "MDS5", "MDS6", "MDS7", "MDS8", "MDS9", "MDS10")
 round(mds_eig$eig*100/sum(mds_eig$eig),1)
 colnames(mds_eig$points) <- c("MDS1", "MDS2", "MDS3", "MDS4", "MDS5", "MDS6", "MDS7", "MDS8", "MDS9", "MDS10")
  return(mds_eig)
}


is_5 <- dge.trimmed[,which(dge.trimmed$samples$concentration==5)]
is_15 <- dge.trimmed[,which(dge.trimmed$samples$concentration==15)]
mds <- make_mds(dge.trimmed)
mds_15 <- make_mds(is_15)
mds_5 <- make_mds(is_5)


pdf("../figs/MDS_all_sqrt.pdf", useDingbats = FALSE, fonts = "ArialMT")

p1 <- ggplot(data = data.frame(mds$points), aes(x = MDS1, y = MDS2)) +
  geom_point( aes(color = dge.trimmed$samples$plant_extract, shape = as.factor(sample_order$concentration)), cex = 3) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p2 <- ggplot(data = data.frame(mds_15$points), aes(x = MDS1, y = MDS2)) +
   geom_point( aes(color = is_15$samples$plant_extract), cex = 3) +
   scale_color_viridis_d() +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p3 <- ggplot(data = data.frame(mds_5$points), aes(x = MDS1, y = MDS2)) +
  geom_point( aes(color = is_5$samples$plant_extract), cex = 3) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot_grid(p1+theme(legend.title = element_blank()) ,p2+theme(legend.title = element_blank()),p3+theme(legend.title = element_blank()), ncol =1)

dev.off()


###########################
# THIS IS WHERE THE REGRESSION COEFFICIENT ESIMATION BEGINS
###########################

# Make the model matrix
design = model.matrix( ~ plant_extract + concentration,data = dge.trimmed$samples) # Create design matrix for glm

# Estimate the GLM
dge = estimateGLMCommonDisp(dge, design)
dge = estimateGLMTagwiseDisp(dge, design)
fit = glmFit(dge,design)

lrt2vs1 <- glmLRT(fit, coef=)
#lrt3vs1 <- glmLRT(fit, coef=3)

# Perform the LRT
lrt3vs2 <- glmLRT(fit, coef = 3)
etp = topTags(lrt3vs2, n=100000)
deGenes <- decideTestsDGE(lrt3vs2, p=0.001)

# dge <- estimateCommonDisp(dgList)
# dge <- estimateTagwiseDisp(dge)
# et <- exactTest(dge, pair = c("33t0", "33")) # This output the comparsion of B - A so genes with positive log-fold change are uregulated in group B
# etp <- topTags(et, n=100000)
etp$table$logFC = -etp$table$logFC


fin_data = etp$table


pdf("../figs/tnseq_diff_p25c2_full_genome.pdf", useDingbats = FALSE, fonts = "ArialMT")

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