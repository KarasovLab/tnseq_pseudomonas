library("GenomicFeatures")
#library("DESeq2")
#BiocManager::install("DESeq2")
library( "Rsamtools" )
library("GenomicAlignments")
library('tidyr')
library('dplyr')
library(edgeR)
library(ggplot2)

#https://wikis.utexas.edu/display/bioiteam/Differential+gene+expression+analysisetp
#and more https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

#sample_info = read.csv("/ebio/abt6_projects8/Pseudomonas_mapping/mneumann_tnseq/TnSeq_Oct19/sampleInfo_TnSeq_Oct2019.csv", header = T, row.names = 1)

#the order in the counts2.gff is the order in *.bam
sample_order = read.table("/ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process/bam_order.txt")

counts = read.delim("/ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process/counts_plaurin.gff", header=F)

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


# Make the model matrix
design = model.matrix(~samp_info$class + samp_info$batch) # Create design matrix for glm

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


pdf("/ebio/abt6_projects8/Pseudomonas_mapping/data/fig_misc/tnseq_diff_plaurin_full_genome.pdf", useDingbats = FALSE, fonts = "ArialMT")

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


write.csv(etp$table, "/ebio/abt6_projects8/Pseudomonas_mapping/data/tnseq/plaurin_talia_process/edgeR-wt-vs-mut.csv")

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