#----------------------------------------------------------------------
# Install packages and load libraries
#----------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2) 
library(ggplot2)

install.packages("tidyverse")
library(tidyverse)

BiocManager::install("phyloseq")
library(phyloseq)

install.packages("devtools")
library(devtools)

install_github("Russel88/MicEco")
library(MicEco)

install.packages("gridExtra")
library(gridExtra)

install.packages("cowplot")
library(cowplot)

install.packages("stringr")
library(stringr)

BiocManager::install("wq")

install.packages("multcompView")
library(multcompView)

#----------------------------------------------------------------------
# Male sexual competitiveness assays
#----------------------------------------------------------------------

# Read in the male competitiveness data
male.comp.data <- data.frame(ID=c(rep(1:6,2),rep(7:12,2),rep(13:18,2)),Group=c(rep("G3 H2O",6),rep("G4 H2O",6),rep("G3 ATC",6),rep("G4 ATC",6),rep("G3 DOX",6),rep("G4 DOX",6)),Prop.mated=c(0.697,0.813,0.774,0.641,0.808,0.583,0.852,0.688,0.759,0.647,0.861,0.75,0.758,0.6,0.816,0.765,0.694,0.788,0.742,0.68,0.769,0.792,0.692,0.71,0.519,0.5,0.5,0.542,0.621,0.581,0.622,0.607,0.618,0.545,0.667,0.593))

# Relevel so that groups are displayed in the desired order on the plot
fct_relevel(factor(male.comp.data$Group),"G3 DOX","G4 DOX","G3 ATC","G4 ATC","G3 H2O","G4 H2O")
male.comp.data$Group <- fct_relevel(factor(male.comp.data$Group),"G3 DOX","G4 DOX","G3 ATC","G4 ATC","G3 H2O","G4 H2O")

# Statistical analysis of results

# First, test for normality using Shapiro-Wilk test
shapiro.test(male.comp.data$Prop.mated) # p-value = 0.3097; we assume data are normally-distributed

# Next, test for homogeneity of variances with Bartlett test
bartlett.test(Prop.mated ~ Group, data=male.comp.data) # p-value = 0.2955

# One-way ANOVA to detect differences in means between all groups
res.aov <- aov(Prop.mated ~ Group, data = male.comp.data)
summary(res.aov)

# Post-hoc testing of pairwise comparisons
# Pairwise t-tests with Benjamini & Hochberg correction for multiple comparisons
pairwise.t.test(male.comp.data$Prop.mated, male.comp.data$Group,
                p.adjust.method = "BH")

# output pairwise t-tests
t.test.output <- pairwise.t.test(male.comp.data$Prop.mated, male.comp.data$Group,
                p.adjust.method = "BH")

# use multcompLetters function to view letters denoting significantly different groups
pvals <- as.data.frame(t.test.output[3][[1]])
pvals <- add_row(pvals,.before = 1)
rownames(pvals)[1] <- "G3 DOX"
pvals <- cbind(pvals,empty_column=NA)
colnames(pvals)[6] <- "G4 H2O"
multcompLetters(Matrix::forceSymmetric(as.matrix(pvals),uplo="L"))

# add a new column to the male competitiveness data with letters
male.comp.data[,4] <- c(rep("b",24),rep("a",12))
colnames(male.comp.data)[4] <- "Letters"

# remove all but one letter from each group (helps with font smoothing on the plot)
male.comp.data[1:5,4] <- c(NA,NA,NA,NA,NA)
male.comp.data[7:11,4] <- c(NA,NA,NA,NA,NA)
male.comp.data[13:17,4] <- c(NA,NA,NA,NA,NA)
male.comp.data[19:23,4] <- c(NA,NA,NA,NA,NA)
male.comp.data[25:29,4] <- c(NA,NA,NA,NA,NA)
male.comp.data[31:35,4] <- c(NA,NA,NA,NA,NA)

# add new "maxes" column to help with positioning the sig letters on the plot
male.comp.data[1:6,5] <- rep(max(male.comp.data[1:6,3]),6)
male.comp.data[7:12,5] <- rep(max(male.comp.data[7:12,3]),6)
male.comp.data[13:18,5] <- rep(max(male.comp.data[13:18,3]),6)
male.comp.data[19:24,5] <- rep(max(male.comp.data[19:24,3]),6)
male.comp.data[25:30,5] <- rep(max(male.comp.data[25:30,3]),6)
male.comp.data[31:36,5] <- rep(max(male.comp.data[31:36,3]),6)
colnames(male.comp.data)[5] <- "Maxes"

ggplot(data = male.comp.data, aes(x = Group, y = Prop.mated)) +
  theme_bw() +
  theme(panel.background = element_blank(), legend.text = element_text(size=12,colour = "black",face="bold"), axis.line = element_line(size = 0.1, colour = "black"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.x = element_text(size=15,face="bold",margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(size=15,face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text.y = element_text(size=12,colour = "black",face="bold"), axis.text.x = element_text(size=11,color="black",face="bold"), strip.text.x= element_text(size=15)) +
  labs(y="Male sexual competitiveness",x="Comparisons") +
  scale_y_continuous(limits = c(0, 1),expand=c(0.008,0.008)) +
  geom_point(alpha=1,shape=20,size=4,color=c(rep("#5CACEE",6),rep("#1874CD",6),rep("#35C47D",6),rep("#199151",6),rep("#CD96CD",6),rep("#8B668B",6))) +
  stat_summary(fun.data=mean_se, colour = c("#CD96CD","#8B668B","#35C47D","#199151","#5CACEE","#1874CD"), geom="errorbar", size=0.7, width=0.35) + 
  stat_summary(fun.data=mean_se, colour = c("#CD96CD","#8B668B","#35C47D","#199151","#5CACEE","#1874CD"), geom="point", size = 9, shape="_") +
  geom_text(aes(label = Letters, y=Maxes), vjust = -4,size=4.5,color="black")

#----------------------------------------------------------------------
# DESeq2 analysis to identify differentially-expressed genes
#----------------------------------------------------------------------

# Import the metadata and count data .csv files into R #
metadata <- read.csv("~/Desktop/RNAseq/RNAseq_metadata_groupall.csv", header=T, row.names=1)
countdata <- read.csv("~/Desktop/RNAseq/RNAseq_genecounts.csv", header=T, row.names=1)

# Run DESeq2

dds_full <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~fullgroup)
dds_full <- DESeq(dds_full)

dds_adults_G3F <- DESeqDataSetFromMatrix(countData = countdata[,c(2,4,6,8,10,12,14,16,18)], colData = metadata[c(2,4,6,8,10,12,14,16,18),], design = ~fullgroup)
dds_adults_G3F <- DESeq(dds_adults_G3F)

dds_adults_G3M <- DESeqDataSetFromMatrix(countData = countdata[,c(1,3,5,7,9,11,13,15,17)], colData = metadata[c(1,3,5,7,9,11,13,15,17),], design = ~fullgroup)
dds_adults_G3M <- DESeq(dds_adults_G3M)

dds_adults_G4F <- DESeqDataSetFromMatrix(countData = countdata[,c(20,22,24,26,28,30,32,34,36)], colData = metadata[c(20,22,24,26,28,30,32,34,36),], design = ~fullgroup)
dds_adults_G4F <- DESeq(dds_adults_G4F)

dds_adults_G4M <- DESeqDataSetFromMatrix(countData = countdata[,c(19,21,23,25,27,29,31,33,35)], colData = metadata[c(19,21,23,25,27,29,31,33,35),], design = ~fullgroup)
dds_adults_G4M <- DESeq(dds_adults_G4M)

dds_larvae <- DESeqDataSetFromMatrix(countData = countdata[,37:72], colData = metadata[37:72,], design = ~fullgroup)
dds_larvae <- DESeq(dds_larvae)

# Next, generate lists of DEGs associated with each comparison we care about
# We will first run the results function with DESeq2 for the 16 comparisons of interest

res_DOXG3_H2OG3_AM <- results(dds_adults_G3M, contrast=c("fullgroup","A_DOX_G3_M","A_H2O_G3_M"), alpha=0.05)
res_ATCG3_H2OG3_AM <- results(dds_adults_G3M, contrast=c("fullgroup","A_ATC_G3_M","A_H2O_G3_M"), alpha=0.05)
res_DOXG4_H2OG4_AM <- results(dds_adults_G4M, contrast=c("fullgroup","A_DOX_G4_M","A_H2O_G4_M"), alpha=0.05)
res_ATCG4_H2OG4_AM <- results(dds_adults_G4M, contrast=c("fullgroup","A_ATC_G4_M","A_H2O_G4_M"), alpha=0.05)
res_DOXG3_H2OG3_AF <- results(dds_adults_G3F, contrast=c("fullgroup","A_DOX_G3_F","A_H2O_G3_F"), alpha=0.05)
res_ATCG3_H2OG3_AF <- results(dds_adults_G3F, contrast=c("fullgroup","A_ATC_G3_F","A_H2O_G3_F"), alpha=0.05)
res_DOXG4_H2OG4_AF <- results(dds_adults_G4F, contrast=c("fullgroup","A_DOX_G4_F","A_H2O_G4_F"), alpha=0.05)
res_ATCG4_H2OG4_AF <- results(dds_adults_G4F, contrast=c("fullgroup","A_ATC_G4_F","A_H2O_G4_F"), alpha=0.05)
res_DOXG3_H2OG3_L <- results(dds_larvae, contrast=c("fullgroup","L_DOX_G3","L_H2O_G3"), alpha=0.05)
res_ATCG3_H2OG3_L <- results(dds_larvae, contrast=c("fullgroup","L_ATC_G3","L_H2O_G3"), alpha=0.05)
res_DOXG4_H2OG4_L <- results(dds_larvae, contrast=c("fullgroup","L_DOX_G4","L_H2O_G4"), alpha=0.05)
res_ATCG4_H2OG4_L <- results(dds_larvae, contrast=c("fullgroup","L_ATC_G4","L_H2O_G4"), alpha=0.05)

res_DOXG3_H2OG3_AM_sig <- subset(res_DOXG3_H2OG3_AM, padj < 0.05)
res_ATCG3_H2OG3_AM_sig <- subset(res_ATCG3_H2OG3_AM, padj < 0.05)
res_DOXG4_H2OG4_AM_sig <- subset(res_DOXG4_H2OG4_AM, padj < 0.05)
res_ATCG4_H2OG4_AM_sig <- subset(res_ATCG4_H2OG4_AM, padj < 0.05)
res_DOXG3_H2OG3_AF_sig <- subset(res_DOXG3_H2OG3_AF, padj < 0.05)
res_ATCG3_H2OG3_AF_sig <- subset(res_ATCG3_H2OG3_AF, padj < 0.05)
res_DOXG4_H2OG4_AF_sig <- subset(res_DOXG4_H2OG4_AF, padj < 0.05)
res_ATCG4_H2OG4_AF_sig <- subset(res_ATCG4_H2OG4_AF, padj < 0.05)
res_DOXG3_H2OG3_L_sig <- subset(res_DOXG3_H2OG3_L, padj < 0.05)
res_ATCG3_H2OG3_L_sig <- subset(res_ATCG3_H2OG3_L, padj < 0.05)
res_DOXG4_H2OG4_L_sig <- subset(res_DOXG4_H2OG4_L, padj < 0.05)
res_ATCG4_H2OG4_L_sig <- subset(res_ATCG4_H2OG4_L, padj < 0.05)

# next, we'll get the number of DEGs with adjusted p-values < 0.05
dim(res_DOXG3_H2OG3_L_sig)
dim(res_ATCG3_H2OG3_L_sig)
dim(res_DOXG4_H2OG4_L_sig)
dim(res_ATCG4_H2OG4_L_sig)
dim(res_DOXG3_H2OG3_AM_sig)
dim(res_DOXG3_H2OG3_AF_sig)
dim(res_ATCG3_H2OG3_AM_sig)
dim(res_ATCG3_H2OG3_AF_sig)
dim(res_DOXG4_H2OG4_AM_sig)
dim(res_DOXG4_H2OG4_AF_sig)
dim(res_ATCG4_H2OG4_AM_sig)
dim(res_ATCG4_H2OG4_AF_sig)

# Export the normalized gene counts
write.csv(counts((dds_full),normalized=TRUE),"~/Desktop/RNAseq/Data/gene_counts_normalized.csv")

### code for filtering out low gene counts after running DESeq2 results function comparing 2 groups
### this will eliminate genes on a comparison-by-comparison basis 
### to be included in downstream analysis, of the 6 samples, at least 4 must contain counts >= 5 (i.e. majority of samples must have counts >= 5)

# Get samples relevant for this comparison
chopping_block_samples_DOXG3_L_v_H2OG3_L <- rownames(metadata[metadata$stage == "Larva" & (metadata$gen_diet == "G3_DOX" | metadata$gen_diet == "G3_H2O"),])
chopping_block_samples_ATCG3_L_v_H2OG3_L <- rownames(metadata[metadata$stage == "Larva" & (metadata$gen_diet == "G3_ATC" | metadata$gen_diet == "G3_H2O"),])
chopping_block_samples_DOXG4_L_v_H2OG4_L <- rownames(metadata[metadata$stage == "Larva" & (metadata$gen_diet == "G4_DOX" | metadata$gen_diet == "G4_H2O"),])
chopping_block_samples_ATCG4_L_v_H2OG4_L <- rownames(metadata[metadata$stage == "Larva" & (metadata$gen_diet == "G4_ATC" | metadata$gen_diet == "G4_H2O"),])
chopping_block_samples_DOXG3_AM_v_H2OG3_AM <- rownames(metadata[metadata$sex_stage == "AM" & (metadata$gen_diet == "G3_DOX" | metadata$gen_diet == "G3_H2O"),])
chopping_block_samples_ATCG3_AM_v_H2OG3_AM <- rownames(metadata[metadata$sex_stage == "AM" & (metadata$gen_diet == "G3_ATC" | metadata$gen_diet == "G3_H2O"),])
chopping_block_samples_DOXG4_AM_v_H2OG4_AM <- rownames(metadata[metadata$sex_stage == "AM" & (metadata$gen_diet == "G4_DOX" | metadata$gen_diet == "G4_H2O"),])
chopping_block_samples_ATCG4_AM_v_H2OG4_AM <- rownames(metadata[metadata$sex_stage == "AM" & (metadata$gen_diet == "G4_ATC" | metadata$gen_diet == "G4_H2O"),])
chopping_block_samples_DOXG3_AF_v_H2OG3_AF <- rownames(metadata[metadata$sex_stage == "AF" & (metadata$gen_diet == "G3_DOX" | metadata$gen_diet == "G3_H2O"),])
chopping_block_samples_ATCG3_AF_v_H2OG3_AF <- rownames(metadata[metadata$sex_stage == "AF" & (metadata$gen_diet == "G3_ATC" | metadata$gen_diet == "G3_H2O"),])
chopping_block_samples_DOXG4_AF_v_H2OG4_AF <- rownames(metadata[metadata$sex_stage == "AF" & (metadata$gen_diet == "G4_DOX" | metadata$gen_diet == "G4_H2O"),])
chopping_block_samples_ATCG4_AF_v_H2OG4_AF <- rownames(metadata[metadata$sex_stage == "AF" & (metadata$gen_diet == "G4_ATC" | metadata$gen_diet == "G4_H2O"),])

# get counts for samples we care about
chopping_block_counts_DOXG3_L_v_H2OG3_L <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_DOXG3_L_v_H2OG3_L]
chopping_block_counts_ATCG3_L_v_H2OG3_L <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_ATCG3_L_v_H2OG3_L]
chopping_block_counts_DOXG4_L_v_H2OG4_L <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_DOXG4_L_v_H2OG4_L]
chopping_block_counts_ATCG4_L_v_H2OG4_L <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_ATCG4_L_v_H2OG4_L]
chopping_block_counts_DOXG3_AM_v_H2OG3_AM <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_DOXG3_AM_v_H2OG3_AM]
chopping_block_counts_ATCG3_AM_v_H2OG3_AM <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_ATCG3_AM_v_H2OG3_AM]
chopping_block_counts_DOXG4_AM_v_H2OG4_AM <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_DOXG4_AM_v_H2OG4_AM]
chopping_block_counts_ATCG4_AM_v_H2OG4_AM <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_ATCG4_AM_v_H2OG4_AM]
chopping_block_counts_DOXG3_AF_v_H2OG3_AF <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_DOXG3_AF_v_H2OG3_AF]
chopping_block_counts_ATCG3_AF_v_H2OG3_AF <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_ATCG3_AF_v_H2OG3_AF]
chopping_block_counts_DOXG4_AF_v_H2OG4_AF <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_DOXG4_AF_v_H2OG4_AF]
chopping_block_counts_ATCG4_AF_v_H2OG4_AF <- counts(dds_full, normalized=TRUE)[,chopping_block_samples_ATCG4_AF_v_H2OG4_AF]

# reduce it to genes we care about
chopping_block_counts_DOXG3_L_v_H2OG3_L <- chopping_block_counts_DOXG3_L_v_H2OG3_L[rownames(res_DOXG3_H2OG3_L_sig),,drop = FALSE]
chopping_block_counts_ATCG3_L_v_H2OG3_L <- chopping_block_counts_ATCG3_L_v_H2OG3_L[rownames(res_ATCG3_H2OG3_L_sig),,drop = FALSE]
chopping_block_counts_DOXG4_L_v_H2OG4_L <- chopping_block_counts_DOXG4_L_v_H2OG4_L[rownames(res_DOXG4_H2OG4_L_sig),,drop = FALSE]
chopping_block_counts_ATCG4_L_v_H2OG4_L <- chopping_block_counts_ATCG4_L_v_H2OG4_L[rownames(res_ATCG4_H2OG4_L_sig),,drop = FALSE]
chopping_block_counts_DOXG3_AM_v_H2OG3_AM <- chopping_block_counts_DOXG3_AM_v_H2OG3_AM[rownames(res_DOXG3_H2OG3_AM_sig),,drop = FALSE]
chopping_block_counts_ATCG3_AM_v_H2OG3_AM <- chopping_block_counts_ATCG3_AM_v_H2OG3_AM[rownames(res_ATCG3_H2OG3_AM_sig),,drop = FALSE]
chopping_block_counts_DOXG4_AM_v_H2OG4_AM <- chopping_block_counts_DOXG4_AM_v_H2OG4_AM[rownames(res_DOXG4_H2OG4_AM_sig),,drop = FALSE]
chopping_block_counts_ATCG4_AM_v_H2OG4_AM <- chopping_block_counts_ATCG4_AM_v_H2OG4_AM[rownames(res_ATCG4_H2OG4_AM_sig),,drop = FALSE]
chopping_block_counts_DOXG3_AF_v_H2OG3_AF <- chopping_block_counts_DOXG3_AF_v_H2OG3_AF[rownames(res_DOXG3_H2OG3_AF_sig),,drop = FALSE]
chopping_block_counts_ATCG3_AF_v_H2OG3_AF <- chopping_block_counts_ATCG3_AF_v_H2OG3_AF[rownames(res_ATCG3_H2OG3_AF_sig),,drop = FALSE]
chopping_block_counts_DOXG4_AF_v_H2OG4_AF <- chopping_block_counts_DOXG4_AF_v_H2OG4_AF[rownames(res_DOXG4_H2OG4_AF_sig),,drop = FALSE]
chopping_block_counts_ATCG4_AF_v_H2OG4_AF <- chopping_block_counts_ATCG4_AF_v_H2OG4_AF[rownames(res_ATCG4_H2OG4_AF_sig),,drop = FALSE]

# Cut out genes where majority of counts for a given comparison are zeroes
chopping_block_counts_DOXG3_L_v_H2OG3_L <- chopping_block_counts_DOXG3_L_v_H2OG3_L[rowSums(chopping_block_counts_DOXG3_L_v_H2OG3_L == 0) < 8,,drop = FALSE]
chopping_block_counts_DOXG3_AM_v_H2OG3_AM <- chopping_block_counts_DOXG3_AM_v_H2OG3_AM[rowSums(chopping_block_counts_DOXG3_AM_v_H2OG3_AM == 0) < 4,,drop = FALSE]
chopping_block_counts_DOXG3_AF_v_H2OG3_AF <- chopping_block_counts_DOXG3_AF_v_H2OG3_AF[rowSums(chopping_block_counts_DOXG3_AF_v_H2OG3_AF == 0) < 4,,drop = FALSE]
chopping_block_counts_ATCG3_L_v_H2OG3_L <- chopping_block_counts_ATCG3_L_v_H2OG3_L[rowSums(chopping_block_counts_ATCG3_L_v_H2OG3_L == 0) < 8,,drop = FALSE]
chopping_block_counts_ATCG3_AM_v_H2OG3_AM <- chopping_block_counts_ATCG3_AM_v_H2OG3_AM[rowSums(chopping_block_counts_ATCG3_AM_v_H2OG3_AM == 0) < 4,,drop = FALSE]
chopping_block_counts_ATCG3_AF_v_H2OG3_AF <- chopping_block_counts_ATCG3_AF_v_H2OG3_AF[rowSums(chopping_block_counts_ATCG3_AF_v_H2OG3_AF == 0) < 4,,drop = FALSE]
chopping_block_counts_DOXG4_L_v_H2OG4_L <- chopping_block_counts_DOXG4_L_v_H2OG4_L[rowSums(chopping_block_counts_DOXG4_L_v_H2OG4_L == 0) < 8,,drop = FALSE]
chopping_block_counts_DOXG4_AM_v_H2OG4_AM <- chopping_block_counts_DOXG4_AM_v_H2OG4_AM[rowSums(chopping_block_counts_DOXG4_AM_v_H2OG4_AM == 0) < 4,,drop = FALSE]
chopping_block_counts_DOXG4_AF_v_H2OG4_AF <- chopping_block_counts_DOXG4_AF_v_H2OG4_AF[rowSums(chopping_block_counts_DOXG4_AF_v_H2OG4_AF == 0) < 4,,drop = FALSE]
chopping_block_counts_ATCG4_L_v_H2OG4_L <- chopping_block_counts_ATCG4_L_v_H2OG4_L[rowSums(chopping_block_counts_ATCG4_L_v_H2OG4_L == 0) < 8,,drop = FALSE]
chopping_block_counts_ATCG4_AM_v_H2OG4_AM <- chopping_block_counts_ATCG4_AM_v_H2OG4_AM[rowSums(chopping_block_counts_ATCG4_AM_v_H2OG4_AM == 0) < 4,,drop = FALSE]
chopping_block_counts_ATCG4_AF_v_H2OG4_AF <- chopping_block_counts_ATCG4_AF_v_H2OG4_AF[rowSums(chopping_block_counts_ATCG4_AF_v_H2OG4_AF == 0) < 4,,drop = FALSE]

# Cut out genes where majority of positive counts are less than five
chopping_block_counts_DOXG3_L_v_H2OG3_L <- chopping_block_counts_DOXG3_L_v_H2OG3_L[rowSums(chopping_block_counts_DOXG3_L_v_H2OG3_L < 5) < 8,,drop = FALSE]
chopping_block_counts_DOXG3_AM_v_H2OG3_AM <- chopping_block_counts_DOXG3_AM_v_H2OG3_AM[rowSums(chopping_block_counts_DOXG3_AM_v_H2OG3_AM < 5) < 4,,drop = FALSE]
chopping_block_counts_DOXG3_AF_v_H2OG3_AF <- chopping_block_counts_DOXG3_AF_v_H2OG3_AF[rowSums(chopping_block_counts_DOXG3_AF_v_H2OG3_AF < 5) < 4,,drop = FALSE]
chopping_block_counts_ATCG3_L_v_H2OG3_L <- chopping_block_counts_ATCG3_L_v_H2OG3_L[rowSums(chopping_block_counts_ATCG3_L_v_H2OG3_L < 5) < 8,,drop = FALSE]
chopping_block_counts_ATCG3_AM_v_H2OG3_AM <- chopping_block_counts_ATCG3_AM_v_H2OG3_AM[rowSums(chopping_block_counts_ATCG3_AM_v_H2OG3_AM < 5) < 4,,drop = FALSE]
chopping_block_counts_ATCG3_AF_v_H2OG3_AF <- chopping_block_counts_ATCG3_AF_v_H2OG3_AF[rowSums(chopping_block_counts_ATCG3_AF_v_H2OG3_AF < 5) < 4,,drop = FALSE]
chopping_block_counts_DOXG4_L_v_H2OG4_L <- chopping_block_counts_DOXG4_L_v_H2OG4_L[rowSums(chopping_block_counts_DOXG4_L_v_H2OG4_L < 5) < 8,,drop = FALSE]
chopping_block_counts_DOXG4_AM_v_H2OG4_AM <- chopping_block_counts_DOXG4_AM_v_H2OG4_AM[rowSums(chopping_block_counts_DOXG4_AM_v_H2OG4_AM < 5) < 4,,drop = FALSE]
chopping_block_counts_DOXG4_AF_v_H2OG4_AF <- chopping_block_counts_DOXG4_AF_v_H2OG4_AF[rowSums(chopping_block_counts_DOXG4_AF_v_H2OG4_AF < 5) < 4,,drop = FALSE]
chopping_block_counts_ATCG4_L_v_H2OG4_L <- chopping_block_counts_ATCG4_L_v_H2OG4_L[rowSums(chopping_block_counts_ATCG4_L_v_H2OG4_L < 5) < 8,,drop = FALSE]
chopping_block_counts_ATCG4_AM_v_H2OG4_AM <- chopping_block_counts_ATCG4_AM_v_H2OG4_AM[rowSums(chopping_block_counts_ATCG4_AM_v_H2OG4_AM < 5) < 4,,drop = FALSE]
chopping_block_counts_ATCG4_AF_v_H2OG4_AF <- chopping_block_counts_ATCG4_AF_v_H2OG4_AF[rowSums(chopping_block_counts_ATCG4_AF_v_H2OG4_AF < 5) < 4,,drop = FALSE]

# Rerun code
res_DOXG3_H2OG3_L_sig[rownames(chopping_block_counts_DOXG3_L_v_H2OG3_L),]
res_DOXG3_H2OG3_AM_sig[rownames(chopping_block_counts_DOXG3_AM_v_H2OG3_AM),]
res_DOXG3_H2OG3_AF_sig[rownames(chopping_block_counts_DOXG3_AF_v_H2OG3_AF),]
res_DOXG4_H2OG4_L_sig[rownames(chopping_block_counts_DOXG4_L_v_H2OG4_L),]
res_DOXG4_H2OG4_AM_sig[rownames(chopping_block_counts_DOXG4_AM_v_H2OG4_AM),]
res_DOXG4_H2OG4_AF_sig[rownames(chopping_block_counts_DOXG4_AF_v_H2OG4_AF),]
res_ATCG3_H2OG3_L_sig[rownames(chopping_block_counts_ATCG3_L_v_H2OG3_L),]
res_ATCG3_H2OG3_AM_sig[rownames(chopping_block_counts_ATCG3_AM_v_H2OG3_AM),]
res_ATCG3_H2OG3_AF_sig[rownames(chopping_block_counts_ATCG3_AF_v_H2OG3_AF),]
res_ATCG4_H2OG4_L_sig[rownames(chopping_block_counts_ATCG4_L_v_H2OG4_L),]
res_ATCG4_H2OG4_AM_sig[rownames(chopping_block_counts_ATCG4_AM_v_H2OG4_AM),]
res_ATCG4_H2OG4_AF_sig[rownames(chopping_block_counts_ATCG4_AF_v_H2OG4_AF),]

# Export 4 files for all 16 group comparisons
# list of all filtered DEGs (p adj < 0.05, gene count filter applied)
# list of downregulated, filtered DEGs
# list of upregulated, filtered DEGs
# unfiltered results from DESeq2 group comparison (all genes)

res_DOXG3_H2OG3_AM_ordered <- res_DOXG3_H2OG3_AM_sig[rownames(chopping_block_counts_DOXG3_AM_v_H2OG3_AM),][order(res_DOXG3_H2OG3_AM_sig[rownames(chopping_block_counts_DOXG3_AM_v_H2OG3_AM),]$log2FoldChange),]
downreg_res_DOXG3_H2OG3_AM_ordered <- res_DOXG3_H2OG3_AM_ordered[res_DOXG3_H2OG3_AM_ordered$log2FoldChange < 0,]
upreg_res_DOXG3_H2OG3_AM_ordered <- res_DOXG3_H2OG3_AM_ordered[res_DOXG3_H2OG3_AM_ordered$log2FoldChange > 0,]
upreg_res_DOXG3_H2OG3_AM_ordered <- upreg_res_DOXG3_H2OG3_AM_ordered[order(upreg_res_DOXG3_H2OG3_AM_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_DOXG3_H2OG3_AM_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_DOXG3_H2OG3_AM.csv")
write.csv(rownames(downreg_res_DOXG3_H2OG3_AM_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_DOXG3_AM.csv")
write.csv(rownames(upreg_res_DOXG3_H2OG3_AM_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_DOXG3_AM.csv")
write.csv(res_DOXG3_H2OG3_AM,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG3_H2OG3_AM.csv")

res_ATCG3_H2OG3_AM_ordered <- res_ATCG3_H2OG3_AM_sig[rownames(chopping_block_counts_ATCG3_AM_v_H2OG3_AM),][order(res_ATCG3_H2OG3_AM_sig[rownames(chopping_block_counts_ATCG3_AM_v_H2OG3_AM),]$log2FoldChange),]
downreg_res_ATCG3_H2OG3_AM_ordered <- res_ATCG3_H2OG3_AM_ordered[res_ATCG3_H2OG3_AM_ordered$log2FoldChange < 0,]
upreg_res_ATCG3_H2OG3_AM_ordered <- res_ATCG3_H2OG3_AM_ordered[res_ATCG3_H2OG3_AM_ordered$log2FoldChange > 0,]
upreg_res_ATCG3_H2OG3_AM_ordered <- upreg_res_ATCG3_H2OG3_AM_ordered[order(upreg_res_ATCG3_H2OG3_AM_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_ATCG3_H2OG3_AM_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_ATCG3_H2OG3_AM.csv")
write.csv(rownames(downreg_res_ATCG3_H2OG3_AM_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_ATCG3_AM.csv")
write.csv(rownames(upreg_res_ATCG3_H2OG3_AM_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_ATCG3_AM.csv")
write.csv(res_ATCG3_H2OG3_AM,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG3_H2OG3_AM.csv")

res_DOXG4_H2OG4_AM_ordered <- res_DOXG4_H2OG4_AM_sig[rownames(chopping_block_counts_DOXG4_AM_v_H2OG4_AM),][order(res_DOXG4_H2OG4_AM_sig[rownames(chopping_block_counts_DOXG4_AM_v_H2OG4_AM),]$log2FoldChange),]
downreg_res_DOXG4_H2OG4_AM_ordered <- res_DOXG4_H2OG4_AM_ordered[res_DOXG4_H2OG4_AM_ordered$log2FoldChange < 0,]
upreg_res_DOXG4_H2OG4_AM_ordered <- res_DOXG4_H2OG4_AM_ordered[res_DOXG4_H2OG4_AM_ordered$log2FoldChange > 0,]
upreg_res_DOXG4_H2OG4_AM_ordered <- upreg_res_DOXG4_H2OG4_AM_ordered[order(upreg_res_DOXG4_H2OG4_AM_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_DOXG4_H2OG4_AM_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_DOXG4_H2OG4_AM.csv")
write.csv(rownames(downreg_res_DOXG4_H2OG4_AM_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_DOXG4_AM.csv")
write.csv(rownames(upreg_res_DOXG4_H2OG4_AM_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_DOXG4_AM.csv")
write.csv(res_DOXG4_H2OG4_AM,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG4_H2OG4_AM.csv")

res_ATCG4_H2OG4_AM_ordered <- res_ATCG4_H2OG4_AM_sig[rownames(chopping_block_counts_ATCG4_AM_v_H2OG4_AM),][order(res_ATCG4_H2OG4_AM_sig[rownames(chopping_block_counts_ATCG4_AM_v_H2OG4_AM),]$log2FoldChange),]
downreg_res_ATCG4_H2OG4_AM_ordered <- res_ATCG4_H2OG4_AM_ordered[res_ATCG4_H2OG4_AM_ordered$log2FoldChange < 0,]
upreg_res_ATCG4_H2OG4_AM_ordered <- res_ATCG4_H2OG4_AM_ordered[res_ATCG4_H2OG4_AM_ordered$log2FoldChange > 0,]
upreg_res_ATCG4_H2OG4_AM_ordered <- upreg_res_ATCG4_H2OG4_AM_ordered[order(upreg_res_ATCG4_H2OG4_AM_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_ATCG4_H2OG4_AM_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_ATCG4_H2OG4_AM.csv")
write.csv(rownames(downreg_res_ATCG4_H2OG4_AM_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_ATCG4_AM.csv")
write.csv(rownames(upreg_res_ATCG4_H2OG4_AM_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_ATCG4_AM.csv")
write.csv(res_ATCG4_H2OG4_AM,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG4_H2OG4_AM.csv")

res_DOXG3_H2OG3_AF_ordered <- res_DOXG3_H2OG3_AF_sig[rownames(chopping_block_counts_DOXG3_AF_v_H2OG3_AF),][order(res_DOXG3_H2OG3_AF_sig[rownames(chopping_block_counts_DOXG3_AF_v_H2OG3_AF),]$log2FoldChange),]
downreg_res_DOXG3_H2OG3_AF_ordered <- res_DOXG3_H2OG3_AF_ordered[res_DOXG3_H2OG3_AF_ordered$log2FoldChange < 0,]
upreg_res_DOXG3_H2OG3_AF_ordered <- res_DOXG3_H2OG3_AF_ordered[res_DOXG3_H2OG3_AF_ordered$log2FoldChange > 0,]
upreg_res_DOXG3_H2OG3_AF_ordered <- upreg_res_DOXG3_H2OG3_AF_ordered[order(upreg_res_DOXG3_H2OG3_AF_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_DOXG3_H2OG3_AF_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_DOXG3_H2OG3_AF.csv")
write.csv(rownames(downreg_res_DOXG3_H2OG3_AF_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_DOXG3_AF.csv")
write.csv(rownames(upreg_res_DOXG3_H2OG3_AF_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_DOXG3_AF.csv")
write.csv(res_DOXG3_H2OG3_AF,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG3_H2OG3_AF.csv")

res_ATCG3_H2OG3_AF_ordered <- res_ATCG3_H2OG3_AF_sig[rownames(chopping_block_counts_ATCG3_AF_v_H2OG3_AF),][order(res_ATCG3_H2OG3_AF_sig[rownames(chopping_block_counts_ATCG3_AF_v_H2OG3_AF),]$log2FoldChange),]
downreg_res_ATCG3_H2OG3_AF_ordered <- res_ATCG3_H2OG3_AF_ordered[res_ATCG3_H2OG3_AF_ordered$log2FoldChange < 0,]
upreg_res_ATCG3_H2OG3_AF_ordered <- res_ATCG3_H2OG3_AF_ordered[res_ATCG3_H2OG3_AF_ordered$log2FoldChange > 0,]
upreg_res_ATCG3_H2OG3_AF_ordered <- upreg_res_ATCG3_H2OG3_AF_ordered[order(upreg_res_ATCG3_H2OG3_AF_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_ATCG3_H2OG3_AF_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_ATCG3_H2OG3_AF.csv")
write.csv(rownames(downreg_res_ATCG3_H2OG3_AF_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_ATCG3_AF.csv")
write.csv(rownames(upreg_res_ATCG3_H2OG3_AF_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_ATCG3_AF.csv")
write.csv(res_ATCG3_H2OG3_AF,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG3_H2OG3_AF.csv")

res_DOXG4_H2OG4_AF_ordered <- res_DOXG4_H2OG4_AF_sig[rownames(chopping_block_counts_DOXG4_AF_v_H2OG4_AF),][order(res_DOXG4_H2OG4_AF_sig[rownames(chopping_block_counts_DOXG4_AF_v_H2OG4_AF),]$log2FoldChange),]
downreg_res_DOXG4_H2OG4_AF_ordered <- res_DOXG4_H2OG4_AF_ordered[res_DOXG4_H2OG4_AF_ordered$log2FoldChange < 0,]
upreg_res_DOXG4_H2OG4_AF_ordered <- res_DOXG4_H2OG4_AF_ordered[res_DOXG4_H2OG4_AF_ordered$log2FoldChange > 0,]
upreg_res_DOXG4_H2OG4_AF_ordered <- upreg_res_DOXG4_H2OG4_AF_ordered[order(upreg_res_DOXG4_H2OG4_AF_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_DOXG4_H2OG4_AF_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_DOXG4_H2OG4_AF.csv")
write.csv(rownames(downreg_res_DOXG4_H2OG4_AF_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_DOXG4_AF.csv")
write.csv(rownames(upreg_res_DOXG4_H2OG4_AF_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_DOXG4_AF.csv")
write.csv(res_DOXG4_H2OG4_AF,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG4_H2OG4_AF.csv")

res_ATCG4_H2OG4_AF_ordered <- res_ATCG4_H2OG4_AF_sig[rownames(chopping_block_counts_ATCG4_AF_v_H2OG4_AF),][order(res_ATCG4_H2OG4_AF_sig[rownames(chopping_block_counts_ATCG4_AF_v_H2OG4_AF),]$log2FoldChange),]
downreg_res_ATCG4_H2OG4_AF_ordered <- res_ATCG4_H2OG4_AF_ordered[res_ATCG4_H2OG4_AF_ordered$log2FoldChange < 0,]
upreg_res_ATCG4_H2OG4_AF_ordered <- res_ATCG4_H2OG4_AF_ordered[res_ATCG4_H2OG4_AF_ordered$log2FoldChange > 0,]
upreg_res_ATCG4_H2OG4_AF_ordered <- upreg_res_ATCG4_H2OG4_AF_ordered[order(upreg_res_ATCG4_H2OG4_AF_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_ATCG4_H2OG4_AF_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_ATCG4_H2OG4_AF.csv")
write.csv(rownames(downreg_res_ATCG4_H2OG4_AF_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_ATCG4_AF.csv")
write.csv(rownames(upreg_res_ATCG4_H2OG4_AF_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_ATCG4_AF.csv")
write.csv(res_ATCG4_H2OG4_AF,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG4_H2OG4_AF.csv")

res_DOXG3_H2OG3_L_ordered <- res_DOXG3_H2OG3_L_sig[rownames(chopping_block_counts_DOXG3_L_v_H2OG3_L),][order(res_DOXG3_H2OG3_L_sig[rownames(chopping_block_counts_DOXG3_L_v_H2OG3_L),]$log2FoldChange),]
downreg_res_DOXG3_H2OG3_L_ordered <- res_DOXG3_H2OG3_L_ordered[res_DOXG3_H2OG3_L_ordered$log2FoldChange < 0,]
upreg_res_DOXG3_H2OG3_L_ordered <- res_DOXG3_H2OG3_L_ordered[res_DOXG3_H2OG3_L_ordered$log2FoldChange > 0,]
upreg_res_DOXG3_H2OG3_L_ordered <- upreg_res_DOXG3_H2OG3_L_ordered[order(upreg_res_DOXG3_H2OG3_L_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_DOXG3_H2OG3_L_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_DOXG3_H2OG3_L.csv")
write.csv(rownames(downreg_res_DOXG3_H2OG3_L_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_DOXG3_L.csv")
write.csv(rownames(upreg_res_DOXG3_H2OG3_L_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_DOXG3_L.csv")
write.csv(res_DOXG3_H2OG3_L,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG3_H2OG3_L.csv")

res_ATCG3_H2OG3_L_ordered <- res_ATCG3_H2OG3_L_sig[rownames(chopping_block_counts_ATCG3_L_v_H2OG3_L),][order(res_ATCG3_H2OG3_L_sig[rownames(chopping_block_counts_ATCG3_L_v_H2OG3_L),]$log2FoldChange),]
downreg_res_ATCG3_H2OG3_L_ordered <- res_ATCG3_H2OG3_L_ordered[res_ATCG3_H2OG3_L_ordered$log2FoldChange < 0,]
upreg_res_ATCG3_H2OG3_L_ordered <- res_ATCG3_H2OG3_L_ordered[res_ATCG3_H2OG3_L_ordered$log2FoldChange > 0,]
upreg_res_ATCG3_H2OG3_L_ordered <- upreg_res_ATCG3_H2OG3_L_ordered[order(upreg_res_ATCG3_H2OG3_L_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_ATCG3_H2OG3_L_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_ATCG3_H2OG3_L.csv")
write.csv(rownames(downreg_res_ATCG3_H2OG3_L_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_ATCG3_L.csv")
write.csv(rownames(upreg_res_ATCG3_H2OG3_L_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_ATCG3_L.csv")
write.csv(res_ATCG3_H2OG3_L,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG3_H2OG3_L.csv")

res_DOXG4_H2OG4_L_ordered <- res_DOXG4_H2OG4_L_sig[rownames(chopping_block_counts_DOXG4_L_v_H2OG4_L),][order(res_DOXG4_H2OG4_L_sig[rownames(chopping_block_counts_DOXG4_L_v_H2OG4_L),]$log2FoldChange),]
downreg_res_DOXG4_H2OG4_L_ordered <- res_DOXG4_H2OG4_L_ordered[res_DOXG4_H2OG4_L_ordered$log2FoldChange < 0,]
upreg_res_DOXG4_H2OG4_L_ordered <- res_DOXG4_H2OG4_L_ordered[res_DOXG4_H2OG4_L_ordered$log2FoldChange > 0,]
upreg_res_DOXG4_H2OG4_L_ordered <- upreg_res_DOXG4_H2OG4_L_ordered[order(upreg_res_DOXG4_H2OG4_L_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_DOXG4_H2OG4_L_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_DOXG4_H2OG4_L.csv")
write.csv(rownames(downreg_res_DOXG4_H2OG4_L_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_DOXG4_L.csv")
write.csv(rownames(upreg_res_DOXG4_H2OG4_L_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_DOXG4_L.csv")
write.csv(res_DOXG4_H2OG4_L,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG4_H2OG4_L.csv")

res_ATCG4_H2OG4_L_ordered <- res_ATCG4_H2OG4_L_sig[rownames(chopping_block_counts_ATCG4_L_v_H2OG4_L),][order(res_ATCG4_H2OG4_L_sig[rownames(chopping_block_counts_ATCG4_L_v_H2OG4_L),]$log2FoldChange),]
downreg_res_ATCG4_H2OG4_L_ordered <- res_ATCG4_H2OG4_L_ordered[res_ATCG4_H2OG4_L_ordered$log2FoldChange < 0,]
upreg_res_ATCG4_H2OG4_L_ordered <- res_ATCG4_H2OG4_L_ordered[res_ATCG4_H2OG4_L_ordered$log2FoldChange > 0,]
upreg_res_ATCG4_H2OG4_L_ordered <- upreg_res_ATCG4_H2OG4_L_ordered[order(upreg_res_ATCG4_H2OG4_L_ordered$log2FoldChange,decreasing=TRUE),]
write.csv(res_ATCG4_H2OG4_L_ordered,"~/Desktop/RNAseq/Data/DEGs/Filtered/All/DEGs_ATCG4_H2OG4_L.csv")
write.csv(rownames(downreg_res_ATCG4_H2OG4_L_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Downregulated/downreg_ATCG4_L.csv")
write.csv(rownames(upreg_res_ATCG4_H2OG4_L_ordered),"~/Desktop/RNAseq/Data/DEGs/Filtered/Upregulated/upreg_ATCG4_L.csv")
write.csv(res_ATCG4_H2OG4_L,"~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG4_H2OG4_L.csv")

################################################################################################
# how many DEGs do the DOX and ATC groups have in common? (i.e. how similar are the 2 antibiotics in terms of their effects on gene expression?)
################################################################################################

downreg_DOXG3_ATCG3_shared <- intersect(DOXG3_downreg,ATCG3_downreg)
length(downreg_DOXG3_ATCG3_shared) # 143 shared downreg genes

upreg_DOXG3_ATCG3_shared <- intersect(DOXG3_upreg,ATCG3_upreg)
length(upreg_DOXG3_ATCG3_shared) # 324 shared upreg genes

downreg_DOXG4_ATCG4_shared <- intersect(DOXG4_downreg,ATCG4_downreg)
length(downreg_DOXG4_ATCG4_shared) # 19 shared downreg genes

upreg_DOXG4_ATCG4_shared <- intersect(DOXG4_upreg,ATCG4_upreg)
length(upreg_DOXG4_ATCG4_shared) # 7 shared upreg genes

### patterns similar across G3 and G4 for a given antibiotic?
downreg_DOXG3_DOXG4_shared <- intersect(DOXG3_downreg,DOXG4_downreg)
length(downreg_DOXG3_DOXG4_shared) # 0 shared downreg genes

upreg_DOXG3_DOXG4_shared <- intersect(DOXG3_upreg,DOXG4_upreg)
length(upreg_DOXG3_DOXG4_shared) # 5 shared downreg genes

downreg_ATCG3_ATCG4_shared <- intersect(ATCG3_downreg,ATCG4_downreg)
length(downreg_ATCG3_ATCG4_shared) # 9 shared downreg genes

upreg_ATCG3_ATCG4_shared <- intersect(ATCG3_upreg,ATCG4_upreg)
length(upreg_ATCG3_ATCG4_shared) # 12 shared downreg genes

# export .csv files
write.csv(upreg_DOXG3_ATCG3_shared,"~/Desktop/RNAseq/Data/DEGs/Filtered/Shared/All/upreg_DOXG3_ATCG3_shared.csv")
write.csv(upreg_DOXG4_ATCG4_shared,"~/Desktop/RNAseq/Data/DEGs/Filtered/Shared/All/upreg_DOXG4_ATCG4_shared.csv")
write.csv(downreg_DOXG3_ATCG3_shared,"~/Desktop/RNAseq/Data/DEGs/Filtered/Shared/All/downreg_DOXG3_ATCG3_shared.csv")
write.csv(downreg_DOXG4_ATCG4_shared,"~/Desktop/RNAseq/Data/DEGs/Filtered/Shared/All/downreg_DOXG4_ATCG4_shared.csv")
write.csv(upreg_DOXG3_DOXG4_shared,"~/Desktop/RNAseq/Data/DEGs/Filtered/Shared/All/upreg_DOXG3_DOXG4_shared.csv")
write.csv(upreg_ATCG3_ATCG4_shared,"~/Desktop/RNAseq/Data/DEGs/Filtered/Shared/All/upreg_ATCG3_ATCG4_shared.csv")
write.csv(downreg_DOXG3_DOXG4_shared,"~/Desktop/RNAseq/Data/DEGs/Filtered/Shared/All/downreg_DOXG3_DOXG4_shared.csv")
write.csv(downreg_ATCG3_ATCG4_shared,"~/Desktop/RNAseq/Data/DEGs/Filtered/Shared/All/downreg_ATCG3_ATCG4_shared.csv")

#----------------------------------------------------------------------
# Figures 3 and 4: Volcano plots for the 12 subgroup comparisons
#----------------------------------------------------------------------

install.packages("ggiraph")
library(ggiraph)

#### DOX G3 larvae #####

# import DESeq2 results for comparison of interest
DOXG3_L_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG3_H2OG3_L.csv")

# remove spurious genes with adj pvals < 0.05
DOXG3_L_volc_diff <- setdiff(rownames(res_DOXG3_H2OG3_L_sig),rownames(chopping_block_counts_DOXG3_L_v_H2OG3_L))
DOXG3_L_volc <- DOXG3_L_volc[!(DOXG3_L_volc$X %in% DOXG3_L_volc_diff),]

# remove unnecessary columns
DOXG3_L_volc <- DOXG3_L_volc[,c(1,3,7)]

# order by padj (low to high)
DOXG3_L_volc <- DOXG3_L_volc[order(DOXG3_L_volc$padj),]

# add new empty column for locus names:
DOXG3_L_volc$name <- rep("NA",nrow(DOXG3_L_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(DOXG3_L_volc[DOXG3_L_volc$padj < 0.05 & !(is.na(DOXG3_L_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_DOXG3_L.csv")

# import gene IDs
GeneIDs_DOXG3_L <- read.csv("GOI_DOXG3_L.csv")
DOXG3_L_volc[1:dim(DOXG3_L_volc[DOXG3_L_volc$padj < 0.05 & !(is.na(DOXG3_L_volc$padj)),])[1],4] <- GeneIDs_DOXG3_L$name

# make initial interactive volc plot to identify GOIs to be labeled
DOXG3_L_volcplot <- DOXG3_L_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  scale_color_identity() +
  scale_alpha_identity() +
  ggiraph::geom_point_interactive(aes(tooltip = sprintf("%s", name), hover_nearest = TRUE), size = 0.5) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="red",linetype="solid",alpha=0.4) +
  scale_y_continuous(limits = c(-0.25, 60), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-10, 10), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")
DOXG3_L_volcplot

ggiraph::girafe(ggobj = DOXG3_L_volcplot,
                options = list(
                  opts_tooltip(use_fill = TRUE),
                  opts_zoom(min = 0.5, max = 5),
                  opts_sizing(rescale = FALSE),
                  opts_toolbar(saveaspng = TRUE, delay_mouseout = 2000)
                )
)

# make locus IDs rownames so we can index by them
rownames(DOXG3_L_volc) <- DOXG3_L_volc$X

# set color and transparency values for genes
DOXG3_L_volc$color <- "gray40"
DOXG3_L_volc$alpha <- 0.5

DOXG3_L_selected <- c("LOC111681819","LOC111681823","LOC111678842","LOC111678840","LOC111682084","LOC111681336","LOC111682376","LOC111679190","LOC111677781","LOC111674949")

DOXG3_L_volc[DOXG3_L_selected,5] <- c(rep("#005eff",5),rep("firebrick1",5))
DOXG3_L_volc[DOXG3_L_selected,6] <- 1

DOXG3_L_volcplot <- DOXG3_L_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0.4)) +
  scale_x_continuous(limits = c(-10, 10), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = DOXG3_L_volc[DOXG3_L_selected,], aes(label = str_wrap(name, 21)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(1,10,18,22,30,-4,4,-5,9,-1),nudge_x=c(-4,-5.5,-4.5,-4.5,0,1,1,3,4,4))
DOXG3_L_volcplot

### DOXG3 adult males ###

# import DESeq2 results for comparison of interest
DOXG3_AM_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG3_H2OG3_AM.csv")

# remove spurious genes with adj pvals < 0.05
DOXG3_AM_volc_diff <- setdiff(rownames(res_DOXG3_H2OG3_AM_sig),rownames(chopping_block_counts_DOXG3_AM_v_H2OG3_AM))
DOXG3_AM_volc <- DOXG3_AM_volc[!(DOXG3_AM_volc$X %in% DOXG3_AM_volc_diff),]

# remove unnecessary columns
DOXG3_AM_volc <- DOXG3_AM_volc[,c(1,3,7)]

# order by padj (low to high)
DOXG3_AM_volc <- DOXG3_AM_volc[order(DOXG3_AM_volc$padj),]

# add new empty column for locus names:
DOXG3_AM_volc$name <- rep("NA",nrow(DOXG3_AM_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(DOXG3_AM_volc[DOXG3_AM_volc$padj < 0.05 & !(is.na(DOXG3_AM_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_DOXG3_AM.csv")

# import gene IDs
GeneIDs_DOXG3_AM <- read.csv("GOI_DOXG3_AM.csv")
DOXG3_AM_volc[1:dim(DOXG3_AM_volc[DOXG3_AM_volc$padj < 0.05 & !(is.na(DOXG3_AM_volc$padj)),])[1],4] <- GeneIDs_DOXG3_AM$name

# make initial interactive volc plot to identify GOIs to be labeled
DOXG3_AM_volcplot <- DOXG3_AM_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  scale_color_identity() +
  scale_alpha_identity() +
  ggiraph::geom_point_interactive(aes(tooltip = sprintf("%s", name), hover_nearest = TRUE), size = 0.5) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="red",linetype="solid",alpha=0.4) +
  scale_y_continuous(limits = c(-0.25, 30), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-10, 10), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")
DOXG3_AM_volcplot

ggiraph::girafe(ggobj = DOXG3_AM_volcplot,
                options = list(
                  opts_tooltip(use_fill = TRUE),
                  opts_zoom(min = 0.5, max = 5),
                  opts_sizing(rescale = FALSE),
                  opts_toolbar(saveaspng = TRUE, delay_mouseout = 2000)
                )
)

# make locus IDs rownames so we can index by them
rownames(DOXG3_AM_volc) <- DOXG3_AM_volc$X

# set color and transparency values for genes
DOXG3_AM_volc$color <- "gray40"
DOXG3_AM_volc$alpha <- 0.5

DOXG3_AM_selected <- c("LOC124419634","LOC111675917","LOC111690825","LOC124419393","LOC111676150","LOC111688454","LOC111688465","LOC111685676","LOC111682084","LOC111689384")

DOXG3_AM_volc[DOXG3_AM_selected,5] <- c(rep("#005eff",5),rep("firebrick1",5))
DOXG3_AM_volc[DOXG3_AM_selected,6] <- 1

# abbreviate as needed
#DOXG3_AM_volc[DOXG3_AM_selected,4][c(2,4,7)] <- c("LOC111681068","dendritic arbor reduction protein 1, transcript var. X2*","phosphoenolpyruvate carboxykinase [GTP]-like, transcript var. X1*")

DOXG3_AM_volcplot <- DOXG3_AM_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(-0.15, 25), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-10, 10), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = DOXG3_AM_volc[DOXG3_AM_selected,], aes(label = str_wrap(name, 25)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(0.5,2.4,1.5,6,4.8,-1,-1,-1,1,1),nudge_x=c(-4,-1,-5,-4,-4.7,0,1,-1.8,3,4))

DOXG3_AM_volcplot

#################### DOXG3_AF ###

# import DESeq2 results for comparison of interest
DOXG3_AF_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG3_H2OG3_AF.csv")

# remove spurious genes with adj pvals < 0.05
DOXG3_AF_volc_diff <- setdiff(rownames(res_DOXG3_H2OG3_AF_sig),rownames(chopping_block_counts_DOXG3_AF_v_H2OG3_AF))
DOXG3_AF_volc <- DOXG3_AF_volc[!(DOXG3_AF_volc$X %in% DOXG3_AF_volc_diff),]

# remove unnecessary columns
DOXG3_AF_volc <- DOXG3_AF_volc[,c(1,3,7)]

# order by padj (low to high)
DOXG3_AF_volc <- DOXG3_AF_volc[order(DOXG3_AF_volc$padj),]

# add new empty column for locus names:
DOXG3_AF_volc$name <- rep("NA",nrow(DOXG3_AF_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(DOXG3_AF_volc[DOXG3_AF_volc$padj < 0.05 & !(is.na(DOXG3_AF_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_DOXG3_AF.csv")

# import gene IDs
GeneIDs_DOXG3_AF <- read.csv("GOI_DOXG3_AF.csv")
DOXG3_AF_volc[1:dim(DOXG3_AF_volc[DOXG3_AF_volc$padj < 0.05 & !(is.na(DOXG3_AF_volc$padj)),])[1],4] <- GeneIDs_DOXG3_AF$name

# make locus IDs rownames so we can index by them
rownames(DOXG3_AF_volc) <- DOXG3_AF_volc$X

# set color and transparency values for genes
DOXG3_AF_volc$color <- "gray40"
DOXG3_AF_volc$alpha <- 0.5

DOXG3_AF_selected <- c("LOC111687783","LOC111675917","LOC111690825","LOC111681456","LOC111690237","LOC111677053","LOC111690454","LOC111686175","LOC111683933","LOC111676083")

DOXG3_AF_volc[DOXG3_AF_selected,5] <- c(rep("#005eff",5),rep("firebrick1",5))
DOXG3_AF_volc[DOXG3_AF_selected,6] <- 1

DOXG3_AF_volcplot <- DOXG3_AF_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(-0.10, 13), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-8, 8), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = DOXG3_AF_volc[DOXG3_AF_selected,], aes(label = str_wrap(name, 25)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(-0.2,1.8,1,5,3.6,2.6,2,2.4,3,1),nudge_x=c(-3,-4,-4,-2.1,-0.5,3,3.5,3,1,5))

DOXG3_AF_volcplot

############################## DOX G4 L

# import DESeq2 results for comparison of interest
DOXG4_L_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG4_H2OG4_L.csv")

# remove spurious genes with adj pvals < 0.05
DOXG4_L_volc_diff <- setdiff(rownames(res_DOXG4_H2OG4_L_sig),rownames(chopping_block_counts_DOXG4_L_v_H2OG4_L))
DOXG4_L_volc <- DOXG4_L_volc[!(DOXG4_L_volc$X %in% DOXG4_L_volc_diff),]

# remove unnecessary columns
DOXG4_L_volc <- DOXG4_L_volc[,c(1,3,7)]

# order by padj (low to high)
DOXG4_L_volc <- DOXG4_L_volc[order(DOXG4_L_volc$padj),]

# add new empty column for locus names:
DOXG4_L_volc$name <- rep("NA",nrow(DOXG4_L_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(DOXG4_L_volc[DOXG4_L_volc$padj < 0.05 & !(is.na(DOXG4_L_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_DOXG4_L.csv")

# import gene IDs
GeneIDs_DOXG4_L <- read.csv("GOI_DOXG4_L.csv")
DOXG4_L_volc[1:dim(DOXG4_L_volc[DOXG4_L_volc$padj < 0.05 & !(is.na(DOXG4_L_volc$padj)),])[1],4] <- GeneIDs_DOXG4_L$name

# make locus IDs rownames so we can index by them
rownames(DOXG4_L_volc) <- DOXG4_L_volc$X

# set color and transparency values for genes
DOXG4_L_volc$color <- "gray40"
DOXG4_L_volc$alpha <- 0.5

DOXG4_L_selected <- c("LOC111678332","LOC111676372","LOC111681706","LOC111690205","LOC111682470","LOC111678847","LOC111683242","LOC124419342","LOC111683610","LOC111676624")

DOXG4_L_volc[DOXG4_L_selected,5] <- c(rep("#005eff",5),rep("firebrick1",5))
DOXG4_L_volc[DOXG4_L_selected,6] <- 1

DOXG4_L_volcplot <- DOXG4_L_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(-0.07, 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-6, 6), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = DOXG4_L_volc[DOXG4_L_selected,], aes(label = str_wrap(name, 31)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(0,1.3,2.5,3.7,5,5.5,2.7,4,0.6,0.6),nudge_x=c(-3,-3,-4,-4,-3,2,2.3,2.5,0.8,0))
DOXG4_L_volcplot

#################### DOXG4_AM ###

# import DESeq2 results for comparison of interest
DOXG4_AM_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG4_H2OG4_AM.csv")

# remove spurious genes with adj pvals < 0.05
DOXG4_AM_volc_diff <- setdiff(rownames(res_DOXG4_H2OG4_AM_sig),rownames(chopping_block_counts_DOXG4_AM_v_H2OG4_AM))
DOXG4_AM_volc <- DOXG4_AM_volc[!(DOXG4_AM_volc$X %in% DOXG4_AM_volc_diff),]

# remove unnecessary columns
DOXG4_AM_volc <- DOXG4_AM_volc[,c(1,3,7)]

# order by padj (low to high)
DOXG4_AM_volc <- DOXG4_AM_volc[order(DOXG4_AM_volc$padj),]

# add new empty column for locus names:
DOXG4_AM_volc$name <- rep("NA",nrow(DOXG4_AM_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(DOXG4_AM_volc[DOXG4_AM_volc$padj < 0.05 & !(is.na(DOXG4_AM_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_DOXG4_AM.csv")

# import gene IDs
GeneIDs_DOXG4_AM <- read.csv("GOI_DOXG4_AM.csv")
DOXG4_AM_volc[1:dim(DOXG4_AM_volc[DOXG4_AM_volc$padj < 0.05 & !(is.na(DOXG4_AM_volc$padj)),])[1],4] <- GeneIDs_DOXG4_AM$name

# make locus IDs rownames so we can index by them
rownames(DOXG4_AM_volc) <- DOXG4_AM_volc$X

# set color and transparency values for genes
DOXG4_AM_volc$color <- "gray40"
DOXG4_AM_volc$alpha <- 0.5

DOXG4_AM_selected <- c("LOC111690922","LOC111675342","LOC111675507","LOC111689926","LOC124419342","LOC124418500","LOC111679046","LOC111690329","LOC111683360","LOC111674637")

DOXG4_AM_volc[DOXG4_AM_selected,5] <- c(rep("#005eff",5),rep("firebrick1",5))
DOXG4_AM_volc[DOXG4_AM_selected,6] <- 1

DOXG4_AM_volcplot <- DOXG4_AM_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="red",linetype="solid",alpha=0.4) +
  scale_y_continuous(limits = c(-0.12, 13), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-10, 10), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = DOXG4_AM_volc[DOXG4_AM_selected,], aes(label = str_wrap(name, 23)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(1,-0.5,0.5,0,1,7,4.3,0.5,4,0.5),nudge_x=c(-4,-3,-4,-2,-1,1,5,3,2,0))
DOXG4_AM_volcplot

#################### DOXG4_AF ###

# import DESeq2 results for comparison of interest
DOXG4_AF_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_DOXG4_H2OG4_AF.csv")

# remove spurious genes with adj pvals < 0.05
DOXG4_AF_volc_diff <- setdiff(rownames(res_DOXG4_H2OG4_AF_sig),rownames(chopping_block_counts_DOXG4_AF_v_H2OG4_AF))
DOXG4_AF_volc <- DOXG4_AF_volc[!(DOXG4_AF_volc$X %in% DOXG4_AF_volc_diff),]

# remove unnecessary columns
DOXG4_AF_volc <- DOXG4_AF_volc[,c(1,3,7)]

# order by padj (low to high)
DOXG4_AF_volc <- DOXG4_AF_volc[order(DOXG4_AF_volc$padj),]

# add new empty column for locus names:
DOXG4_AF_volc$name <- rep("NA",nrow(DOXG4_AF_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(DOXG4_AF_volc[DOXG4_AF_volc$padj < 0.05 & !(is.na(DOXG4_AF_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_DOXG4_AF.csv")

# import gene IDs
GeneIDs_DOXG4_AF <- read.csv("GOI_DOXG4_AF.csv")
DOXG4_AF_volc[1:dim(DOXG4_AF_volc[DOXG4_AF_volc$padj < 0.05 & !(is.na(DOXG4_AF_volc$padj)),])[1],4] <- GeneIDs_DOXG4_AF$name

# make locus IDs rownames so we can index by them
rownames(DOXG4_AF_volc) <- DOXG4_AF_volc$X

# set color and transparency values for genes
DOXG4_AF_volc$color <- "gray40"
DOXG4_AF_volc$alpha <- 0.5

DOXG4_AF_selected <- c("LOC124421152","LOC111688178","LOC111680476")

DOXG4_AF_volc[DOXG4_AF_selected,5] <- c(rep("firebrick1",3))
DOXG4_AF_volc[DOXG4_AF_selected,6] <- 1

# abbreviate as needed

DOXG4_AF_volcplot <- DOXG4_AF_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(-0.07, 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-10, 10), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = DOXG4_AF_volc[DOXG4_AF_selected,], aes(label = str_wrap(name, 25)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(3,2,1),nudge_x=c(2,4,2))

DOXG4_AF_volcplot

################### ATC G3 L


# import DESeq2 results for comparison of interest
ATCG3_L_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG3_H2OG3_L.csv")

# remove spurious genes with adj pvals < 0.05
ATCG3_L_volc_diff <- setdiff(rownames(res_ATCG3_H2OG3_L_sig),rownames(chopping_block_counts_ATCG3_L_v_H2OG3_L))
ATCG3_L_volc <- ATCG3_L_volc[!(ATCG3_L_volc$X %in% ATCG3_L_volc_diff),]

# remove unnecessary columns
ATCG3_L_volc <- ATCG3_L_volc[,c(1,3,7)]

# order by padj (low to high)
ATCG3_L_volc <- ATCG3_L_volc[order(ATCG3_L_volc$padj),]

# add new empty column for locus names:
ATCG3_L_volc$name <- rep("NA",nrow(ATCG3_L_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(ATCG3_L_volc[ATCG3_L_volc$padj < 0.05 & !(is.na(ATCG3_L_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_ATCG3_L.csv")

# import gene IDs
GeneIDs_ATCG3_L <- read.csv("GOI_ATCG3_L.csv")
ATCG3_L_volc[1:dim(ATCG3_L_volc[ATCG3_L_volc$padj < 0.05 & !(is.na(ATCG3_L_volc$padj)),])[1],4] <- GeneIDs_ATCG3_L$name

# make locus IDs rownames so we can index by them
rownames(ATCG3_L_volc) <- ATCG3_L_volc$X

# set color and transparency values for genes
ATCG3_L_volc$color <- "gray40"
ATCG3_L_volc$alpha <- 0.5

ATCG3_L_selected <- c("LOC111681823","LOC111687047","LOC111678842","LOC124418485","LOC111678840","LOC111682079","LOC111681373","LOC111683436","LOC111677781","LOC111681336")

ATCG3_L_volc[ATCG3_L_selected,5] <- c(rep("#005eff",5),rep("firebrick1",5))
ATCG3_L_volc[ATCG3_L_selected,6] <- 1

ATCG3_L_volcplot <- ATCG3_L_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(-0.3, 40), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-10, 10), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = ATCG3_L_volc[ATCG3_L_selected,], aes(label = str_wrap(name, 25)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(3,14,8,1,20,23,1.2,6,6.5,4),nudge_x=c(-4,-3,-4.8,-5,-0.5,-0.5,4,3.5,3,1))

ATCG3_L_volcplot

#################### ATCG3_AM ###

# import DESeq2 results for comparison of interest
ATCG3_AM_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG3_H2OG3_AM.csv")

# remove spurious genes with adj pvals < 0.05
ATCG3_AM_volc_diff <- setdiff(rownames(res_ATCG3_H2OG3_AM_sig),rownames(chopping_block_counts_ATCG3_AM_v_H2OG3_AM))
ATCG3_AM_volc <- ATCG3_AM_volc[!(ATCG3_AM_volc$X %in% ATCG3_AM_volc_diff),]

# remove unnecessary columns
ATCG3_AM_volc <- ATCG3_AM_volc[,c(1,3,7)]

# order by padj (low to high)
ATCG3_AM_volc <- ATCG3_AM_volc[order(ATCG3_AM_volc$padj),]

# add new empty column for locus names:
ATCG3_AM_volc$name <- rep("NA",nrow(ATCG3_AM_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(ATCG3_AM_volc[ATCG3_AM_volc$padj < 0.05 & !(is.na(ATCG3_AM_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_ATCG3_AM.csv")

# import gene IDs
GeneIDs_ATCG3_AM <- read.csv("GOI_ATCG3_AM.csv")
ATCG3_AM_volc[1:dim(ATCG3_AM_volc[ATCG3_AM_volc$padj < 0.05 & !(is.na(ATCG3_AM_volc$padj)),])[1],4] <- GeneIDs_ATCG3_AM$name

# make locus IDs rownames so we can index by them
rownames(ATCG3_AM_volc) <- ATCG3_AM_volc$X

# set color and transparency values for genes
ATCG3_AM_volc$color <- "gray40"
ATCG3_AM_volc$alpha <- 0.5

ATCG3_AM_selected <- c("LOC111677764","LOC111687136","LOC111675758","LOC111686359","LOC111677302","LOC111681060","LOC111687045","LOC111685676","LOC111688465","LOC111688454")
ATCG3_AM_volc[ATCG3_AM_selected,5] <- c(rep("#005eff",5),rep("firebrick1",5))
ATCG3_AM_volc[ATCG3_AM_selected,6] <- 1

# abbreviate as needed
#ATCG3_AM_volc[ATCG3_AM_selected,4][c(1,4,5,6,9,10)] <- c("eukaryotic translation init. factor 6*","RNA Pol. I-specific transcription init. factor RRN3*","LOC111677744","LOC111681969","LOC111676913","LOC124419431, transcript var. X2*")

ATCG3_AM_volcplot <- ATCG3_AM_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(-0.15, 20), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-8, 8), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = ATCG3_AM_volc[ATCG3_AM_selected,], aes(label = str_wrap(name, 25)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(4.8,1,1,3,5.3,0,0.8,-1,0.7,1),nudge_x=c(-5,-1,-5,-5,-1.4,2,1.5,1.5,1,1))

ATCG3_AM_volcplot

#################### ATCG3_AF ###

# import DESeq2 results for comparison of interest
ATCG3_AF_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG3_H2OG3_AF.csv")

# remove spurious genes with adj pvals < 0.05
ATCG3_AF_volc_diff <- setdiff(rownames(res_ATCG3_H2OG3_AF_sig),rownames(chopping_block_counts_ATCG3_AF_v_H2OG3_AF))
ATCG3_AF_volc <- ATCG3_AF_volc[!(ATCG3_AF_volc$X %in% ATCG3_AF_volc_diff),]

# remove unnecessary columns
ATCG3_AF_volc <- ATCG3_AF_volc[,c(1,3,7)]

# order by padj (low to high)
ATCG3_AF_volc <- ATCG3_AF_volc[order(ATCG3_AF_volc$padj),]

# add new empty column for locus names:
ATCG3_AF_volc$name <- rep("NA",nrow(ATCG3_AF_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(ATCG3_AF_volc[ATCG3_AF_volc$padj < 0.05 & !(is.na(ATCG3_AF_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_ATCG3_AF.csv")

# import gene IDs
GeneIDs_ATCG3_AF <- read.csv("GOI_ATCG3_AF.csv")
ATCG3_AF_volc[1:dim(ATCG3_AF_volc[ATCG3_AF_volc$padj < 0.05 & !(is.na(ATCG3_AF_volc$padj)),])[1],4] <- GeneIDs_ATCG3_AF$name

# make locus IDs rownames so we can index by them
rownames(ATCG3_AF_volc) <- ATCG3_AF_volc$X

# set color and transparency values for genes
ATCG3_AF_volc$color <- "gray40"
ATCG3_AF_volc$alpha <- 0.5

ATCG3_AF_selected <- c("LOC111682291","LOC111689988","LOC111690347","LOC111682435","LOC111680128")

ATCG3_AF_volc[ATCG3_AF_selected,5] <- c("#005eff",rep("firebrick1",4))
ATCG3_AF_volc[ATCG3_AF_selected,6] <- 1

ATCG3_AF_volcplot <- ATCG3_AF_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(-0.06, 8), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-12, 12), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = ATCG3_AF_volc[ATCG3_AF_selected,], aes(label = str_wrap(name, 25)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(0.5,1,-0.2,0.5,0.7),nudge_x=c(-1,2,4,4,-0.7))

ATCG3_AF_volcplot

############################## ATC G4 L

# import DESeq2 results for comparison of interest
ATCG4_L_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG4_H2OG4_L.csv")

# remove spurious genes with adj pvals < 0.05
ATCG4_L_volc_diff <- setdiff(rownames(res_ATCG4_H2OG4_L_sig),rownames(chopping_block_counts_ATCG4_L_v_H2OG4_L))
ATCG4_L_volc <- ATCG4_L_volc[!(ATCG4_L_volc$X %in% ATCG4_L_volc_diff),]

# remove unnecessary columns
ATCG4_L_volc <- ATCG4_L_volc[,c(1,3,7)]

# order by padj (low to high)
ATCG4_L_volc <- ATCG4_L_volc[order(ATCG4_L_volc$padj),]

# add new empty column for locus names:
ATCG4_L_volc$name <- rep("NA",nrow(ATCG4_L_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(ATCG4_L_volc[ATCG4_L_volc$padj < 0.05 & !(is.na(ATCG4_L_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_ATCG4_L.csv")

# import gene IDs
GeneIDs_ATCG4_L <- read.csv("GOI_ATCG4_L.csv")
ATCG4_L_volc[1:dim(ATCG4_L_volc[ATCG4_L_volc$padj < 0.05 & !(is.na(ATCG4_L_volc$padj)),])[1],4] <- GeneIDs_ATCG4_L$name

# make locus IDs rownames so we can index by them
rownames(ATCG4_L_volc) <- ATCG4_L_volc$X

# set color and transparency values for genes
ATCG4_L_volc$color <- "gray40"
ATCG4_L_volc$alpha <- 0.5

ATCG4_L_volcplot <- ATCG4_L_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(-0.07, 8), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-8, 8), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") 
#+ geom_label_repel(data = ATCG4_L_volc[ATCG4_L_selected,], aes(label = str_wrap(name, 25)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(0.5,0.5,0.5,0.5,-0.1,0,0.9,1.3,0.5,0.6),nudge_x=c(0,0,0,-0.7,-2,2.2,-0.1,-1.2,1.3,1))

ATCG4_L_volcplot

#################### ATCG4_AM ###

# import DESeq2 results for comparison of interest
ATCG4_AM_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG4_H2OG4_AM.csv")

# remove spurious genes with adj pvals < 0.05
ATCG4_AM_volc_diff <- setdiff(rownames(res_ATCG4_H2OG4_AM_sig),rownames(chopping_block_counts_ATCG4_AM_v_H2OG4_AM))
ATCG4_AM_volc <- ATCG4_AM_volc[!(ATCG4_AM_volc$X %in% ATCG4_AM_volc_diff),]

# remove unnecessary columns
ATCG4_AM_volc <- ATCG4_AM_volc[,c(1,3,7)]

# order by padj (low to high)
ATCG4_AM_volc <- ATCG4_AM_volc[order(ATCG4_AM_volc$padj),]

# add new empty column for locus names:
ATCG4_AM_volc$name <- rep("NA",nrow(ATCG4_AM_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(ATCG4_AM_volc[ATCG4_AM_volc$padj < 0.05 & !(is.na(ATCG4_AM_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_ATCG4_AM.csv")

# import gene IDs
GeneIDs_ATCG4_AM <- read.csv("GOI_ATCG4_AM.csv")
ATCG4_AM_volc[1:dim(ATCG4_AM_volc[ATCG4_AM_volc$padj < 0.05 & !(is.na(ATCG4_AM_volc$padj)),])[1],4] <- GeneIDs_ATCG4_AM$name

# make locus IDs rownames so we can index by them
rownames(ATCG4_AM_volc) <- ATCG4_AM_volc$X

# set color and transparency values for genes
ATCG4_AM_volc$color <- "gray40"
ATCG4_AM_volc$alpha <- 0.5

ATCG4_AM_selected <- c("LOC124419621","LOC111689926","LOC111676477","LOC111686208","LOC111682763","LOC111679313","LOC124418474","LOC111681467","LOC111674639","LOC124419761")

ATCG4_AM_volc[ATCG4_AM_selected,5] <- c(rep("#005eff",5),rep("firebrick1",5))
ATCG4_AM_volc[ATCG4_AM_selected,6] <- 1

ATCG4_AM_volcplot <- ATCG4_AM_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(-0.07, 9), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-6, 6), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = ATCG4_AM_volc[ATCG4_AM_selected,], aes(label = str_wrap(name, 31)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(0.5,1,0.4,1.7,1.1,1.8,0.8,1.2,0.3,1),nudge_x=c(-0.4,-0.4,-2,-3,-3,3,3,2.2,3,0.5))
ATCG4_AM_volcplot

#################### ATCG4_AF ###

# import DESeq2 results for comparison of interest
ATCG4_AF_volc <- read.csv("~/Desktop/RNAseq/Data/DEGs/Unfiltered/DEGs_ATCG4_H2OG4_AF.csv")

# remove spurious genes with adj pvals < 0.05
ATCG4_AF_volc_diff <- setdiff(rownames(res_ATCG4_H2OG4_AF_sig),rownames(chopping_block_counts_ATCG4_AF_v_H2OG4_AF))
ATCG4_AF_volc <- ATCG4_AF_volc[!(ATCG4_AF_volc$X %in% ATCG4_AF_volc_diff),]

# remove unnecessary columns
ATCG4_AF_volc <- ATCG4_AF_volc[,c(1,3,7)]

# order by padj (low to high)
ATCG4_AF_volc <- ATCG4_AF_volc[order(ATCG4_AF_volc$padj),]

# add new empty column for locus names:
ATCG4_AF_volc$name <- rep("NA",nrow(ATCG4_AF_volc))

# export the genes with p adj < 0.05 as a .csv
write.csv(ATCG4_AF_volc[ATCG4_AF_volc$padj < 0.05 & !(is.na(ATCG4_AF_volc$padj)),],"~/Desktop/RNAseq/Data/DEGs/Unfiltered/GOI/GOI_ATCG4_AF.csv")

# import gene IDs
GeneIDs_ATCG4_AF <- read.csv("GOI_ATCG4_AF.csv")
ATCG4_AF_volc[1:dim(ATCG4_AF_volc[ATCG4_AF_volc$padj < 0.05 & !(is.na(ATCG4_AF_volc$padj)),])[1],4] <- GeneIDs_ATCG4_AF$name

# make locus IDs rownames so we can index by them
rownames(ATCG4_AF_volc) <- ATCG4_AF_volc$X

# set color and transparency values for genes
ATCG4_AF_volc$color <- "gray40"
ATCG4_AF_volc$alpha <- 0.5

ATCG4_AF_selected <- c("LOC111689384","LOC111676136","LOC111674796","LOC111682256","LOC111677994","LOC124421090","LOC111681708","LOC111680380","LOC111688492","LOC111685972")

ATCG4_AF_volc[ATCG4_AF_selected,5] <- c(rep("#005eff",5),rep("firebrick1",5))
ATCG4_AF_volc[ATCG4_AF_selected,6] <- 1

ATCG4_AF_volcplot <- ATCG4_AF_volc %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),color=color,alpha=alpha)) +
  scale_color_identity() +
  scale_alpha_identity() +
  geom_point(size=2) +  
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.background=element_blank(), panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.text.y = element_text(size=12,colour = "black",face="bold"),axis.title.y = element_text(size=12,color="black",face="bold"),axis.title.x = element_text(size=12,color="black",face="bold"),axis.text.x = element_text(size=12,color="black",face="bold")) +
  geom_vline(xintercept = 0,color="black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.05),color="orange",linetype="solid",alpha=0.8) +
  scale_y_continuous(limits = c(-0.06, 8), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-13, 13), expand = c(0, 0)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_label_repel(data = ATCG4_AF_volc[ATCG4_AF_selected,], aes(label = str_wrap(name, 30)),force = 0,min.segment.length = 0,label.size = NA,fontface="bold",fill=NA,nudge_y=c(0.4,1.8,2,0.3,0.5,1.3,0.5,0.4,1.2,0.7),nudge_x=c(-3.5,-5.5,-2,-6,-2,0.5,7,1.5,6,2))
ATCG4_AF_volcplot
# Export plots

DOXG3_L_volcplot
ATCG3_L_volcplot
DOXG4_L_volcplot
ATCG4_L_volcplot

DOXG3_AM_volcplot
ATCG3_AM_volcplot
DOXG3_AF_volcplot
ATCG3_AF_volcplot

DOXG4_AM_volcplot
ATCG4_AM_volcplot
DOXG4_AF_volcplot
ATCG4_AF_volcplot


plot_grid(DOXG3_L_volcplot,DOXG3_AM_volcplot,DOXG3_AF_volcplot,ATCG3_L_volcplot,ATCG3_AM_volcplot,ATCG3_AF_volcplot,DOXG4_L_volcplot,DOXG4_AM_volcplot,DOXG4_AF_volcplot,ATCG4_L_volcplot,ATCG4_AM_volcplot,ATCG4_AF_volcplot,ncol=3,align="vh")

plot_grid(DOXG3_L_volcplot,DOXG4_L_volcplot,ATCG3_L_volcplot,ATCG4_L_volcplot,ncol=4,align="vh")
plot_grid(DOXG3_AM_volcplot,DOXG3_AF_volcplot,DOXG4_AM_volcplot,DOXG4_AF_volcplot,ATCG3_AM_volcplot,ATCG3_AF_volcplot,ATCG4_AM_volcplot,ATCG4_AF_volcplot,ncol=4,align="vh")

#----------------------------------------------------------------------
# PCA (Figure S2) showing clustering of RNAseq samples
#----------------------------------------------------------------------

# VST-transform the data for PCA
vsd <- vst(dds_full, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))

PCA_data <- plotPCA(vsd, intgroup=c("sex", "stage","fullgroup"), returnData=TRUE)
PCA_data$group <- as.character(PCA_data$group)
PCA_data$group <- gsub("(.*):(.*)", "\\2 \\1", PCA_data$group)
percentVar_stage <- round(100 * attr(PCA_data , "percentVar"))

PCA_all_plot <-
  ggplot(PCA_data, aes(PC1, PC2, color=fullgroup)) +
  geom_point(size=2.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position="right",axis.text.y = element_text(size=10),axis.text.x = element_text(size=10),legend.text = element_text(size=8), legend.title = element_text(size=10)) +
  xlab(paste0("PC1: ",percentVar_stage[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_stage[2],"% variance")) +
  scale_color_manual(values = c("A_H2O_G4_F"="steelblue1", "A_DOX_G4_F"="hotpink", "A_ATC_G4_F"="green")) +
  #scale_shape_manual(values = c("Adult Male"=1, "Adult Female"=1, "Larva Male"=5, "Larva Female"=5)) +
  coord_fixed(ylim=c(-75,75), xlim=c(-75,75))

PCA_all_plot <- PCA_all_plot + 
  theme(legend.title=element_blank())

PCA_all_plot

#----------------------------------------------------------------------
# Barplot showing DEGs across all samples (Figure 3)
#----------------------------------------------------------------------

# make the grouped barplot (4 main comparisons)

comparison=c(rep("G3 DOX\nvs.\nG3 H2O", 2), rep("G3 ATC\nvs.\nG3 H2O",2), rep("G4 DOX\nvs.\nG4 H2O",2), rep("G4 ATC\nvs.\nG4 H2O",2))
change=rep(c("downregulated","upregulated"),4)
numbers=c(length(DOXG3_downreg),length(DOXG3_upreg),length(ATCG3_downreg),length(ATCG3_upreg),length(DOXG4_downreg),length(DOXG4_upreg),length(ATCG4_downreg),length(ATCG4_upreg))

data <- data.frame(comparison, change, numbers)

fct_relevel(factor(data$comparison),"G3 DOX\nvs.\nG3 H2O","G3 ATC\nvs.\nG3 H2O","G4 DOX\nvs.\nG4 H2O","G4 ATC\nvs.\nG4 H2O")

baseplot <- ggplot(data, aes(fill=change, y=numbers, x=fct_relevel(factor(comparison),"G3 DOX\nvs.\nG3 H2O","G3 ATC\nvs.\nG3 H2O","G4 DOX\nvs.\nG4 H2O","G4 ATC\nvs.\nG4 H2O"))) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label=numbers),position=position_dodge(0.9),vjust=-0.5,fontface="bold") +
  scale_fill_manual(values=c("#31a4a9","#f9c132")) +
  labs(y="Number of genes",x="Comparisons") +
  scale_y_continuous(limits = c(0, 1200), expand = c(0, 0)) +
  theme(panel.background = element_blank(), legend.text = element_text(size=12,colour = "black",face="bold"), axis.line = element_line(size = 0.5, colour = "black"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.x = element_text(size=15,face="bold",margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(size=15,face="bold",margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.text.y = element_text(size=12,colour = "black",face="bold"), axis.text.x = element_text(size=10,color="black",face="bold"), strip.text.x= element_text(size=15)) 

baseplot + theme(legend.position="inside",legend.justification = c("right", "top"),legend.title=element_blank())

#----------------------------------------------------------------------
# Euler plots: DOX vs. ATC in male comp groups
#----------------------------------------------------------------------

all_loci_upreg <- union(ATCG4_upreg, union(DOXG3_upreg, union(DOXG4_upreg, ATCG3_upreg)))
all_loci_upreg

up_DOXG3_H2OG3_LF <- as.numeric(all_loci_upreg %in% rownames(upreg_res_DOXG3_H2OG3_LF_ordered))
up_DOXG3_H2OG3_L <- as.numeric(all_loci_upreg %in% rownames(upreg_res_DOXG3_H2OG3_L_ordered))
up_DOXG3_H2OG3_AF <- as.numeric(all_loci_upreg %in% rownames(upreg_res_DOXG3_H2OG3_AF_ordered))
up_DOXG3_H2OG3_AM <- as.numeric(all_loci_upreg %in% rownames(upreg_res_DOXG3_H2OG3_AM_ordered))
up_ATCG3_H2OG3_LF <- as.numeric(all_loci_upreg %in% rownames(upreg_res_ATCG3_H2OG3_LF_ordered))
up_ATCG3_H2OG3_L <- as.numeric(all_loci_upreg %in% rownames(upreg_res_ATCG3_H2OG3_L_ordered))
up_ATCG3_H2OG3_AF <- as.numeric(all_loci_upreg %in% rownames(upreg_res_ATCG3_H2OG3_AF_ordered))
up_ATCG3_H2OG3_AM <- as.numeric(all_loci_upreg %in% rownames(upreg_res_ATCG3_H2OG3_AM_ordered))
up_DOXG4_H2OG4_LF <- as.numeric(all_loci_upreg %in% rownames(upreg_res_DOXG4_H2OG4_LF_ordered))
up_DOXG4_H2OG4_L <- as.numeric(all_loci_upreg %in% rownames(upreg_res_DOXG4_H2OG4_L_ordered))
up_DOXG4_H2OG4_AF <- as.numeric(all_loci_upreg %in% rownames(upreg_res_DOXG4_H2OG4_AF_ordered))
up_DOXG4_H2OG4_AM <- as.numeric(all_loci_upreg %in% rownames(upreg_res_DOXG4_H2OG4_AM_ordered))
up_ATCG4_H2OG4_LF <- as.numeric(all_loci_upreg %in% rownames(upreg_res_ATCG4_H2OG4_LF_ordered))
up_ATCG4_H2OG4_L <- as.numeric(all_loci_upreg %in% rownames(upreg_res_ATCG4_H2OG4_L_ordered))
up_ATCG4_H2OG4_AF <- as.numeric(all_loci_upreg %in% rownames(upreg_res_ATCG4_H2OG4_AF_ordered))
up_ATCG4_H2OG4_AM <- as.numeric(all_loci_upreg %in% rownames(upreg_res_ATCG4_H2OG4_AM_ordered))

otumat = matrix(c(up_DOXG3_H2OG3_LF,up_DOXG3_H2OG3_L,up_DOXG3_H2OG3_AF,up_DOXG3_H2OG3_AM,up_ATCG3_H2OG3_LF,up_ATCG3_H2OG3_L,up_ATCG3_H2OG3_AF,up_ATCG3_H2OG3_AM,up_DOXG4_H2OG4_LF,up_DOXG4_H2OG4_L,up_DOXG4_H2OG4_AF,up_DOXG4_H2OG4_AM,up_ATCG4_H2OG4_LF,up_ATCG4_H2OG4_L,up_ATCG4_H2OG4_AF,up_ATCG4_H2OG4_AM), nrow = length(all_loci_upreg), ncol = 16)
otumat

rownames(otumat) <- all_loci_upreg
colnames(otumat) <- c("up_DOXG3_H2OG3_LF","up_DOXG3_H2OG3_L","up_DOXG3_H2OG3_AF","up_DOXG3_H2OG3_AM","up_ATCG3_H2OG3_LF","up_ATCG3_H2OG3_L","up_ATCG3_H2OG3_AF","up_ATCG3_H2OG3_AM","up_DOXG4_H2OG4_LF","up_DOXG4_H2OG4_L","up_DOXG4_H2OG4_AF","up_DOXG4_H2OG4_AM","up_ATCG4_H2OG4_LF","up_ATCG4_H2OG4_L","up_ATCG4_H2OG4_AF","up_ATCG4_H2OG4_AM")
otumat

taxmat = matrix(sample(letters, length(all_loci_upreg) * 7, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU

physeq = phyloseq(OTU, TAX)

gendiet <- c(rep("DOXG3_H2OG3",4),rep("ATCG3_H2OG3",4),rep("DOXG4_H2OG4",4),rep("ATCG4_H2OG4",4))
sexstage <- rep(c("LF","L","AF","AM"),4)
stage <- rep(c("L", "L","A", "A"),4)
sampledata = sample_data(data.frame(gendiet=gendiet,sexstage=sexstage,stage=stage))
rownames(sampledata) <- sample_names(physeq)

physeq_upreg = merge_phyloseq(physeq, sampledata)

# 12 genes shared
upreg_list <- ps_euler(physeq_upreg, "stage", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","yellow"), main = "G3 Larvae", shape="ellipse",plot=FALSE)
upreg_plot <- ps_euler(physeq_upreg, "stage", quantities = list(type=c("counts")), fraction = 0, labels = list(cex=0), col = "black", fill = c("darkseagreen1", "darkseagreen", "thistle1","plum"), main = "Downregulated genes", shape="ellipse",plot=TRUE)
ps_DOX_G3_upreg <- phyloseq::subset_samples(physeq_upreg, gendiet=="DOXG3_H2OG3")
ps_DOX_G3_upreg_list <- ps_euler(ps_DOX_G3_upreg, "stage", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","yellow"), main = "G3 Larvae", shape="ellipse",plot=FALSE)
ps_DOX_G3_upreg_plot <- ps_euler(ps_DOX_G3_upreg, "stage", quantities = list(type=c("counts")), fraction = 0, labels = list(cex=0), col = "black", fill = c("darkseagreen1", "darkseagreen", "thistle1","plum"), main = "Downregulated genes", shape="ellipse",plot=TRUE)

ps_DOX_G4_upreg <- phyloseq::subset_samples(physeq_upreg, gendiet=="DOXG4_H2OG4")
ps_DOX_G4_upreg_list <- ps_euler(ps_DOX_G4_upreg, "stage", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","yellow"), main = "G3 Larvae", shape="ellipse",plot=FALSE)
ps_DOX_G4_upreg_plot <- ps_euler(ps_DOX_G4_upreg, "stage", quantities = list(type=c("counts")), fraction = 0, labels = list(cex=0), col = "black", fill = c("darkseagreen1", "darkseagreen", "thistle1","plum"), main = "Downregulated genes", shape="ellipse",plot=TRUE)
ps_DOX_G4_upreg_list
ps_DOX_G4_upreg_plot

ps_ATC_G4_upreg <- phyloseq::subset_samples(physeq_upreg, gendiet=="ATCG4_H2OG4")
ps_ATC_G4_upreg_list <- ps_euler(ps_ATC_G4_upreg, "stage", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","yellow"), main = "G3 Larvae", shape="ellipse",plot=FALSE)
ps_ATC_G4_upreg_plot <- ps_euler(ps_ATC_G4_upreg, "stage", quantities = list(type=c("counts")), fraction = 0, labels = list(cex=0), col = "black", fill = c("darkseagreen1", "darkseagreen", "thistle1","plum"), main = "Downregulated genes", shape="ellipse",plot=TRUE)
ps_ATC_G4_upreg_list
ps_ATC_G4_upreg_plot

# 1 shared upreg gene
ps_ATC_G3_upreg <- phyloseq::subset_samples(physeq_upreg, gendiet=="ATCG3_H2OG3")
ps_ATC_G3_upreg_list <- ps_euler(ps_ATC_G3_upreg, "stage", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","yellow"), main = "G3 Larvae", shape="ellipse",plot=FALSE)
ps_ATC_G3_upreg_plot <- ps_euler(ps_ATC_G3_upreg, "stage", quantities = list(type=c("counts")), fraction = 0, labels = list(cex=0), col = "black", fill = c("darkseagreen1", "darkseagreen", "thistle1","plum"), main = "Downregulated genes", shape="ellipse",plot=TRUE)
ps_ATC_G3_upreg_list
ps_ATC_G3_upreg_plot

ps_AM_upreg <- phyloseq::subset_samples(physeq_upreg, sexstage=="AM")
ps_ATCG3_H2OG3_upreg <- phyloseq::subset_samples(physeq_upreg, gendiet=="ATCG3_H2OG3")
ps_DOXG4_H2OG4_upreg <- phyloseq::subset_samples(physeq_upreg, gendiet=="DOXG4_H2OG4")
ps_ATCG4_H2OG4_upreg <- phyloseq::subset_samples(physeq_upreg, gendiet=="ATCG4_H2OG4")

all_loci_downreg <- union(ATCG4_downreg, union(DOXG3_downreg, union(DOXG4_downreg, ATCG3_downreg)))
all_loci_downreg

down_DOXG3_H2OG3_LF <- as.numeric(all_loci_downreg %in% rownames(downreg_res_DOXG3_H2OG3_LF_ordered))
down_DOXG3_H2OG3_L <- as.numeric(all_loci_downreg %in% rownames(downreg_res_DOXG3_H2OG3_L_ordered))
down_DOXG3_H2OG3_AF <- as.numeric(all_loci_downreg %in% rownames(downreg_res_DOXG3_H2OG3_AF_ordered))
down_DOXG3_H2OG3_AM <- as.numeric(all_loci_downreg %in% rownames(downreg_res_DOXG3_H2OG3_AM_ordered))
down_ATCG3_H2OG3_LF <- as.numeric(all_loci_downreg %in% rownames(downreg_res_ATCG3_H2OG3_LF_ordered))
down_ATCG3_H2OG3_L <- as.numeric(all_loci_downreg %in% rownames(downreg_res_ATCG3_H2OG3_L_ordered))
down_ATCG3_H2OG3_AF <- as.numeric(all_loci_downreg %in% rownames(downreg_res_ATCG3_H2OG3_AF_ordered))
down_ATCG3_H2OG3_AM <- as.numeric(all_loci_downreg %in% rownames(downreg_res_ATCG3_H2OG3_AM_ordered))
down_DOXG4_H2OG4_LF <- as.numeric(all_loci_downreg %in% rownames(downreg_res_DOXG4_H2OG4_LF_ordered))
down_DOXG4_H2OG4_L <- as.numeric(all_loci_downreg %in% rownames(downreg_res_DOXG4_H2OG4_L_ordered))
down_DOXG4_H2OG4_AF <- as.numeric(all_loci_downreg %in% rownames(downreg_res_DOXG4_H2OG4_AF_ordered))
down_DOXG4_H2OG4_AM <- as.numeric(all_loci_downreg %in% rownames(downreg_res_DOXG4_H2OG4_AM_ordered))
down_ATCG4_H2OG4_LF <- as.numeric(all_loci_downreg %in% rownames(downreg_res_ATCG4_H2OG4_LF_ordered))
down_ATCG4_H2OG4_L <- as.numeric(all_loci_downreg %in% rownames(downreg_res_ATCG4_H2OG4_L_ordered))
down_ATCG4_H2OG4_AF <- as.numeric(all_loci_downreg %in% rownames(downreg_res_ATCG4_H2OG4_AF_ordered))
down_ATCG4_H2OG4_AM <- as.numeric(all_loci_downreg %in% rownames(downreg_res_ATCG4_H2OG4_AM_ordered))

otumat = matrix(c(down_DOXG3_H2OG3_LF,down_DOXG3_H2OG3_L,down_DOXG3_H2OG3_AF,down_DOXG3_H2OG3_AM,down_ATCG3_H2OG3_LF,down_ATCG3_H2OG3_L,down_ATCG3_H2OG3_AF,down_ATCG3_H2OG3_AM,down_DOXG4_H2OG4_LF,down_DOXG4_H2OG4_L,down_DOXG4_H2OG4_AF,down_DOXG4_H2OG4_AM,down_ATCG4_H2OG4_LF,down_ATCG4_H2OG4_L,down_ATCG4_H2OG4_AF,down_ATCG4_H2OG4_AM), nrow = length(all_loci_downreg), ncol = 16)
otumat

rownames(otumat) <- all_loci_downreg
colnames(otumat) <- c("down_DOXG3_H2OG3_LF","down_DOXG3_H2OG3_L","down_DOXG3_H2OG3_AF","down_DOXG3_H2OG3_AM","down_ATCG3_H2OG3_LF","down_ATCG3_H2OG3_L","down_ATCG3_H2OG3_AF","down_ATCG3_H2OG3_AM","down_DOXG4_H2OG4_LF","down_DOXG4_H2OG4_L","down_DOXG4_H2OG4_AF","down_DOXG4_H2OG4_AM","down_ATCG4_H2OG4_LF","down_ATCG4_H2OG4_L","down_ATCG4_H2OG4_AF","down_ATCG4_H2OG4_AM")
otumat

taxmat = matrix(sample(letters, length(all_loci_downreg) * 7, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU

physeq = phyloseq(OTU, TAX)

gendiet <- c(rep("DOXG3_H2OG3",4),rep("ATCG3_H2OG3",4),rep("DOXG4_H2OG4",4),rep("ATCG4_H2OG4",4))
sexstage <- rep(c("LF","L","AF","AM"),4)
stage <- rep(c("L", "L","A", "A"),4)
sampledata = sample_data(data.frame(gendiet=gendiet,sexstage=sexstage,stage=stage))
rownames(sampledata) <- sample_names(physeq)

physeq_downreg = merge_phyloseq(physeq, sampledata)

ps_ATC_G3_downreg <- phyloseq::subset_samples(physeq_downreg, gendiet=="ATCG3_H2OG3")
ps_ATC_G3_downreg_list <- ps_euler(ps_ATC_G3_downreg, "stage", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","yellow"), main = "G3 Larvae", shape="ellipse",plot=FALSE)
ps_ATC_G3_downreg_plot <- ps_euler(ps_ATC_G3_downreg, "stage", quantities = list(type=c("counts")), fraction = 0, labels = list(cex=0), col = "black", fill = c("darkseagreen1", "darkseagreen", "thistle1","plum"), main = "Downregulated genes", shape="ellipse",plot=TRUE)
ps_ATC_G3_downreg_list
ps_ATC_G3_downreg_plot

# 3 genes shared
ps_DOX_G3_downreg <- phyloseq::subset_samples(physeq_downreg, gendiet=="DOXG3_H2OG3")
ps_DOX_G3_downreg_list <- ps_euler(ps_DOX_G3_downreg, "stage", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","yellow"), main = "G3 Larvae", shape="ellipse",plot=FALSE)
ps_DOX_G3_downreg_plot <- ps_euler(ps_DOX_G3_downreg, "stage", quantities = list(type=c("counts")), fraction = 0, labels = list(cex=0), col = "black", fill = c("darkseagreen1", "darkseagreen", "thistle1","plum"), main = "Downregulated genes", shape="ellipse",plot=TRUE)
ps_DOX_G3_downreg_list
ps_DOX_G3_downreg_plot

ps_AM_downreg <- phyloseq::subset_samples(physeq_downreg, sexstage=="AM")
ps_L_downreg <- phyloseq::subset_samples(physeq_downreg, sexstage=="L")
ps_L_downreg <- phyloseq::subset_samples(physeq_downreg, stage=="L")
ps_A_downreg <- phyloseq::subset_samples(physeq_downreg, stage=="A")

# Code for plotting Euler diagrams of shared and unique downreg genes in adult males

downreg_plot <- ps_euler(ps_AM_downreg, "gendiet", quantities = list(type=c("counts")), fraction = 1, labels = list(cex=0), col = "black", fill = c("darkseagreen1", "darkseagreen", "thistle1","plum"), main = "Downregulated genes", shape="ellipse",plot=TRUE)
downreg_list <- ps_euler(ps_AM_downreg, "gendiet", fraction = 0.5, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","yellow"), main = "G3 Larvae", shape="ellipse",plot=FALSE)

downreg_list <- ps_euler(physeq_downreg, "stage", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","yellow"), main = "G3 Larvae", shape="ellipse",plot=FALSE)
downreg_plot <- ps_euler(physeq_downreg, "stage", quantities = list(type=c("counts")), fraction = 0, labels = list(cex=0), col = "black", fill = c("darkseagreen1", "darkseagreen", "thistle1","plum"), main = "Downregulated genes", shape="ellipse",plot=TRUE)

# Plot
downreg_plot
# List of genes
downreg_list

## Code for Euler diagrams of upregulated genes in adult male samples
upreg_plot <- ps_euler(ps_AM_upreg, "gendiet", quantities = list(type=c("counts")), fraction = 1, labels = list(cex = 0), col = "black", fill = c("darkseagreen1","darkseagreen","thistle1","plum"), main = "Upregulated genes", shape="ellipse",plot=TRUE)
upreg_list <- ps_euler(ps_AM_upreg, "gendiet", fraction = 0.5, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","yellow"), main = "G3 Larvae", shape="ellipse",plot=FALSE)

# plot
upreg_plot
# list of genes
upreg_list

# plot the Euler diagrams side-by-side, roughly scaled to match circle sizes
plot_grid(upreg_plot,downreg_plot,nrow=1,rel_widths = c(1.7, 1))

# export the downreg and upreg lists for annotation with gProfiler
write.csv(upreg_res_DOXG3_H2OG3_AM_ordered[rownames(upreg_res_DOXG3_H2OG3_AM_ordered) %in% upreg_list$DOXG3_H2OG3,],"forgPro_upreg_res_DOXG3_H2OG3_AM_ordered.csv")
write.csv(upreg_res_DOXG4_H2OG4_AM_ordered[rownames(upreg_res_DOXG4_H2OG4_AM_ordered) %in% upreg_list$DOXG4_H2OG4,],"forgPro_upreg_res_DOXG4_H2OG4_AM_ordered.csv")
write.csv(downreg_res_DOXG3_H2OG3_AM_ordered[rownames(downreg_res_DOXG3_H2OG3_AM_ordered) %in% downreg_list$DOXG3_H2OG3,],"forgPro_downreg_res_DOXG3_H2OG3_AM_ordered.csv")
write.csv(downreg_res_DOXG4_H2OG4_AM_ordered[rownames(downreg_res_DOXG4_H2OG4_AM_ordered) %in% downreg_list$DOXG4_H2OG4,],"forgPro_downreg_res_DOXG4_H2OG4_AM_ordered.csv")

write.csv(ps_DOX_G3_upreg_list$A__L,"forgPro_shared_upreg.csv")
write.csv(downreg_list$A__L,"forgPro_shared_downreg.csv")    

# for personal use -- do not run
# sed -nE 's/.*gene_id "([a-zA-Z0-9_-]+)".*gene_biotype "([a-zA-Z_]+)".*/\1 \2/p' genomic.gtf >test.txt 
