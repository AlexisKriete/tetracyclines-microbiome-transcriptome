#----------------------------------------------------------------------
# Install packages and load libraries
#----------------------------------------------------------------------

install.packages("BiocManager")
library(BiocManager)
BiocManager::install("dada2")
library(dada2)
BiocManager::install("Biostrings")
library(Biostrings)
BiocManager::install("DECIPHER")
library(DECIPHER)
BiocManager::install("phyloseq")
library(phyloseq)
BiocManager::install("microbiome")
library(microbiome)
remotes::install_github("vmikk/metagMisc",force=TRUE)
library(metagMisc)
install.packages("SRS")
library(SRS)
install.packages("dplyr")
library(dplyr)
install.packages("tidyverse")
library(tidyverse)
install.packages("ape")
library(ape)
install.packages("picante")
library(picante)
BiocManager::install("MicEco")
library(MicEco)
remotes::install_github("jbisanz/qiime2R",force=TRUE)
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
library(microViz)

BiocManager::install("phyloseq")
BiocManager::install("microbiome")
install.packages("vegan")
install.packages("ggplot2")
BiocManager::install("MicEco")
BiocManager::install("ANCOMBC")
BiocManager::install("dada2")
install.packages("magrittr")
install.packages("Cairo")
install.packages("remotes")
install.packages("forcats")
install.packages("SRS")
install.packages("gridExtra")
install.packages("metacoder")
install.packages("ape")
install.packages("ggimage")
install.packages("microViz")
install.packages("multcompView")
install.packages("picante")
install.packages("cowplot")
install.packages("ggh4x")
remotes::install_github("jbisanz/qiime2R",force=TRUE)
install.packages("ggforce")

library(DECIPHER)
library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(ggplot2)
library(MicEco)
library(ANCOMBC)
library(dada2)
library(magrittr)
library(Cairo)
library(remotes)
remotes::install_github("vmikk/metagMisc",force=TRUE)
remotes::install_github("schuyler-smith/phylosmith",force=TRUE)
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(phylosmith)
library(metagMisc)
library(pairwiseAdonis)
library(forcats)
library(SRS)
library(gridExtra)
library(metacoder)
library(ape)
library(ggimage)
library(microViz)
library(multcompView)
library(picante)
library(cowplot)
library(ggh4x)
library(qiime2R)
library(ggforce)

#----------------------------------------------------------------------
# Import reads and filter with DADA2
#----------------------------------------------------------------------


## Define the path where the raw reads are located

path <- "/Users/alexiskriete/Desktop/Microbiome/Lc_antibiotic_study/New/raw-reads-demuxed"

list.files(path)

## set up reads for trimming/QC 

fnFs <- sort(list.files(path, pattern="1.fq", full.names = TRUE))

fnRs <- sort(list.files(path, pattern="2.fq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

sample.names

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names

names(filtRs) <- sample.names

## trim and filter reads: remove first 28 bases, truncate at base 195, remove unknown bases etc.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(28,28), truncLen=c(195,195),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

head(out)

plotQualityProfile(filtFs[1])

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
# dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

plot(table(nchar(getSequences(seqtab)))) 

# Note reads of non-target length at 191bp and ~300bp
# BLAST search of the 191bp seqs showed they are from mitochondria
# The first ~250bp of the 300bp seqs map to bacteria; the remaining ~50bp have no hits or map to weird taxa (plants etc.)
# Therefore, I used the following command to restrict read length 

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 245:260]

# We now have fewer ASVs
length(getSequences(seqtab))
length(getSequences(seqtab2))

# Now we remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

# After chimera removal, 79% of ASVs and 97% of reads remain.
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)

# View how many reads are retained at each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

# get the total # of reads in the filtered dataset, along with minimum, maximum and mean # of reads/sample
sum(track[,6]) # 8,164,337
min(track[,6]) # 70,501
max(track[,6]) # 256,874
mean(track[,6]) # 151,191.4

#----------------------------------------------------------------------
# Assign taxonomy using IDTAXA
#----------------------------------------------------------------------

load("/Users/alexiskriete/Desktop/Microbiome/Lc_antibiotic_study/New/raw-reads-demuxed/SILVA_SSU_r138_2019.RData")

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

ids <- list(1:length(getSequences(seqtab.nochim)))

for (x in 1:length(getSequences(seqtab.nochim)))
{if (x %% 100 == 0) {print(c("On sequence",x))}
  dna <- DNAStringSet(getSequences(seqtab.nochim)[x:x])
  ids[x] <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE)[1]}

taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)

taxa <- taxid

# read in sample metadata
sample_data <- utils::read.csv("/Users/alexiskriete/Desktop/Microbiome/Lc_antibiotic_study/New/raw-reads-demuxed/sample-metadata.csv", header=TRUE, row.names=1, stringsAsFactors = TRUE)

# make phyloseq object from the ASV table, sample metadata table, and taxonomy table
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sample_data), 
               tax_table(taxa))

# check the ps object
ps

# Add refseq slot to our phyloseq object that shows the DNA sequence associated with each ASV
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#### Filtering the ASV dataset ####

# remove singletons (i.e., ASVs with an abundance of 1 read, which must necessarily be only found in 1 sample)
# remove any ASVs classified as mitochondria or chloroplast (e.g. cyanobacteria)

# Check how many domains are represented in our dataset
get_taxa_unique(ps, "domain")

# See how many Bacteria and Archaea are in our dataset
ps_bacteria <- phyloseq::subset_taxa(ps, domain=="Bacteria")
phyloseq::ntaxa(ps_bacteria) # 23995 ASVs assigned to Bacteria (78% of total ASVs)
ps_archaea <- phyloseq::subset_taxa(ps, domain=="Archaea")
phyloseq::ntaxa(ps_archaea) # 48 ASVs assigned to Archaea (0.2% of total ASVs)

# 23% of ASVs failed to be assigned to anything (domain = NA)
1 - (phyloseq::ntaxa(ps_archaea) + phyloseq::ntaxa(ps_bacteria))/phyloseq::ntaxa(ps)

# Let's see how many singletons are in this dataset, along with some other stats
summarize_phyloseq(ps) # we have 65 singletons (0.21% of ASVs)

# Remove singletons
ps_nosingle <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 1, ps)

# Verify singletons were removed
summarize_phyloseq(ps_nosingle) 

# Remove ASVs assigned to mitochondria (Family)
ps_nomito <- phyloseq::subset_taxa(ps_nosingle, family!="Mitochondria"|is.na(family))
summarize_phyloseq(ps_nomito)

# Remove ASVs assigned to chloroplast (Order)
ps_nochloro <- phyloseq::subset_taxa(ps_nomito, order!="Chloroplast"|is.na(order))

# Check how many taxa have been retained so far after removing singletons, mitochondria and chloroplasts (99.5%)
phyloseq::ntaxa(ps_nochloro)/phyloseq::ntaxa(ps)

# Check out the Archaea
ps_filt_archaea <- phyloseq::subset_taxa(ps_nochloro, domain=="Archaea")
write.csv(as.data.frame(otu_table(ps_filt_archaea)), file="taxa table archaea.csv")

#### Exporting ASV sequences for cross-checking taxonomic assignments w/BLAST ####

# export .csv file that shows the sequences of the top 50 most abundant ASVs
write.csv(as.data.frame(refseq(ps_nochloro)[1:50,]), file="newASVseqs_top50.csv")

# export .csv file with the taxonomic assignments of the top 50 most abundant ASVs
write.csv(tax_table(ps_nochloro)[1:50],file="newASVnames_top50.csv")

# after cross-referencing by BLASTing to the 16S rRNA database, we can narrow down some taxonomic assignments
# can only find genus-level info at best; no species-level assignments are possible due to the short region sequenced
# the next few lines of code replace the old data for the top 50 most abundant ASVs with new dta

newNames <- as.matrix(read.csv("ASV_names_top50.csv"))
newNames2 <- newNames[,-1]
rownames(newNames2) <- newNames[,1]
tax_table(ps_nochloro)[1:50,] <- tax_table(newNames2)

# Next, we remove NAs (i.e., ASVs that could not be classified to any taxonomic level) and Archaea
# we can do this by keeping only ASVs assigned to Bacteria
ps_noNAs<- phyloseq::subset_taxa(ps_nochloro, domain=="Bacteria")

# we see that this removes ~23% of all ASVs, leaving us with 23877 ASVs total
ps_noNAs
ps_nochloro

#### Abundance/prevalence filtering ####

# first, we will remove ASVs found in only 1 sample (prevalence filtering)
ps_noNAs_prevfilt <- filter_taxa(ps_noNAs, function(x){sum(x > 0) > 1}, prune = TRUE)

# next, we will remove ASVs whose total abundance (summed across all 54 samples) is less than 79 reads
# 79 reads = 0.001% of taxa
# we can find this number using:
phyloseq_summary(ps_noNAs_prevfilt)[5,2] # 7,882,708 total reads
phyloseq_summary(ps_noNAs_prevfilt)[5,2]*.00001

# abundance filtering step
ps_final <- prune_taxa(taxa_sums(ps_noNAs_prevfilt) >= 79, ps_noNAs_prevfilt)

# we can see how many reads and ASVs are left over at each step
# we have lost 92% of ASVs by applying these 2 filtering steps!
# but, we retained 97% of reads

phyloseq_summary(ps_noNAs) # 8,020,203 reads
ps_noNAs # 23,877 taxa
phyloseq_summary(ps_noNAs_prevfilt) # 7,882708 reads
ps_noNAs_prevfilt # 5,122 taxa (21%)
phyloseq_summary(ps_final) # 7,794,147 reads
ps_final # 1,933 taxa (8%)

# Exporting phyloseq data to a csv file for safekeeping -- ASV table
OTU1 = as(otu_table(ps_final), "matrix")
write.csv(as.data.frame(OTU1), file="ps_final_ASV_table_raw.csv")

# Exporting phyloseq object data to a csv file for safekeeping -- taxa table
taxa1 = as(tax_table(ps_final), "matrix")
write.csv(as.data.frame(taxa1), file="ps_final_taxa_table_raw.csv")

# Exporting phyloseq object data to a csv file for safekeeping -- refseq
write.csv(as.data.frame(refseq(ps_final)), file="ps_final_refseq.csv")

# let's find the minimum, maximum, and mean number of reads per sample in the final dataset
ps_final_reads <- 1:54
for (x in 1:54) {ps_final_reads[x] <- sum(otu_table(ps_final)[x,])}
min(ps_final_reads) # 65,135
max(ps_final_reads) # 254,106
mean(ps_final_reads) # 144,336

ps_adults_norm <- phyloseq::subset_samples(ps_norm, stage=="Adult" & diet=="H2O")
ps_larvae_norm <- phyloseq::subset_samples(ps_norm, stage=="Larva" & diet=="H2O")
names(sort(colSums(otu_table(ps_adults_norm)),TRUE)[1:14])
names(sort(colSums(otu_table(ps_larvae_norm)),TRUE)[1:14])
# note that the 3 major dietary groups have similar total # of reads (sampling depth)
# same story for generation 3 samples vs. generation 4 samples
# we will still normalize our data, because some experimental groups have uneven sampling depths

ps_H2O_final <- phyloseq::subset_samples(ps_final, diet=="H2O")
sum(otu_table(ps_H2O_final))
ps_ATC_final <- phyloseq::subset_samples(ps_final, diet=="ATC")
sum(otu_table(ps_ATC_final))
ps_DOX_final <- phyloseq::subset_samples(ps_final, diet=="DOX")
sum(otu_table(ps_DOX_final))

ps_larvae_final <- phyloseq::subset_samples(ps_final, stage=="Larva")
ps_adultM_final <- phyloseq::subset_samples(ps_final, stage=="Adult" & sex=="M")
ps_adultF_final <- phyloseq::subset_samples(ps_final, stage=="Adult" & sex=="F")

#### Core microbiome and shared/unique taxa ####

# we will now turn to analyzing presence/absence data, such as shared/unique ASVs displayed with Euler plots
# using an un-normalized dataset is risky, because some samples might have more ASVs than others simply as a consequence of having a higher sampling depth
# we will use SRS (scaling with ranked subsampling) to normalize our samples to even sample size (Cmin)
# note that ASVs with count values of 1 *can* be converted to 0 values during downsampling; in other words, a presence can become an absence

SRS.shiny.app(as.data.frame(t(otu_table(ps_final)))) #Cmin = 65,135 reads

# make the normalized OTU table and put column names back
normalized_otu_table <- t(SRS(as.data.frame(t(otu_table(ps_final))),65135))
colnames(normalized_otu_table) <- colnames(otu_table(ps_final))

# make new phyloseq object with normalized ASV counts
ps_norm <- phyloseq(otu_table(normalized_otu_table, taxa_are_rows=FALSE), sample_data(sample_data), tax_table(tax_table(ps_final)))

# let's see how many ASVs were identified to genus, family, etc. level
1-(sum(is.na(tax_table(ps_norm)[,6]))/1936) # 46% IDed to genus
1-(sum(is.na(tax_table(ps_norm)[,5]))/1936) # 71% IDed to family
1-(sum(is.na(tax_table(ps_norm)[,4]))/1936) # 80% IDed to order
1-(sum(is.na(tax_table(ps_norm)[,3]))/1936) # 90% IDed to class
1-(sum(is.na(tax_table(ps_norm)[,2]))/1936) # 93% IDed to phylum

# this code will tell us how many genera, families etc. are in our final dataset
# 1 domain, 27 phyla, 62 classes, 136 orders, 189 families, 314 genera
tax_table(ps_norm) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  gather("Rank", "Name", rank_names(ps)) %>%
  na.omit() %>% # remove rows with NA value
  group_by(Rank) %>%
  summarize(ntaxa = length(unique(Name))) %>% # compute number of unique taxa
  mutate(Rank = factor(Rank, rank_names(ps))) %>%
  arrange(Rank)

ps_G3_A <- phyloseq::subset_samples(ps_norm, generation=="G3" & stage=="Adult")
ps_G4_A <- phyloseq::subset_samples(ps_norm, generation=="G4" & stage=="Adult")
ps_G3_L <- phyloseq::subset_samples(ps_norm, generation=="G3" & stage=="Larva")
ps_G4_L <- phyloseq::subset_samples(ps_norm, generation=="G4" & stage=="Larva")
ps_G3_AM <- phyloseq::subset_samples(ps_norm, generation=="G3" & stage=="Adult" & sex=="M")
ps_G4_AM <- phyloseq::subset_samples(ps_norm, generation=="G4" & stage=="Adult" & sex=="M")
#----------------------------------------------------------------------
# Euler plots to visualize shared and unique taxa
#----------------------------------------------------------------------

# code for generating Euler plots
# fraction = 0
G3L_euler <- ps_euler(ps_G3_L, "fullgroup", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1"), shape="ellipse",lty=1,plot=TRUE)
G3L_euler
G3A_euler <- ps_euler(ps_G3_A, "fullgroup", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","darkseagreen2", "thistle1", "lightskyblue1"),lty = c(3,5), shape="ellipse",plot=TRUE)
G3A_euler
G4L_euler <- ps_euler(ps_G4_L, "fullgroup", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1"), shape="ellipse",lty=1,plot=TRUE)
G4L_euler
G4A_euler <- ps_euler(ps_G4_A, "fullgroup", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","darkseagreen2", "thistle1", "lightskyblue1"),lty = c(3,5), shape="ellipse",plot=TRUE)
G4A_euler
G3Acomb_euler <- ps_euler(ps_G3_A, "group", fraction = 0.2, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1"),lty = 1, shape="ellipse",plot=TRUE)
G3Acomb_euler
G4Acomb_euler <- ps_euler(ps_G4_A, "group", fraction = 0.2, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1"),lty = 1, shape="ellipse",plot=TRUE)
G4Acomb_euler


### code for scaling Euler plots
# first, make vectors of the sums of ASVs in each Euler plot
# then, get scaling factor by taking the square root

list_G3L_euler <- ps_euler(ps_G3_L, "fullgroup", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1"), shape="ellipse",plot=FALSE)
list_G3A_euler <- ps_euler(ps_G3_A, "fullgroup", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","darkseagreen2", "thistle1", "lightskyblue1"),lty = 1:2, shape="ellipse",plot=FALSE)
list_G4L_euler <- ps_euler(ps_G4_L, "fullgroup", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1"), shape="ellipse",plot=FALSE)
list_G4A_euler <- ps_euler(ps_G4_A, "fullgroup", fraction = 0, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","darkseagreen2", "thistle1", "lightskyblue1"),lty = 1:2, shape="ellipse",plot=FALSE)
list_G3Acomb_euler <- ps_euler(ps_G3_A, "group", fraction = 0.2, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1"), shape="ellipse",plot=FALSE)
list_G4Acomb_euler <- ps_euler(ps_G4_A, "group", fraction = 0.2, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1"), shape="ellipse",plot=FALSE)

list_G3AM_euler <- ps_euler(ps_G3_AM, "fullgroup", fraction = 1, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","darkseagreen2", "thistle1", "lightskyblue1"),lty = 1:2, shape="ellipse",plot=FALSE)
list_G4AM_euler <- ps_euler(ps_G4_AM, "fullgroup", fraction = 1, labels = list(cex = 0), col = "black", fill = c("darkseagreen2", "thistle1", "lightskyblue1","darkseagreen2", "thistle1", "lightskyblue1"),lty = 1:2, shape="ellipse",plot=FALSE)


sums_for_scaling <- c(sum(lengths(list_G3L_euler)),sum(lengths(list_G3A_euler)),sum(lengths(list_G4L_euler)),sum(lengths(list_G4A_euler)),sum(lengths(list_G3Acomb_euler)),sum(lengths(list_G4Acomb_euler)))
scaling_factor <- sqrt(sums_for_scaling/sums_for_scaling[1])

# make the scaled plots
gridExtra::grid.arrange(G3L_euler, G3A_euler, G4L_euler, G4A_euler,G3Acomb_euler,G4Acomb_euler, widths=scaling_factor)

# get info for labeling the Euler plots
lengths(list_G3L_euler)
sum(lengths(list_G3L_euler))
lengths(list_G3L_euler)/sum(lengths(list_G3L_euler))

lengths(list_G4L_euler)
sum(lengths(list_G4L_euler))
lengths(list_G4L_euler)/sum(lengths(list_G4L_euler))

lengths(list_G3Acomb_euler)
sum(lengths(list_G3Acomb_euler))
lengths(list_G3Acomb_euler)/sum(lengths(list_G3Acomb_euler))

lengths(list_G4Acomb_euler)
sum(lengths(list_G4Acomb_euler))
lengths(list_G4Acomb_euler)/sum(lengths(list_G4Acomb_euler))

#----------------------------------------------------------------------
# Alpha diversity
#----------------------------------------------------------------------

### get Faith PD, Pielou Evenness and Shannon Index values into dataframe with sample info for plotting

# we will calculate Faith PD, so we need a phylogenetic tree, which we will make with the 'ape' package
# create the tree with "rtree", and use set.seed() to ensure reproducibility

set.seed(2)
random_tree = rtree(ntaxa(ps_norm), rooted=TRUE, tip.label=taxa_names(ps_norm))

# make "alpha" data frame with sample info and Shannon index data
alpha <- plot_richness(ps_norm,x="group",measures="Shannon",color="group")$data

# rename the "value" columns to ShIndex
colnames(alpha)[10] <- "ShIndex"

# add a column for the Faith PD values
alpha[,12] <- as.vector(pd(samp = otu_table(ps_norm), tree = random_tree)[1])

# add a column for the Pielou evenness values
alpha[,13] <- evenness(ps_norm, index = "all", zeroes = TRUE, detection = 0)[,2]

# rename the column
colnames(alpha)[13] <- "Evenness"

#### Hypothesis testing (Shannon Index) ####

# global test for differences among any experimental group
kruskal.test(alpha[,10] ~ alpha$group, data = alpha)

# pairwise tests with Benjamini-Hochberg correction for multiple comparisons
Shannon_pairwise <- pairwise.wilcox.test(alpha[,10], alpha$group ,p.adjust.method = "BH")

# p-value matrix
Shannon_pairwise$p.value

# we can export the p-value matrix to a .csv file
write.csv(Shannon_pairwise$p.value,"Shannon_pvals_pairs.csv")
# p-values < 0.05
Shannon_pairwise$p.value <0.05

#### Hypothesis testing (Faith PD) ####

# global test
kruskal.test(alpha[,12] ~ alpha$group, data = alpha)
# pairwise tests with Benjamini-Hochberg correction for multiple comparisons
PD_pairwise <- pairwise.wilcox.test(alpha[,12], alpha$group ,p.adjust.method = "BH")
# p-value matrix
PD_pairwise$p.value
# we can export the p-value matrix to a .csv file
write.csv(PD_pairwise$p.value,"PD_pvals_pairs.csv")
# p-values < 0.05
PD_pairwise$p.value <0.05

##### Hypothesis testing (Evenness) ####

# global test
kruskal.test(alpha[,13] ~ alpha$group, data = alpha)
# pairwise tests with Benjamini-Hochberg correction for multiple comparisons
Evenness_pairwise <- pairwise.wilcox.test(alpha[,13], alpha$group ,p.adjust.method = "BH")
# p-value matrix
Evenness_pairwise$p.value
# we can export the p-value matrix to a .csv file
write.csv(Evenness_pairwise$p.value,"Evenness_pvals_pairs.csv")
# p-values < 0.05
Evenness_pairwise$p.value <0.05

#### Alpha diversity plots ####

# we will first figure out which groups are significantly different from each other (p < 0.05)
# since we have so many comparisons, we will designate groups by letters
# let's first get data into appropriate format
Shannon_pairwise_df <- as.data.frame(Shannon_pairwise$p.value)
PD_pairwise_df <- as.data.frame(PD_pairwise$p.value)
Evenness_pairwise_df <- as.data.frame(Evenness_pairwise$p.value)

Shannon_pairwise_df <- add_row(Shannon_pairwise_df,.before = 1)
PD_pairwise_df <- add_row(PD_pairwise_df,.before = 1)
Evenness_pairwise_df <- add_row(Evenness_pairwise_df,.before = 1)

rownames(Shannon_pairwise_df)[1] <- "G3 ATC adults"
rownames(PD_pairwise_df)[1] <- "G3 ATC adults"
rownames(Evenness_pairwise_df)[1] <- "G3 ATC adults"

Shannon_pairwise_df <- cbind(Shannon_pairwise_df,empty_column=NA)
PD_pairwise_df <- cbind(PD_pairwise_df,empty_column=NA)
Evenness_pairwise_df <- cbind(Evenness_pairwise_df,empty_column=NA)

colnames(Shannon_pairwise_df)[12] <- "G4 H2O larvae"
colnames(PD_pairwise_df)[12] <- "G4 H2O larvae"
colnames(Evenness_pairwise_df)[12] <- "G4 H2O larvae"

# make the matrices of p-values symmetrical
Shannon_pairwise_df.sym <- Matrix::forceSymmetric(as.matrix(Shannon_pairwise_df),uplo="L")
PD_pairwise_df.sym <- Matrix::forceSymmetric(as.matrix(PD_pairwise_df),uplo="L")
Evenness_pairwise_df.sym <- Matrix::forceSymmetric(as.matrix(Evenness_pairwise_df),uplo="L")

# add new columns to the alpha data
alpha[,14:16] <- c(NA,NA,NA)

# get letters denoting significantly different groups (p<0.05)
letters <- multcompLetters(Shannon_pairwise_df.sym)
for (x in 1:54) { alpha[x, 14] <- letters$Letters[[alpha$group[x]]]}
letters <- multcompLetters(PD_pairwise_df.sym)
for (x in 1:54) { alpha[x, 15] <- letters$Letters[[alpha$group[x]]]}
letters <- multcompLetters(Evenness_pairwise_df.sym)
for (x in 1:54) { alpha[x, 16] <- letters$Letters[[alpha$group[x]]]}

# make new columns in alpha with max values -- this will help us plot the letters above data points on the graph
for (x in 1:54) { alpha[x, 17] <- max(alpha[alpha$group == alpha$group[x],10])}
for (x in 1:54) { alpha[x, 18] <- max(alpha[alpha$group == alpha$group[x],12])}
for (x in 1:54) { alpha[x, 19] <- max(alpha[alpha$group == alpha$group[x],13])}

# rename columns for easier plotting
colnames(alpha)[14:19] <- c("letters_ShIndex","letters_PD","letters_Evenness","maxes_ShIndex","maxes_PD","maxes_Evenness") 

# rename the stages from singular to plural
levels(alpha$stage) <- c("Adults","Larvae")

### Code for alpha diversity plots ###
# first, set colors for the stage panels
strip <- strip_themed(background_x = elem_list_rect(fill = c("#FFE7BA","#DBDBDB")))

# Shannon Index plot (article version)
Sh_plot_final <- ggplot(alpha, aes(forcats::fct_relevel(diet, "DOX", "ATC", "H2O"), ShIndex, group=generation, color=group)) + 
scale_color_manual(values = c("#35C47D","#35C47D","#CD96CD","#CD96CD","#5CACEE","#5CACEE","#199151","#199151","#8B668B","#8B668B","#1874CD","#1874CD")) + 
geom_point(size=2.5,shape=16,alpha=1,position = position_jitterdodge(jitter.width = 0, dodge.width = .7)) + 
theme(strip.background = element_rect(colour = "black"),panel.border = element_rect(linetype = "solid", fill = NA), panel.background = element_rect(fill = "white", colour = "grey30"), panel.grid.minor=element_blank(), axis.title.x = element_text(size=15),  axis.title.y = element_text(size=15), axis.text.y = element_text(size=14), axis.text.x = element_text(size=12, angle = 0, vjust = 0, hjust=0.5), strip.text.x= element_text(size=20), legend.text = element_text(size=10)) + 
scale_y_continuous(limits = c(0, 7)) + 
stat_summary(fun.data=mean_se, color=c(rep(c("#CD96CD","#35C47D","#5CACEE","#8B668B","#199151","#1874CD"),2)),position=position_dodge(width=0.7), size=0.6, width=0.4, geom="errorbar") + 
stat_summary(fun.data=mean_se, color=c(rep(c("#CD96CD","#35C47D","#5CACEE","#8B668B","#199151","#1874CD"),2)),geom="point",position=position_dodge(width=0.7), size=2.5, shape=15) +
xlab("G1-G3 Diet") + 
ylab("Shannon Index") + 
geom_text(aes(label = str_trim(letters_ShIndex), y = rep(6.4,54)), vjust = -0.5, position=position_dodge(width=0.7),size=4,color="black") + 
theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
facet_wrap2(~forcats::fct_relevel(stage, "Larvae", "Adults"), strip = strip)

Sh_plot_final

# Faith PD plot
PD_plot_final <- ggplot(alpha, aes(forcats::fct_relevel(diet, "DOX", "ATC", "H2O"), PD, group=generation, color=group)) + 
  scale_color_manual(values = c("#35C47D","#35C47D","#CD96CD","#CD96CD","#5CACEE","#5CACEE","#199151","#199151","#8B668B","#8B668B","#1874CD","#1874CD")) + 
  geom_point(size=2.5,shape=16,alpha=1,position = position_jitterdodge(jitter.width = 0, dodge.width = .7)) + 
  theme(strip.background = element_blank(),panel.border = element_rect(linetype = "solid", fill = NA), panel.background = element_rect(fill = "white", colour = "grey30"), panel.grid.minor=element_blank(), axis.title.x = element_text(size=15),  axis.title.y = element_text(size=15), axis.text.y = element_text(size=14), axis.text.x = element_text(size=12, angle = 0, vjust = 0, hjust=0.5), strip.text.x= element_blank(), legend.text = element_text(size=10)) + 
  scale_y_continuous(limits = c(0, 2000)) + 
  stat_summary(fun.data=mean_se, color=c(rep(c("#CD96CD","#35C47D","#5CACEE","#8B668B","#199151","#1874CD"),2)),position=position_dodge(width=0.7), size=0.6, width=0.4, geom="errorbar") + 
  stat_summary(fun.data=mean_se, color=c(rep(c("#CD96CD","#35C47D","#5CACEE","#8B668B","#199151","#1874CD"),2)),geom="point",position=position_dodge(width=0.7), size=2.5, shape=15) +
  xlab("G1-G3 Diet") + 
  ylab("Faith PD") + 
  geom_text(aes(label = str_trim(letters_PD), y = rep(1900,54)), vjust = -0.5, position=position_dodge(width=0.7),size=4,color="black") + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  facet_wrap2(~forcats::fct_relevel(stage, "Larvae", "Adults"))
PD_plot_final

# Evenness plot
Evenness_plot_final <- ggplot(alpha, aes(forcats::fct_relevel(diet, "DOX", "ATC", "H2O"), Evenness, group=generation, color=group)) + 
  scale_color_manual(values = c("#35C47D","#35C47D","#CD96CD","#CD96CD","#5CACEE","#5CACEE","#199151","#199151","#8B668B","#8B668B","#1874CD","#1874CD")) + 
  geom_point(size=2.5,shape=16,alpha=1,position = position_jitterdodge(jitter.width = 0, dodge.width = .7)) + 
  theme(strip.background = element_blank(),panel.border = element_rect(linetype = "solid", fill = NA), panel.background = element_rect(fill = "white", colour = "grey30"), panel.grid.minor=element_blank(), axis.title.x = element_text(size=15),  axis.title.y = element_text(size=15), axis.text.y = element_text(size=14), axis.text.x = element_text(size=12, angle = 0, vjust = 0, hjust=0.5), strip.text.x= element_blank(), legend.text = element_text(size=10)) + 
  scale_y_continuous(limits = c(0, 1)) + 
  stat_summary(fun.data=mean_se, color=c(rep(c("#CD96CD","#35C47D","#5CACEE","#8B668B","#199151","#1874CD"),2)),position=position_dodge(width=0.7), size=0.6, width=0.4, geom="errorbar") + 
  stat_summary(fun.data=mean_se, color=c(rep(c("#CD96CD","#35C47D","#5CACEE","#8B668B","#199151","#1874CD"),2)),geom="point",position=position_dodge(width=0.7), size=2.5, shape=15) +
  xlab("G1-G3 Diet") + ylab("Evenness") + 
  geom_text(aes(label = str_trim(letters_Evenness), y = rep(0.93,54)), vjust = -0.5, position=position_dodge(width=0.7),size=4,color="black") + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap2(~forcats::fct_relevel(stage, "Larvae", "Adults"), strip = strip)
Evenness_plot_final
# Arrange alpha diversity plots in a vertical stack
# export as PDF: portrait mode, 15"x8"
alpha_div_vertical <- plot_grid(Sh_plot_final, PD_plot_final, Evenness_plot_final, ncol=1, align="v")
alpha_div_vertical

######## Beta diversity and ordination ########

### Use the phyloseq2qiime2 function to get our phyloseq object, ps_final, into a qiime2-readable format
# author credit: Christian Edwardson (https://github.com/cedwardson4)
# the function's code can be viewed at: https://github.com/cedwardson4/phyloseq-extras/blob/master/phyloseq2QIIME2.R

phyloseq2qiime2<-function(physeq){
  #take a phyloseq object,check for individual parts, write to files ready for qiime2 upload
  library(biomformat)
  if(packageVersion("biomformat") < "1.7") {
    stop("This will only work with biomformat version > 1.7")
  }
  ps_name <-deparse(substitute(physeq))
  taxa_are_rows_logical<-taxa_are_rows(physeq)
  #write OTU table to biom file
  if(is.null(access(physeq,"otu_table"))==FALSE){
    if(taxa_are_rows_logical==TRUE) {
      otu<-as(otu_table(physeq),"matrix")
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom,biom_file=paste0(ps_name,"_features-table.biom"))
      print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
    } else if (taxa_are_rows_logical==FALSE) {
      otu<-t(as(otu_table(physeq),"matrix"))
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom,biom_file=paste0(ps_name,"_feature-table.biom"))
      print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
    }
  }
  #write sample data (metadata) to tsv
  if(is.null(access(physeq,"sam_data"))==FALSE){
    write.table(data.frame(sample_data(physeq)),file=paste0(ps_name,"_sample-metadata.txt"), 
                sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    print(paste0("Writing sample metadata to ",ps_name,"_sample-metadata.txt"))
  }
  #write taxonomy table to qiime2 formatted taxonomy
  if(is.null(access(physeq,"tax_table"))==FALSE){
    tax<-as(tax_table(physeq),"matrix")
    tax_cols <- colnames(tax)
    tax<-as.data.frame(tax)
    tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
    for(co in tax_cols) tax[co]<-NULL
    write.table(tax, file=paste0(ps_name,"_tax.txt"), quote=FALSE, col.names=FALSE, sep="\t")
    print(paste0("Writing taxonomy table to ",ps_name,"_tax.txt"))
  }
  #write phylogenetic tree to newick formwat
  if(is.null(access(physeq,"phy_tree"))==FALSE){
    if(is.rooted(phy_tree(physeq))==TRUE) {
      ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-rooted.newick"))
      print(paste0("Writing rooted tree to ",ps_name,"_tree-rooted.newick"))
    } else if (is.rooted(phy_tree(physeq))==FALSE) {
      ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-unrooted.newick"))
      print(paste0("Writing unrooted tree to ",ps_name,"_tree-unrooted.newick"))
    }      
  }
  #write representative sequences to fasta format
  if(is.null(access(physeq,"refseq"))==FALSE){
    writeXStringSet(refseq(physeq),filepath=paste0(ps_name,"_ref-seqs.fasta"))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
  } else if (taxa_are_rows_logical==FALSE && unique(grepl("[^ATCG]",colnames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
    uniquesToFasta(t(otu), fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))      
  } else if (taxa_are_rows_logical==TRUE && unique(grepl("[^ATCG]",rownames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
    uniquesToFasta(otu, fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))    
  }
}

# add new "id" column to the sample metadata prior to exporting the metadata as .txt file

ps_final <- prune_taxa(taxa_sums(ps_noNAs_prevfilt) >= 79, ps_noNAs_prevfilt)
ps_larvae_final <- phyloseq::subset_samples(ps_final, stage=="Larva")
sample_data(ps_larvae_final) <- cbind(rownames(sample_data(ps_larvae_final)),sample_data(ps_larvae_final))
colnames(sample_data(ps_larvae_final))[1] <- "id"

ps_final <- prune_taxa(taxa_sums(ps_noNAs_prevfilt) >= 79, ps_noNAs_prevfilt)
ps_adults_final <- phyloseq::subset_samples(ps_final, stage=="Adult")
sample_data(ps_adults_final) <- cbind(rownames(sample_data(ps_adults_final)),sample_data(ps_adults_final))
colnames(sample_data(ps_adults_final))[1] <- "id"

ps_final <- prune_taxa(taxa_sums(ps_noNAs_prevfilt) >= 79, ps_noNAs_prevfilt)
sample_data(ps_final)$fullgroup <- paste(sample_data(ps_final)$sex,sample_data(ps_final)$group)
ps_adults_sexed_final <- phyloseq::subset_samples(ps_final, stage=="Adult")
sample_data(ps_adults_sexed_final) <- cbind(rownames(sample_data(ps_adults_sexed_final)),sample_data(ps_adults_sexed_final))
colnames(sample_data(ps_adults_sexed_final))[1] <- "id"

ps_final <- prune_taxa(taxa_sums(ps_noNAs_prevfilt) >= 79, ps_noNAs_prevfilt)
sample_data(ps_final) <- cbind(rownames(sample_data(ps_final)),sample_data(ps_final))
colnames(sample_data(ps_final))[1] <- "id"

# run the phyloseq2qiime2 functions to convert our phyloseq objects into qiime-readable data
# note: the output files are automatically created in our working directory
# we will make a new folder, "QIIME2", and put the files in this directory
setwd("Data/QIIME2")
phyloseq2qiime2(ps_final)
phyloseq2qiime2(ps_larvae_final)
phyloseq2qiime2(ps_adults_final)
phyloseq2qiime2(ps_adults_sexed_final)

### Command line code for running DEICODE in qiime2 (ordination + beta diversity pairwise comparisons with PERMANOVA) ###

# We need to install qiime2 and deicode first; installation method depends on your device
# See instructions for qiime2 here: https://docs.qiime2.org/2023.7/
# Instructions for installing DEICODE here: https://library.qiime2.org/plugins/deicode/19/
# Make sure to activate qiime2 environment before running code below!! This code will NOT run in R.

# The qiime commands below will run DEICODE (rclr-transforming our data, calculating distance matrix, and ordinating for biplot visualization)
# We can view the biplots by dragging the "biplot.qzv" files into qiime2view: https://view.qiime2.org/
# The top 20 differentiating features (ASVs) will be plotted on the biplot as arrows
# We also run pairwise comparisons of all groups using PERMANOVA, which tells us which groups differ from each other in terms of beta diversity
# To view the PERMANOVA data, including q values (aka p values after correcting for multiple comparisons using BH/FDR method), see .csv files in supplement 

# for all 54 samples (larvae and adults)
qiime tools import --input-path ps_final_feature-table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path feature-table-all.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path ps_final_tax.txt --output-path taxonomy-all.qza
qiime deicode rpca --i-table feature-table-all.qza --o-biplot ordination-all.qza --o-distance-matrix distance-all.qza
qiime emperor biplot --i-biplot ordination-all.qza --m-sample-metadata-file ps_final_sample-metadata.txt --m-feature-metadata-file taxonomy-all.qza --o-visualization biplot-all.qzv --p-number-of-features 20
qiime diversity beta-group-significance --i-distance-matrix distance-all.qza --m-metadata-file ps_final_sample-metadata.txt --m-metadata-column group --p-method permanova --p-pairwise --o-visualization group_significance-all.qzv

# for larvae only (N=18 samples)
qiime tools import --input-path ps_larvae_final_feature-table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path feature-table-larvae.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path ps_larvae_final_tax.txt --output-path taxonomy-larvae.qza
qiime deicode rpca --i-table feature-table-larvae.qza --o-biplot ordination_larvae.qza --o-distance-matrix distance_larvae.qza
qiime emperor biplot --i-biplot ordination_larvae.qza --m-sample-metadata-file ps_larvae_final_sample-metadata.txt --m-feature-metadata-file taxonomy-larvae.qza --o-visualization biplot-larvae.qzv --p-number-of-features 20
qiime diversity beta-group-significance --i-distance-matrix distance_larvae.qza --m-metadata-file ps_larvae_final_sample-metadata.txt --m-metadata-column group --p-method permanova --p-pairwise --o-visualization group_significance-larvae.qzv

# for adults only (N=36 samples)
qiime tools import --input-path ps_adults_sexed_final_feature-table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path feature-table-adults-sexed.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path ps_adults_sexed_final_tax.txt --output-path taxonomy-adults-sexed.qza
qiime deicode rpca --i-table feature-table-adults-sexed.qza --o-biplot ordination_adults-sexed.qza --o-distance-matrix distance_adults-sexed.qza
qiime emperor biplot --i-biplot ordination_adults-sexed.qza --m-sample-metadata-file ps_adults_sexed_final_sample-metadata.txt --m-feature-metadata-file taxonomy-adults-sexed.qza --o-visualization biplot-adults-sexed.qzv --p-number-of-features 20

# PERMANOVA by sex (test statistic = 0.02913, p value = 0.987)
qiime diversity beta-group-significance --i-distance-matrix distance_adults-sexed.qza --m-metadata-file ps_adults_sexed_final_sample-metadata.txt --m-metadata-column fullgroup --p-method permanova --p-pairwise --o-visualization sex_significance-adults.qzv
qiime diversity beta-group-significance --i-distance-matrix distance_larvae.qza --m-metadata-file ps_larvae_final_sample-metadata.txt --m-metadata-column group --p-method permanova --p-pairwise --o-visualization significance-larvae.qzv


### From command line, run qurro to identify differentiating features

# all samples
#qiime qurro loading-plot --i-table feature-table-all.qza --i-ranks ordination-all.qza --m-sample-metadata-file ps_final_sample-metadata.txt --m-feature-metadata-file taxonomy-all.qza --o-visualization qurro/qurro-plot-all.qzv 
# larvae
#qiime qurro loading-plot --i-table feature-table-larvae.qza --i-ranks ordination_larvae.qza --m-sample-metadata-file ps_larvae_final_sample-metadata.txt --m-feature-metadata-file taxonomy-larvae.qza --o-visualization qurro/qurro-plot-larvae.qzv
# adults
#qiime qurro loading-plot --i-table feature-table-adults.qza --i-ranks ordination_adults.qza --m-sample-metadata-file ps_adults_final_sample-metadata.txt --m-feature-metadata-file taxonomy-adults.qza --o-visualization qurro/qurro-plot-adults.qzv 

###### Graphics for DEICODE PCoA plot ######

# first, import fly and larva icons for PCoA plot
myvec <- sample_data(ps_final)$generation=="G3" & sample_data(ps_final)$stage=="Larva"
vec1 <- replace(myvec,myvec==TRUE,'larva.png')
vec2 <- replace(vec1,sample_data(ps_final)$generation=="G4" & sample_data(ps_final)$stage=="Larva",'larva_border.png')
vec3 <- replace(vec2,sample_data(ps_final)$generation=="G3" & sample_data(ps_final)$stage=="Adult",'adult.png')
vec4 <- replace(vec3,sample_data(ps_final)$generation=="G4" & sample_data(ps_final)$stage=="Adult",'adult_border.png')

# PCoA with all 54 samples
ord<-read_qza("ordination-all.qza")
tax<-read_qza("taxonomy-all.qza")$data %>% dplyr::rename(FeatureID = Feature.ID)
meta<-read_tsv(file = "ps_final_sample-metadata.txt") %>% 
  dplyr::rename(SampleID=id) %>%
  filter(SampleID!="#q2:types")

# find proportion of variance explained by principal components axis 1 (66.36%)
round(100*ord$data$ProportionExplained[1],2)

# find proportion of variance explained by principal components axis 2 (32.15%)
round(100*ord$data$ProportionExplained[2],2) 

PCoA_all_plot <-
  ggplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position="none") +
  xlab("Axis 1 (66.36 %)") +
  ylab("Axis 2 (32.15 %)") + 
  geom_image(data=ord$data$Vectors %>% left_join(meta), aes(x=PC1, y=PC2, image=vec4, color=group),size=0.05,alpha=0.8) +
  scale_color_manual(values = c("G3 ATC larvae"="limegreen", "G3 DOX larvae"="#f352e4", "G3 H2O larvae"="steelblue1", "G3 ATC adults"="limegreen", "G3 DOX adults"="#f352e4", "G3 H2O adults"="steelblue1","G4 ATC larvae"="limegreen", "G4 DOX larvae"="#f352e4", "G4 H2O larvae"="steelblue1", "G4 ATC adults"="limegreen", "G4 DOX adults"="#f352e4", "G4 H2O adults"="steelblue1"))

PCoA_all_plot

######## Differential abundance testing with ANCOM-BC2 #######

# run ANCOMBC-2
ancombc2_allcomps$res_pair <- ancombc2(data=ps_final,fix_formula="group",group="group",pairwise=TRUE,p_adj_method = "bonferroni",struc_zero=FALSE, prv_cut=0)
# Export ANCOM-BC results to .csv files
write.csv(ancombc2_allcomps$res_pair,"ANCOMBC2_comps_pairwise_allgroups.csv")

# Run ANCOMBC2
ps_adult_males <- phyloseq::subset_samples(ps_final, stage=="Adult" & sex=="M")
ancom_adult_males <- ancombc2(data=ps_adult_males,fix_formula="group",group="group",pairwise=TRUE,p_adj_method = "holm",struc_zero=TRUE, prv_cut=0)

ps_adult_females <- phyloseq::subset_samples(ps_final, stage=="Adult" & sex=="F")
ancom_adult_females <- ancombc2(data=ps_adult_females,fix_formula="group",group="group",pairwise=TRUE,p_adj_method = "holm",struc_zero=TRUE, prv_cut=0)

ps_larvae <- phyloseq::subset_samples(ps_final, stage=="Larva")
ancom_larvae <- ancombc2(data=ps_larvae,fix_formula="group",group="group",pairwise=TRUE,p_adj_method = "holm",struc_zero=TRUE, prv_cut=0)

# DOX G3 L vs. H2O G3 L (5, 188, 347)
ancom_larvae$res_pair[ancom_larvae$res_pair$`diff_groupG3 H2O larvae_groupG3 DOX larvae`,c("taxon","lfc_groupG3 H2O larvae_groupG3 DOX larvae","se_groupG3 H2O larvae_groupG3 DOX larvae","q_groupG3 H2O larvae_groupG3 DOX larvae")]
dim(ancom_larvae$zero_ind[ancom_larvae$zero_ind$`structural_zero (group = G3 DOX larvae)` & !(ancom_larvae$zero_ind$`structural_zero (group = G3 H2O larvae)`),])
dim(ancom_larvae$zero_ind[ancom_larvae$zero_ind$`structural_zero (group = G3 H2O larvae)` & !(ancom_larvae$zero_ind$`structural_zero (group = G3 DOX larvae)`),])
# ATC G3 L vs. H2O G3 L (6, 174, 325)
ancom_larvae$res_pair[ancom_larvae$res_pair$`diff_groupG3 H2O larvae`,c("taxon","lfc_groupG3 H2O larvae","se_groupG3 H2O larvae","q_groupG3 H2O larvae")]
dim(ancom_larvae$zero_ind[ancom_larvae$zero_ind$`structural_zero (group = G3 ATC larvae)` & !(ancom_larvae$zero_ind$`structural_zero (group = G3 H2O larvae)`),])
dim(ancom_larvae$zero_ind[ancom_larvae$zero_ind$`structural_zero (group = G3 H2O larvae)` & !(ancom_larvae$zero_ind$`structural_zero (group = G3 ATC larvae)`),])
# DOX G4 L vs. H2O G4 L (1, 184, 344)
ancom_larvae$res_pair[ancom_larvae$res_pair$`diff_groupG4 H2O larvae_groupG4 DOX larvae`,c("taxon","lfc_groupG4 H2O larvae_groupG4 DOX larvae","se_groupG4 H2O larvae_groupG4 DOX larvae","q_groupG4 H2O larvae_groupG4 DOX larvae")]
dim(ancom_larvae$zero_ind[ancom_larvae$zero_ind$`structural_zero (group = G4 DOX larvae)` & !(ancom_larvae$zero_ind$`structural_zero (group = G4 H2O larvae)`),])
dim(ancom_larvae$zero_ind[ancom_larvae$zero_ind$`structural_zero (group = G4 H2O larvae)` & !(ancom_larvae$zero_ind$`structural_zero (group = G4 DOX larvae)`),])
# ATC G4 L vs. H2O G4 L (0, 327, 117)
ancom_larvae$res_pair[ancom_larvae$res_pair$`diff_groupG4 H2O larvae_groupG4 ATC larvae`,c("taxon","lfc_groupG4 H2O larvae_groupG4 ATC larvae","se_groupG4 H2O larvae_groupG4 ATC larvae","q_groupG4 H2O larvae_groupG4 ATC larvae")]
dim(ancom_larvae$zero_ind[ancom_larvae$zero_ind$`structural_zero (group = G4 ATC larvae)` & !(ancom_larvae$zero_ind$`structural_zero (group = G4 H2O larvae)`),])
dim(ancom_larvae$zero_ind[ancom_larvae$zero_ind$`structural_zero (group = G4 H2O larvae)` & !(ancom_larvae$zero_ind$`structural_zero (group = G4 ATC larvae)`),])

### adult male comparisons
# DOX G3 A vs. H2O G3 A (0, 150, 200)
ancom_adult_males$res_pair[ancom_adult_males$res_pair$`diff_groupG3 H2O adults_groupG3 DOX adults`,c("taxon","lfc_groupG3 H2O adults_groupG3 DOX adults","se_groupG3 H2O adults_groupG3 DOX adults","q_groupG3 H2O adults_groupG3 DOX adults")]
dim(ancom_adult_males$zero_ind[ancom_adult_males$zero_ind$`structural_zero (group = G3 DOX adults)` & !(ancom_adult_males$zero_ind$`structural_zero (group = G3 H2O adults)`),])
dim(ancom_adult_males$zero_ind[ancom_adult_males$zero_ind$`structural_zero (group = G3 H2O adults)` & !(ancom_adult_males$zero_ind$`structural_zero (group = G3 DOX adults)`),])
# ATC G3 A vs. H2O G3 A (0, 155, 236)
ancom_adult_males$res_pair[ancom_adult_males$res_pair$`diff_groupG3 H2O adults`,c("taxon","lfc_groupG3 H2O adults","se_groupG3 H2O adults","q_groupG3 H2O adults")]
dim(ancom_adult_males$zero_ind[ancom_adult_males$zero_ind$`structural_zero (group = G3 ATC adults)` & !(ancom_adult_males$zero_ind$`structural_zero (group = G3 H2O adults)`),])
dim(ancom_adult_males$zero_ind[ancom_adult_males$zero_ind$`structural_zero (group = G3 H2O adults)` & !(ancom_adult_males$zero_ind$`structural_zero (group = G3 ATC adults)`),])
# DOX G4 A vs. H2O G4 A (7, 1044, 236)
dim(ancom_adult_males$res_pair[ancom_adult_males$res_pair$`diff_groupG4 H2O adults_groupG4 DOX adults`,c("taxon","lfc_groupG4 H2O adults_groupG4 DOX adults","se_groupG4 H2O adults_groupG4 DOX adults","q_groupG4 H2O adults_groupG4 DOX adults")])
dim(ancom_adult_males$zero_ind[ancom_adult_males$zero_ind$`structural_zero (group = G4 DOX adults)` & !(ancom_adult_males$zero_ind$`structural_zero (group = G4 H2O adults)`),])
dim(ancom_adult_males$zero_ind[ancom_adult_males$zero_ind$`structural_zero (group = G4 H2O adults)` & !(ancom_adult_males$zero_ind$`structural_zero (group = G4 DOX adults)`),])
# ATC G4 A vs. H2O G4 A (0, 866, 69)
ancom_adult_males$res_pair[ancom_adult_males$res_pair$`diff_groupG4 H2O adults_groupG4 ATC adults`,c("taxon","lfc_groupG4 H2O adults_groupG4 ATC adults","se_groupG4 H2O adults_groupG4 ATC adults","q_groupG4 H2O adults_groupG4 ATC adults")]
dim(ancom_adult_males$zero_ind[ancom_adult_males$zero_ind$`structural_zero (group = G4 ATC adults)` & !(ancom_adult_males$zero_ind$`structural_zero (group = G4 H2O adults)`),])
dim(ancom_adult_males$zero_ind[ancom_adult_males$zero_ind$`structural_zero (group = G4 H2O adults)` & !(ancom_adult_males$zero_ind$`structural_zero (group = G4 ATC adults)`),])

### adult female comparisons
# DOX G3 A vs. H2O G3 A (0, 263, 166)
ancom_adult_females$res_pair[ancom_adult_females$res_pair$`diff_groupG3 H2O adults_groupG3 DOX adults`,c("taxon","lfc_groupG3 H2O adults_groupG3 DOX adults","se_groupG3 H2O adults_groupG3 DOX adults","q_groupG3 H2O adults_groupG3 DOX adults")]
dim(ancom_adult_females$zero_ind[ancom_adult_females$zero_ind$`structural_zero (group = G3 DOX adults)` & !(ancom_adult_females$zero_ind$`structural_zero (group = G3 H2O adults)`),])
dim(ancom_adult_females$zero_ind[ancom_adult_females$zero_ind$`structural_zero (group = G3 H2O adults)` & !(ancom_adult_females$zero_ind$`structural_zero (group = G3 DOX adults)`),])
# ATC G3 A vs. H2O G3 A (0, 143, 279)
ancom_adult_females$res_pair[ancom_adult_females$res_pair$`diff_groupG3 H2O adults`,c("taxon","lfc_groupG3 H2O adults","se_groupG3 H2O adults","q_groupG3 H2O adults")]
dim(ancom_adult_females$zero_ind[ancom_adult_females$zero_ind$`structural_zero (group = G3 ATC adults)` & !(ancom_adult_females$zero_ind$`structural_zero (group = G3 H2O adults)`),])
dim(ancom_adult_females$zero_ind[ancom_adult_females$zero_ind$`structural_zero (group = G3 H2O adults)` & !(ancom_adult_females$zero_ind$`structural_zero (group = G3 ATC adults)`),])
# DOX G4 A vs. H2O G4 A (13, 766, 272)
dim(ancom_adult_females$res_pair[ancom_adult_females$res_pair$`diff_groupG4 H2O adults_groupG4 DOX adults`,c("taxon","lfc_groupG4 H2O adults_groupG4 DOX adults","se_groupG4 H2O adults_groupG4 DOX adults","q_groupG4 H2O adults_groupG4 DOX adults")])
dim(ancom_adult_females$zero_ind[ancom_adult_females$zero_ind$`structural_zero (group = G4 DOX adults)` & !(ancom_adult_females$zero_ind$`structural_zero (group = G4 H2O adults)`),])
dim(ancom_adult_females$zero_ind[ancom_adult_females$zero_ind$`structural_zero (group = G4 H2O adults)` & !(ancom_adult_females$zero_ind$`structural_zero (group = G4 DOX adults)`),])
# ATC G4 A vs. H2O G4 A (1, 537, 234)
ancom_adult_females$res_pair[ancom_adult_females$res_pair$`diff_groupG4 H2O adults_groupG4 ATC adults`,c("taxon","lfc_groupG4 H2O adults_groupG4 ATC adults","se_groupG4 H2O adults_groupG4 ATC adults","q_groupG4 H2O adults_groupG4 ATC adults")]
dim(ancom_adult_females$zero_ind[ancom_adult_females$zero_ind$`structural_zero (group = G4 ATC adults)` & !(ancom_adult_females$zero_ind$`structural_zero (group = G4 H2O adults)`),])
dim(ancom_adult_females$zero_ind[ancom_adult_females$zero_ind$`structural_zero (group = G4 H2O adults)` & !(ancom_adult_females$zero_ind$`structural_zero (group = G4 ATC adults)`),])

# Export differentially abundant taxa
G3DOXL <- ancom_larvae$res_pair[ancom_larvae$res_pair$`diff_groupG3 H2O larvae_groupG3 DOX larvae`,c("taxon","lfc_groupG3 H2O larvae_groupG3 DOX larvae","se_groupG3 H2O larvae_groupG3 DOX larvae","q_groupG3 H2O larvae_groupG3 DOX larvae")]
G3DOXL[,5:8] <- tax_table(ps_final)[ancom_larvae$res_pair[ancom_larvae$res_pair$`diff_groupG3 H2O larvae_groupG3 DOX larvae`,c("taxon","lfc_groupG3 H2O larvae_groupG3 DOX larvae","se_groupG3 H2O larvae_groupG3 DOX larvae","q_groupG3 H2O larvae_groupG3 DOX larvae")][,1],3:6]
write.csv(G3DOXL,"G3DOXL.csv")
G3ATCL <- ancom_larvae$res_pair[ancom_larvae$res_pair$`diff_groupG3 H2O larvae`,c("taxon","lfc_groupG3 H2O larvae","se_groupG3 H2O larvae","q_groupG3 H2O larvae")]
G3ATCL[,5:8] <- tax_table(ps_final)[ancom_larvae$res_pair[ancom_larvae$res_pair$`diff_groupG3 H2O larvae`,c("taxon","lfc_groupG3 H2O larvae","se_groupG3 H2O larvae","q_groupG3 H2O larvae")][,1],3:6]
write.csv(G3ATCL,"G3ATCL.csv")
G4DOXL <- ancom_larvae$res_pair[ancom_larvae$res_pair$`diff_groupG4 H2O larvae_groupG4 DOX larvae`,c("taxon","lfc_groupG4 H2O larvae_groupG4 DOX larvae","se_groupG4 H2O larvae_groupG4 DOX larvae","q_groupG4 H2O larvae_groupG4 DOX larvae")]
G4DOXL[,5:8] <- tax_table(ps_final)[ancom_larvae$res_pair[ancom_larvae$res_pair$`diff_groupG4 H2O larvae_groupG4 DOX larvae`,c("taxon","lfc_groupG4 H2O larvae_groupG4 DOX larvae","se_groupG4 H2O larvae_groupG4 DOX larvae","q_groupG4 H2O larvae_groupG4 DOX larvae")][,1],3:6]
write.csv(G4DOXL,"G4DOXL.csv")
G4DOXA_males <- ancom_adult_males$res_pair[ancom_adult_males$res_pair$`diff_groupG4 H2O adults_groupG4 DOX adults`,c("taxon","lfc_groupG4 H2O adults_groupG4 DOX adults","se_groupG4 H2O adults_groupG4 DOX adults","q_groupG4 H2O adults_groupG4 DOX adults")]
G4DOXA_males[,5:8] <- tax_table(ps_final)[ancom_adult_males$res_pair[ancom_adult_males$res_pair$`diff_groupG4 H2O adults_groupG4 DOX adults`,c("taxon","lfc_groupG4 H2O adults_groupG4 DOX adults","se_groupG4 H2O adults_groupG4 DOX adults","q_groupG4 H2O adults_groupG4 DOX adults")][,1],3:6]
write.csv(G4DOXA_males,"G4DOXA_males.csv")
G4DOXA_females <- ancom_adult_females$res_pair[ancom_adult_females$res_pair$`diff_groupG4 H2O adults_groupG4 DOX adults`,c("taxon","lfc_groupG4 H2O adults_groupG4 DOX adults","se_groupG4 H2O adults_groupG4 DOX adults","q_groupG4 H2O adults_groupG4 DOX adults")]
G4DOXA_females[,5:8] <- tax_table(ps_final)[ancom_adult_females$res_pair[ancom_adult_females$res_pair$`diff_groupG4 H2O adults_groupG4 DOX adults`,c("taxon","lfc_groupG4 H2O adults_groupG4 DOX adults","se_groupG4 H2O adults_groupG4 DOX adults","q_groupG4 H2O adults_groupG4 DOX adults")][,1],3:6]
write.csv(G4DOXA_females,"G4DOXA_females.csv")
G4ATCA_females <- ancom_adult_females$res_pair[ancom_adult_females$res_pair$`diff_groupG4 H2O adults_groupG4 ATC adults`,c("taxon","lfc_groupG4 H2O adults_groupG4 ATC adults","se_groupG4 H2O adults_groupG4 ATC adults","q_groupG4 H2O adults_groupG4 ATC adults")]
G4ATCA_females[,5:8] <- tax_table(ps_final)[ancom_adult_females$res_pair[ancom_adult_females$res_pair$`diff_groupG4 H2O adults_groupG4 ATC adults`,c("taxon","lfc_groupG4 H2O adults_groupG4 ATC adults","se_groupG4 H2O adults_groupG4 ATC adults","q_groupG4 H2O adults_groupG4 ATC adults")][,1],3:6]
write.csv(G4ATCA_females,"G4ATCA_females.csv")

# export ANCOMBC2 results
write.csv(ancom_larvae$res_pair,"ANCOMBC2_larvae_respair.csv")
write.csv(ancom_larvae$zero_ind,"ANCOMBC2_larvae_structzeroes.csv")
write.csv(ancom_adult_males$res_pair,"ANCOMBC2_adult_males_respair.csv")
write.csv(ancom_adult_males$zero_ind,"ANCOMBC2_adult_males_structzeroes.csv")
write.csv(ancom_adult_females$res_pair,"ANCOMBC2_adult_females_respair.csv")
write.csv(ancom_adult_females$zero_ind,"ANCOMBC2_adult_females_structzeroes.csv")

# export final datasets (ASV table, taxonomy table and reference sequences)
write.csv(otu_table(ps_final),"ASV_table.csv")
write.csv(tax_table(ps_final),"taxonomy_table.csv")
write.csv(refseq(ps_final),"refseq.csv")

######## BAR PLOTS FOR DIFFERENTIAL ABUNDANCE TESTING ########

# first, we will agglomerate taxa at genus level, which requires us to normalize the dataset afterwards
ps_glom <- tax_glom(ps_final, taxrank="genus", NArm=FALSE)
SRS.shiny.app(as.data.frame(t(otu_table(ps_glom)))) #Cmin = 65,135 reads
normalized_otu_table_glom <- t(SRS(as.data.frame(t(otu_table(ps_glom))),65135))
colnames(normalized_otu_table_glom) <- colnames(otu_table(ps_glom))
ps_norm_glom <- phyloseq(otu_table(normalized_otu_table_glom, taxa_are_rows=FALSE), sample_data(sample_data), tax_table(tax_table(ps_glom)))

### first, we make a dataframe with proportions of each of the top 12 ASVs in larvae and adults at the genus level

ps_norm_glom_H2OL <- phyloseq::subset_samples(ps_norm_glom, stage=="Larva" & diet=="H2O")
ps_norm_glom_H2OA <- phyloseq::subset_samples(ps_norm_glom, stage=="Adult" & diet=="H2O")
ps_norm_glom_L <- phyloseq::subset_samples(ps_norm_glom, stage=="Larva")
ps_norm_glom_A <- phyloseq::subset_samples(ps_norm_glom, stage=="Adult")

top_10_L <- (otu_table(ps_norm_glom_L)[,rownames(tax_table(tax_sort(ps_norm_glom_H2OL,by=mean))[c(1:9,11),5:6])]/65135)
colnames(top_10_L) <- (tax_table(tax_sort(ps_norm_glom_H2OL,by=mean))[c(1:9,11),6])
top_10_L <- cbind.data.frame(top_10_L,as.vector(sample_data(ps_norm_glom_L)$group))
top_10_L <- cbind.data.frame(top_10_L,as.vector(sample_data(ps_norm_glom_L)$diet))
top_10_L <- cbind.data.frame(top_10_L,as.vector(sample_data(ps_norm_glom_L)$generation))

top_10_A <- (otu_table(ps_norm_glom_A)[,rownames(tax_table(tax_sort(ps_norm_glom_H2OA,by=mean))[c(1:7,9:11),5:6])]/65135)
colnames(top_10_A) <- (tax_table(tax_sort(ps_norm_glom_H2OA,by=mean))[c(1:7,9:11),6])
top_10_A <- cbind.data.frame(top_10_A,as.vector(sample_data(ps_norm_glom_A)$group))
top_10_A <- cbind.data.frame(top_10_A,as.vector(sample_data(ps_norm_glom_A)$diet))
top_10_A <- cbind.data.frame(top_10_A,as.vector(sample_data(ps_norm_glom_A)$generation))
top_10_A <- cbind.data.frame(top_10_A,as.vector(sample_data(ps_norm_glom_A)$fullgroup))

colnames(top_10_L)[11:13] <- c("group","diet","generation")
colnames(top_10_A)[11:14] <- c("group","diet","generation","fullgroup")

write.csv(as.data.frame(top_10_L), file="top10L_genera.csv")

vec <- rep(1,60)
for (x in 1:10) { for (y in 1:6) { vec[(x-1)*6+y] = colnames(top_10_L)[x] }}
vec <- unlist(vec)
vec <- cbind.data.frame(vec,rep(levels(sample_data(ps_norm_glom_L)$group),10))
vec[,3] <- 1:60
for (taxonid in 1:10) {
  averages = aggregate(top_10_L[,taxonid], list(top_10_L$group), FUN=mean)
  for (groupid in 1:6) { vec[(taxonid-1)*6+groupid, 3] <- averages[groupid, 2] }
}

vec_Ag <- rep(1,60)
for (x in 1:10) { for (y in 1:6) { vec_Ag[(x-1)*6+y] = colnames(top_10_A)[x] }}
vec_Ag <- unlist(vec_Ag)
vec_Ag <- cbind.data.frame(vec_Ag,rep(levels(sample_data(ps_norm_glom_A)$group),10))
vec_Ag[,3] <- 1:60
for (taxonid in 1:10) {
  averages = aggregate(top_10_A[,taxonid], list(top_10_A$group), FUN=mean)
  for (groupid in 1:6) { vec_Ag[(taxonid-1)*6+groupid, 3] <- averages[groupid, 2] }
}

colnames(vec) <- c("taxon","group","mean.prop")
colnames(vec_Ag) <- c("taxon","group","mean.prop")

# define a good color palette for plotting
palette_L <- c("navy", "skyblue2", "darkcyan", "aquamarine2", "limegreen", "greenyellow", "khaki1", "orange", "#df4a53","orchid1")
palette_A <- c("skyblue2", "darkcyan", "aquamarine2","khaki1","#df4a53","orchid1","gray80","gray60","gray40","gray20")
palette_all <- c("navy", "skyblue2", "darkcyan", "aquamarine2", "limegreen", "greenyellow", "khaki1", "orange", "#df4a53","orchid1","gray80","gray60","gray40","gray20")
c("Myroides","Morganella","Providencia","Acinetobacter", "Bacteroides", "Vagococcus", "Pseudomonas", "Peptoniphilus", "Staphylococcus", "Sphingobacterium","Mycoplasma","Psychrilyobacter","Undibacterium","Gaiella")

#### plot the genus-level proportional data

# plot larvae
larvae_barplot_genus <- ggplot(data=vec, aes(forcats::fct_relevel(group,"G3 DOX larvae","G4 DOX larvae","G3 ATC larvae","G4 ATC larvae","G3 H2O larvae","G4 H2O larvae"), y=mean.prop, fill=forcats::fct_relevel(taxon,rev(c("Myroides","Morganella","Providencia","Acinetobacter", "Bacteroides", "Vagococcus", "Pseudomonas", "Peptoniphilus", "Staphylococcus", "Sphingobacterium"))))) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=rev(palette_L)) +
  theme(panel.background = element_blank(), axis.line = element_line(size = 0.5, colour = "gray30"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), axis.text.x = element_blank(), strip.text.x= element_text(size=15)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(aspect.ratio = 1.2, legend.position="none")

larvae_barplot_genus

# plot adults (fullgroup)
#adults_fullgroup_barplot_genus <- ggplot(data=vec_A, aes(forcats::fct_relevel(group,"M G3 DOX adults","F G3 DOX adults","M G4 DOX adults","F G4 DOX adults","M G3 ATC adults","F G3 ATC adults","M G4 ATC adults","F G4 ATC adults","M G3 H2O adults","F G3 H2O adults","M G4 H2O adults","F G4 H2O adults"), y=mean.prop, fill=forcats::fct_relevel(taxon,rev(c("Providencia","Pseudomonas","Morganella","Staphylococcus", "Acinetobacter", "Mycoplasma", "Sphingobacterium", "Psychrilyobacter", "Undibacterium", "Gaiella"))))) + 
 # geom_bar(stat="identity") + 
  #scale_fill_manual(values=rev(palette_A)) +
  #theme(panel.background = element_blank(), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.line = element_line(size = 0.5, colour = "gray30"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), axis.text.x = element_blank(), strip.text.x= element_text(size=15)) + 
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  #guides(fill = guide_legend(title = "Genus")) +
  #theme(aspect.ratio = 1.2, legend.position="none")
#adults_fullgroup_barplot_genus

# plot adults (group)
adults_barplot_genus <- ggplot(data=vec_Ag, aes(forcats::fct_relevel(group,"G3 DOX adults","G4 DOX adults","G3 ATC adults","G4 ATC adults","G3 H2O adults","G4 H2O adults"), y=mean.prop, fill=forcats::fct_relevel(taxon,rev(c("Morganella","Providencia","Acinetobacter","Pseudomonas","Staphylococcus","Sphingobacterium","Mycoplasma","Psychrilyobacter","Undibacterium","Gaiella"))))) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=rev(palette_A)) +
  theme(panel.background = element_blank(), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.line = element_line(size = 0.5, colour = "gray30"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), axis.text.x = element_blank(), strip.text.x= element_text(size=15)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  guides(fill = guide_legend(title = "Genus")) +
  theme(aspect.ratio = 1.2, legend.position="none")

adults_barplot_genus

# plot legend

combination_top <- rbind(vec,vec_Ag[c(31:36,43:60),])

legend_genus <- ggplot(data=combination_top, aes(forcats::fct_relevel(group,"G3 DOX adults","G4 DOX adults","G3 ATC adults","G4 ATC adults","G3 H2O adults","G4 H2O adults"), y=mean.prop, fill=forcats::fct_relevel(taxon,rev(c("Myroides","Morganella","Providencia","Acinetobacter", "Bacteroides", "Vagococcus", "Pseudomonas", "Peptoniphilus", "Staphylococcus", "Sphingobacterium","Mycoplasma","Psychrilyobacter","Undibacterium","Gaiella"))))) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=rev(palette_all)) +
  theme(panel.background = element_blank(), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.line = element_line(size = 0.5, colour = "gray30"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), axis.text.x = element_blank(), strip.text.x= element_text(size=15)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  guides(fill = guide_legend(title = "Genus")) +
  theme(aspect.ratio = 0.5, legend.position="left")

legend_genus

#### next, we repeat this at the family level

vecL_g <- vec
vecAg_g <- vecAg

# first, we will agglomerate taxa at family level, which requires us to normalize the dataset afterwards
ps_glom <- tax_glom(ps_final, taxrank="family", NArm=FALSE)
SRS.shiny.app(as.data.frame(t(otu_table(ps_glom)))) #Cmin = 65,135 reads
normalized_otu_table_glom <- t(SRS(as.data.frame(t(otu_table(ps_glom))),65135))
colnames(normalized_otu_table_glom) <- colnames(otu_table(ps_glom))
ps_norm_glom <- phyloseq(otu_table(normalized_otu_table_glom, taxa_are_rows=FALSE), sample_data(sample_data), tax_table(tax_table(ps_glom)))

### first, we make a dataframe with proportions of each of the top 12 ASVs in larvae and adults at the genus level

ps_norm_glom_H2OL <- phyloseq::subset_samples(ps_norm_glom, stage=="Larva" & diet=="H2O")
ps_norm_glom_H2OA <- phyloseq::subset_samples(ps_norm_glom, stage=="Adult" & diet=="H2O")
ps_norm_glom_L <- phyloseq::subset_samples(ps_norm_glom, stage=="Larva")
ps_norm_glom_A <- phyloseq::subset_samples(ps_norm_glom, stage=="Adult")

tax_table(tax_sort(ps_norm_glom_H2OL,by=mean))[1:18,5]
tax_table(tax_sort(ps_norm_glom_H2OA,by=mean))[1:18,5]

top_10_L <- (otu_table(ps_norm_glom_L)[,rownames(tax_table(tax_sort(ps_norm_glom_H2OL,by=mean))[c(1:6,8:11),5])]/65135)
colnames(top_10_L) <- (tax_table(tax_sort(ps_norm_glom_H2OL,by=mean))[c(1:6,8:11),5])
top_10_L <- cbind.data.frame(top_10_L,as.vector(sample_data(ps_norm_glom_L)$group))
top_10_L <- cbind.data.frame(top_10_L,as.vector(sample_data(ps_norm_glom_L)$diet))
top_10_L <- cbind.data.frame(top_10_L,as.vector(sample_data(ps_norm_glom_L)$generation))

top_10_A <- (otu_table(ps_norm_glom_A)[,rownames(tax_table(tax_sort(ps_norm_glom_H2OA,by=mean))[c(1:6,8:11),5])]/65135)
colnames(top_10_A) <- (tax_table(tax_sort(ps_norm_glom_H2OA,by=mean))[c(1:6,8:11),5])
top_10_A <- cbind.data.frame(top_10_A,as.vector(sample_data(ps_norm_glom_A)$group))
top_10_A <- cbind.data.frame(top_10_A,as.vector(sample_data(ps_norm_glom_A)$diet))
top_10_A <- cbind.data.frame(top_10_A,as.vector(sample_data(ps_norm_glom_A)$generation))
top_10_A <- cbind.data.frame(top_10_A,as.vector(sample_data(ps_norm_glom_A)$fullgroup))

colnames(top_10_L)[11:13] <- c("group","diet","generation")
colnames(top_10_A)[11:14] <- c("group","diet","generation","fullgroup")

vec <- rep(1,60)
for (x in 1:10) { for (y in 1:6) { vec[(x-1)*6+y] = colnames(top_10_L)[x] }}
vec <- unlist(vec)
vec <- cbind.data.frame(vec,rep(levels(sample_data(ps_norm_glom_L)$group),10))
vec[,3] <- 1:60
for (taxonid in 1:10) {
  averages = aggregate(top_10_L[,taxonid], list(top_10_L$group), FUN=mean)
  for (groupid in 1:6) { vec[(taxonid-1)*6+groupid, 3] <- averages[groupid, 2] }
}

vec_Ag <- rep(1,60)
for (x in 1:10) { for (y in 1:6) { vec_Ag[(x-1)*6+y] = colnames(top_10_A)[x] }}
vec_Ag <- unlist(vec_Ag)
vec_Ag <- cbind.data.frame(vec_Ag,rep(levels(sample_data(ps_norm_glom_A)$group),10))
vec_Ag[,3] <- 1:60
for (taxonid in 1:10) {
  averages = aggregate(top_10_A[,taxonid], list(top_10_A$group), FUN=mean)
  for (groupid in 1:6) { vec_Ag[(taxonid-1)*6+groupid, 3] <- averages[groupid, 2] }
}

colnames(vec) <- c("taxon","group","mean.prop")
colnames(vec_Ag) <- c("taxon","group","mean.prop")

# define a good color palette for plotting
palette_L <- c("navy","darkcyan","aquamarine2","limegreen","greenyellow","khaki1","#df4a53","bisque2","orchid1","purple")
palette_A <- c("darkcyan","aquamarine2","khaki1","#df4a53","orchid1","purple","gray80","gray40","gray60","sienna")
palette_all <- c("navy","darkcyan","aquamarine2","limegreen","greenyellow","khaki1","#df4a53","orchid1","purple","gray80","gray60","gray40","bisque2","sienna")

# plot larvae
larvae_barplot_family <- ggplot(data=vec, aes(forcats::fct_relevel(group,"G3 DOX larvae","G4 DOX larvae","G3 ATC larvae","G4 ATC larvae","G3 H2O larvae","G4 H2O larvae"), y=mean.prop, fill=forcats::fct_relevel(taxon,rev(c("Flavobacteriaceae","Morganellaceae","Moraxellaceae","Bacteroidaceae","Vagococcaceae","Pseudomonadaceae","Staphylococcaceae","Hafniaceae","Sphingobacteriaceae","Xanthomonadaceae"))))) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=rev(palette_L)) +
  theme(panel.background = element_blank(), axis.line = element_line(size = 0.5, colour = "gray30"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), axis.text.x = element_blank(), strip.text.x= element_text(size=15)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(aspect.ratio = 1.2, legend.position="none")

larvae_barplot_family

# plot adults (group)
adults_barplot_family <- ggplot(data=vec_Ag, aes(forcats::fct_relevel(group,"G3 DOX adults","G4 DOX adults","G3 ATC adults","G4 ATC adults","G3 H2O adults","G4 H2O adults"), y=mean.prop, fill=forcats::fct_relevel(taxon,rev(c("Morganellaceae","Moraxellaceae","Pseudomonadaceae","Staphylococcaceae","Sphingobacteriaceae","Xanthomonadaceae","Mycoplasmataceae","Oxalobacteraceae","Fusobacteriaceae","Nitrosomonadaceae"))))) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=rev(palette_A)) +
  theme(panel.background = element_blank(), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.line = element_line(size = 0.5, colour = "gray30"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), axis.text.x = element_blank(), strip.text.x= element_text(size=15)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  guides(fill = guide_legend(title = "Genus")) +
  theme(aspect.ratio = 1.2, legend.position="none")

adults_barplot_family

# plot legend

combination_top <- rbind(vec,vec_Ag[c(31:42,49:60),])

legend_family <- ggplot(data=combination_top, aes(forcats::fct_relevel(group,"G3 DOX adults","G4 DOX adults","G3 ATC adults","G4 ATC adults","G3 H2O adults","G4 H2O adults"), y=mean.prop, fill=forcats::fct_relevel(taxon,rev(c("Flavobacteriaceae","Morganellaceae","Moraxellaceae","Bacteroidaceae","Vagococcaceae","Pseudomonadaceae","Staphylococcaceae","Sphingobacteriaceae","Hafniaceae","Mycoplasmataceae","Fusobacteriaceae","Oxalobacteraceae","Xanthomonadaceae","Nitrosomonadaceae"))))) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=rev(palette_all)) +
  theme(panel.background = element_blank(), legend.text = element_text(size=8), legend.title = element_text(size=10), axis.line = element_line(size = 0.5, colour = "gray30"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), axis.text.x = element_blank(), strip.text.x= element_text(size=15)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  guides(fill = guide_legend(title = "Family")) +
  theme(aspect.ratio = 0.5, legend.position="left")

legend_family

barplots_vertical <- plot_grid(larvae_barplot_genus, adults_barplot_genus, legend_genus, larvae_barplot_family, adults_barplot_family, legend_family, ncol=3, align="h")
barplots_vertical

#### Plotting the abundance of selected taxa #####

# let's get the taxon abundance comps from ANCOMBC2 into the appropriate format for plotting
props <- otu_table(ps_norm)[1:54,]/65258
props <- cbind(sample_data(ps_norm),props)
v1 <- strsplit(names(ancombc2_allcomps$res_pair[1,277:331]), "_group", fixed = TRUE)
v2 <- sapply(doof, `[[`, 2)
v3 <- sapply(doof, `[[`, 3)
comps <- paste(v2, v3, sep = "-")
comps <- c("G3 H2O larvae-G3 ATC adults","G3 H2O larvae-G3 ATC larvae","G3 H2O larvae-G3 DOX adults","G3 H2O larvae-G3 DOX larvae","G3 H2O larvae-G3 H2O adults","G3 H2O larvae-G4 ATC adults","G3 H2O larvae-G4 ATC larvae","G3 H2O larvae-G4 DOX adults","G3 H2O larvae-G4 DOX larvae","G3 H2O larvae-G4 H2O adults","G3 H2O larvae-G4 H2O larvae",comps)

# for ASV2, etc
ASV28props <- cbind(props[,1:6],props[,"ASV28"])
compsvec_ASV28 <- ancombc2_allcomps$res_pair[ancombc2_allcomps$res_pair$taxon=="ASV28",266:331]
names(compsvec_ASV28) <- comps
letters <- multcompLetters(as.vector(compsvec_ASV28))$Letters
for (x in 1:54) {ASV28props[x, 8] <- letters[which(ASV28props$group[x] == names(letters))]}
for (x in 1:54) {ASV28props[x, 9] <- max(ASV28props[ASV28props$group == ASV28props$group[x],7])}
colnames(ASV28props)[7:9] <- c("prop","letters","maxes")

ASV28plot <- ggplot(ASV28props, aes(forcats::fct_relevel(diet, "DOX", "ATC", "H2O"), prop, group=generation, color=generation)) +
  scale_color_manual(values = c("G3" = "salmon2", "G4" = "#5764ba")) + 
  geom_point(size=2.5,shape=16,alpha=1,position = position_jitterdodge(jitter.width = 0, dodge.width = .7)) + 
  theme(strip.background = element_rect(colour = "black"),panel.border = element_rect(linetype = "solid", fill = NA), panel.background = element_rect(fill = "white", colour = "grey30"), panel.grid.minor=element_blank(), axis.title.x = element_text(size=15),  axis.title.y = element_text(size=15), axis.text.y = element_text(size=18), axis.text.x = element_text(size=12, angle = 0, vjust = 0, hjust=0.5), strip.text.x= element_text(size=20), legend.text = element_text(size=10)) + 
  scale_y_continuous(limits = c(0, 0.01)) + 
  stat_summary(fun.data=mean_se, color=c("salmon2","salmon2","salmon2","#5764ba","#5764ba","#5764ba","salmon2","salmon2","salmon2","#5764ba","#5764ba","#5764ba"),position=position_dodge(width=0.7), size=0.6, width=0.4, geom="errorbar") + 
  stat_summary(fun.data=mean_se, color=c("salmon2","salmon2","salmon2","#5764ba","#5764ba","#5764ba","salmon2","salmon2","salmon2","#5764ba","#5764ba","#5764ba"),geom="point",position=position_dodge(width=0.7), size=2.5, shape=15) + 
  geom_text(aes(label = letters, y = rep(0.01,54)), vjust = -0.5, position=position_dodge(width=0.7),size=4,color="black") + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  facet_wrap2(~forcats::fct_relevel(stage, "Larva", "Adult"), strip = strip)

ASV28plot 

## find out which (and how many) taxa differ significantly between different groups

# G4 ATC adults vs. G3 ATC adults (120 taxa)
taxa_G4ATCAvG3ATCA <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 ATC adults_groupG3 ATC adults"])]
lfc_G4ATCAvG3ATCA <- ancombc2_allcomps$res_pair$"lfc_groupG4 ATC adults_groupG3 ATC adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 ATC adults_groupG3 ATC adults"])]
top30taxa_G4ATCAvG3ATCA <- rev(tail(ancombc2_allcomps$res_pair[order(abs(ancombc2_allcomps$res_pair$"lfc_groupG4 ATC adults_groupG3 ATC adults")),], 30)$taxon)
top30lfc_G4ATCAvG3ATCA <- ancombc2_allcomps$res_pair$"lfc_groupG4 ATC adults_groupG3 ATC adults"[ancombc2_allcomps$res_pair$taxon%in%(top30taxa_G4ATCAvG3ATCA)]
se_G4ATCAvG3ATCA  <- ancombc2_allcomps$res_pair$"se_groupG4 ATC adults_groupG3 ATC adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 ATC adults_groupG3 ATC adults"])]

# G4 DOX adults vs. G3 DOX adults (270 taxa)
taxa_G4DOXAvG3DOXA <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX adults_groupG3 DOX adults"])]
G4vG3_DOX_adult_difftaxa_top10lfc <- rev(tail(ancombc2_allcomps$res_pair[order(abs(ancombc2_allcomps$res_pair$"lfc_groupG4 DOX adults_groupG3 DOX adults")),], 10)$taxon)
# get log fold changes for these differentially abundant taxa
lfc_G4DOXAvG3DOXA <- ancombc2_allcomps$res_pair$"lfc_groupG4 DOX adults_groupG3 DOX adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX adults_groupG3 DOX adults"])]
top30taxa_G4DOXAvG3DOXA <- rev(tail(ancombc2_allcomps$res_pair[order(abs(ancombc2_allcomps$res_pair$"lfc_groupG4 DOX adults_groupG3 DOX adults")),], 30)$taxon)
top30lfc_G4DOXAvG3DOXA <- ancombc2_allcomps$res_pair$"lfc_groupG4 DOX adults_groupG3 DOX adults"[ancombc2_allcomps$res_pair$taxon%in%(top30taxa_G4DOXAvG3DOXA)]
se_G4DOXAvG3DOXA <- ancombc2_allcomps$res_pair$"se_groupG4 DOX adults_groupG3 DOX adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX adults_groupG3 DOX adults"])]

# shared taxa between these 2 groups (69 taxa)
G3vG4_ATC_adult_difftaxa_all[G3vG4_ATC_adult_difftaxa_all %in% G3vG4_DOX_adult_difftaxa_all]
# shared taxa among top 10 lfc (1 taxon, Acinetobacter)
G3vG4_DOX_adult_difftaxa_top10lfc[G3vG4_DOX_adult_difftaxa_top10lfc %in% G3vG4_ATC_adult_difftaxa_top10lfc]

# G4 H2O adults vs. G3 H2O adults (5 taxa)
ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O adults_groupG3 H2O adults"])]
taxa_G4H2OAvG3H2OA <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O adults_groupG3 H2O adults"])]
lfc_G4H2OAvG3H2OA <- ancombc2_allcomps$res_pair$"lfc_groupG4 H2O adults_groupG3 H2O adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O adults_groupG3 H2O adults"])]
se_G4H2OAvG3H2OA  <- ancombc2_allcomps$res_pair$"se_groupG4 H2O adults_groupG3 H2O adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O adults_groupG3 H2O adults"])]


# G4 DOX larvae vs. G3 DOX larvae (2 taxa)
taxa_G4DOXLvG3DOXL <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX larvae_groupG3 DOX larvae"])]
lfc_G4DOXLvG3DOXL <- ancombc2_allcomps$res_pair$"lfc_groupG4 DOX larvae_groupG3 DOX larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX larvae_groupG3 DOX larvae"])]
se_G4DOXLvG3DOXL   <- ancombc2_allcomps$res_pair$"se_groupG4 DOX larvae_groupG3 DOX larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX larvae_groupG3 DOX larvae"])]

# G4 ATC larvae vs. G3 ATC larvae (4 taxa)
taxa_G4ATCLvG3ATCL <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 ATC larvae_groupG3 ATC larvae"])]
lfc_G4ATCLvG3ATCL <- ancombc2_allcomps$res_pair$"lfc_groupG4 ATC larvae_groupG3 ATC larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 ATC larvae_groupG3 ATC larvae"])]
se_G4ATCLvG3ATCL   <- ancombc2_allcomps$res_pair$"se_groupG4 ATC larvae_groupG3 ATC larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 ATC larvae_groupG3 ATC larvae"])]

# G4 H2O larvae vs. G3 H2O larvae (3 taxa)
taxa_G4H2OLvG3H2OL <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O larvae"])]
lfc_G4H2OLvG3H2OL <- ancombc2_allcomps$res_pair$"lfc_groupG4 H2O larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O larvae"])]
se_G4H2OLvG3H2OL   <- ancombc2_allcomps$res_pair$"se_groupG4 H2O larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O larvae"])]

#### G3 vs. G3 comparisons ####

# G3 DOX larvae vs. G3 H2O larvae (1 taxon)
taxa_G3DOXLvG3H2OL <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG3 DOX larvae"])]
lfc_G3DOXLvG3H2OL <- ancombc2_allcomps$res_pair$"lfc_groupG3 DOX larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 DOX larvae"])]
se_G3DOXLvG3H2OL <- ancombc2_allcomps$res_pair$"se_groupG3 DOX larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 DOX larvae"])]

# G3 ATC larvae vs. G3 H2O larvae (2 taxa)
taxa_G3ATCLvG3H2OL <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG3 ATC larvae"])]
lfc_G3ATCLvG3H2OL <- ancombc2_allcomps$res_pair$"lfc_groupG3 ATC larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 ATC larvae"])]
se_G3ATCLvG3H2OL <- ancombc2_allcomps$res_pair$"se_groupG3 ATC larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 ATC larvae"])]

# G3 DOX larvae vs. G3 ATC larvae (0 taxa)
taxa_G3DOXLvG3ATCL <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG3 DOX larvae_groupG3 ATC larvae"])]
lfc_G3DOXLvG3ATCL <- ancombc2_allcomps$res_pair$"lfc_groupG3 DOX larvae_groupG3 ATC larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 DOX larvae_groupG3 ATC larvae"])]
se_G3DOXLvG3ATCL <- ancombc2_allcomps$res_pair$"se_groupG3 DOX larvae_groupG3 ATC larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 DOX larvae_groupG3 ATC larvae"])]

# G3 H2O adults vs. G3 DOX adults (0 taxa)
taxa_G3DOXAvG3H2OA <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG3 H2O adults_groupG3 DOX adults"])]
lfc_G3DOXAvG3H2OA <- ancombc2_allcomps$res_pair$"lfc_diff_groupG3 H2O adults_groupG3 DOX adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 H2O adults_groupG3 DOX adults"])]
se_G3DOXAvG3H2OA <- ancombc2_allcomps$res_pair$"se_diff_groupG3 H2O adults_groupG3 DOX adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 H2O adults_groupG3 DOX adults"])]

# G3 H2O adults vs. G3 ATC adults (0 taxa)
taxa_G3ATCAvG3H2OA <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG3 H2O adults_groupG3 ATC adults"])]
lfc_G3ATCAvG3H2OA <- ancombc2_allcomps$res_pair$"lfc_groupG3 H2O adults_groupG3 ATC adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 H2O adults_groupG3 ATC adults"])]
se_G3ATCAvG3H2OA <- ancombc2_allcomps$res_pair$"se_groupG3 H2O adults_groupG3 ATC adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 H2O adults_groupG3 ATC adults"])]

# G3 DOX adults vs. G3 ATC adults (1 taxon)
taxa_G3DOXAvG3ATCA <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG3 DOX adults_groupG3 ATC adults"])]
lfc_G3DOXAvG3ATCA <- ancombc2_allcomps$res_pair$"lfc_groupG3 DOX adults_groupG3 ATC adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 DOX adults_groupG3 ATC adults"])]
se_G3DOXAvG3ATCA <- ancombc2_allcomps$res_pair$"se_groupG3 DOX adults_groupG3 ATC adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG3 DOX adults_groupG3 ATC adults"])]

##### G4 vs. G4 comparisons #####

# G4 DOX larvae vs. G4 H2O larvae (1 taxon)
taxa_G4DOXLvG4H2OL <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O larvae_groupG4 DOX larvae"])]
lfc_G4DOXLvG4H2OL <- -ancombc2_allcomps$res_pair$"lfc_groupG4 H2O larvae_groupG4 DOX larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O larvae_groupG4 DOX larvae"])]
se_G4DOXLvG4H2OL <- ancombc2_allcomps$res_pair$"se_groupG4 H2O larvae_groupG4 DOX larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O larvae_groupG4 DOX larvae"])]

# G4 ATC larvae vs. G4 H2O larvae (0 taxa)
taxa_G4ATCLvG4H2OL <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O larvae_groupG4 ATC larvae"])]
lfc_G4ATCLvG4H2OL <- ancombc2_allcomps$res_pair$"lfc_groupG4 H2O larvae_groupG4 ATC larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O larvae_groupG4 ATC larvae"])]
se_G4ATCLvG4H2OL <- ancombc2_allcomps$res_pair$"se_groupG4 H2O larvae_groupG4 ATC larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O larvae_groupG4 ATC larvae"])]

# G4 DOX larvae vs. G4 ATC larvae (2 taxa)
taxa_G4DOXLvG4ATCL <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX larvae_groupG4 ATC larvae"])]
lfc_G4DOXLvG4ATCL <- ancombc2_allcomps$res_pair$"lfc_groupG4 DOX larvae_groupG4 ATC larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX larvae_groupG4 ATC larvae"])]
se_G4DOXLvG4ATCL <- ancombc2_allcomps$res_pair$"se_groupG4 DOX larvae_groupG4 ATC larvae"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX larvae_groupG4 ATC larvae"])]

# G4 DOX adults vs. G4 H2O adults (191 taxa)
taxa_G4DOXAvG4H2OA <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O adults_groupG4 DOX adults"])]
lfc_G4DOXAvG4H2OA <- -ancombc2_allcomps$res_pair$"lfc_groupG4 H2O adults_groupG4 DOX adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O adults_groupG4 DOX adults"])]
se_G4DOXAvG4H2OA <- ancombc2_allcomps$res_pair$"se_groupG4 H2O adults_groupG4 DOX adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O adults_groupG4 DOX adults"])]

# G4 ATC adults vs. G4 H2O adults (13 taxa)
taxa_G4ATCAvG4H2OA <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O adults_groupG4 ATC adults"])]
lfc_G4ATCAvG4H2OA <- -ancombc2_allcomps$res_pair$"lfc_groupG4 H2O adults_groupG4 ATC adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O adults_groupG4 ATC adults"])]
se_G4ATCAvG4H2OA <- ancombc2_allcomps$res_pair$"se_groupG4 H2O adults_groupG4 ATC adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 H2O adults_groupG4 ATC adults"])]

# G4 DOX adults vs. G4 ATC adults (111 taxa)
taxa_G4DOXAvG4ATCA <- ancombc2_allcomps$res_pair$taxon[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX adults_groupG4 ATC adults"])]
lfc_G4DOXAvG4ATCA <- ancombc2_allcomps$res_pair$"lfc_groupG4 DOX adults_groupG4 ATC adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX adults_groupG4 ATC adults"])]
se_G4DOXAvG4ATCA <- ancombc2_allcomps$res_pair$"se_groupG4 DOX adults_groupG4 ATC adults"[which(ancombc2_allcomps$res_pair[,"diff_groupG4 DOX adults_groupG4 ATC adults"])]


## takeaways: 2 ASVs are different between G3vG3 ATC larvae, and G3vG4 H2O larvae, both Acinetobacter
## these are probably temporal effect, good to have control group!

# make the G3 vs. G3 plot
G3vG3 <- data_frame(taxa=c(taxa_G3DOXLvG3H2OL,taxa_G3ATCLvG3H2OL,taxa_G3DOXLvG3ATCL,taxa_G3DOXAvG3H2OA,taxa_G3ATCAvG3H2OA,taxa_G3DOXAvG3ATCA),lfc=c(lfc_G3DOXLvG3H2OL,lfc_G3ATCLvG3H2OL,lfc_G3DOXLvG3ATCL,lfc_G3DOXAvG3H2OA,lfc_G3ATCAvG3H2OA,lfc_G3DOXAvG3ATCA),comp=c(rep("G3DOXLvG3H2OL",length(taxa_G3DOXLvG3H2OL)),rep("G3ATCLvG3H2OL",length(taxa_G3ATCLvG3H2OL)),rep("G3DOXLvG3ATCL",length(taxa_G3DOXLvG3ATCL)),rep("G3DOXAvG3H2OA",length(taxa_G3DOXAvG3H2OA)),rep("G3ATCAvG3H2OA",length(taxa_G3ATCAvG3H2OA)),rep("G3DOXAvG3ATCA",length(taxa_G3DOXAvG3ATCA))))
G3vG3 <- as.data.frame(G3vG3)
G3vG3 <- G3vG3[order(G3vG3[,3],-G3vG3[,2]),]
G3vG3[,4] <- 1:length(G3vG3[,1])
G3vG3[,5] <- 1:length(G3vG3[,1])
for (x in 1:length(G3vG3[,4])) {G3vG3[x,5] = tax_table(ps_final)[G3vG3[x,"taxa"],"genus"]}
for (x in 1:length(G3vG3[,4])) {if (is.na(G3vG3[x,5])) {G3vG3[x,5] = tax_table(ps_final)[G3vG3[x,"taxa"],"family"]}}
anyNA(G3vG3[,5])

G3vG3[4,5] <- "Stenotrophomonas (1)"

G3vG3[,6] <- c(se_G3DOXLvG3H2OL,se_G3ATCLvG3H2OL,se_G3DOXLvG3ATCL,se_G3DOXAvG3H2OA,se_G3ATCAvG3H2OA,se_G3DOXAvG3ATCA)
for (x in 1:length(G3vG3[,1])) {if (G3vG3[x,2] < 0) {G3vG3[x,6] = G3vG3[x,6] * -1}}
for (x in 1:length(G3vG3[,1])) {if (G3vG3[x,2] > 0) {G3vG3[x,6] = abs(G3vG3[x,6])}}

colnames(G3vG3)[4:6] <- c("order","ID","se")
orders <- ordered(G3vG3$comp, levels = c("G3DOXLvG3H2OL", "G3ATCLvG3H2OL", "G3DOXLvG3ATCL","G3DOXAvG3H2OA","G3ATCAvG3H2OA","G3DOXAvG3ATCA"))
G3vG3 <- G3vG3[order(orders),]

# make the G4 vs. G4 plot

G4vG4 <- data_frame(taxa=c(taxa_G4DOXLvG4H2OL,taxa_G4ATCLvG4H2OL,taxa_G4DOXLvG4ATCL,taxa_G4DOXAvG4H2OA,taxa_G4ATCAvG4H2OA,taxa_G4DOXAvG4ATCA),lfc=c(lfc_G4DOXLvG4H2OL,lfc_G4ATCLvG4H2OL,lfc_G4DOXLvG4ATCL,lfc_G4DOXAvG4H2OA,lfc_G4ATCAvG4H2OA,lfc_G4DOXAvG4ATCA),comp=c(rep("G4DOXLvG4H2OL",length(taxa_G4DOXLvG4H2OL)),rep("G4ATCLvG4H2OL",length(taxa_G4ATCLvG4H2OL)),rep("G4DOXLvG4ATCL",length(taxa_G4DOXLvG4ATCL)),rep("G4DOXAvG4H2OA",length(taxa_G4DOXAvG4H2OA)),rep("G4ATCAvG4H2OA",length(taxa_G4ATCAvG4H2OA)),rep("G4DOXAvG4ATCA",length(taxa_G4DOXAvG4ATCA))))
G4vG4 <- as.data.frame(G4vG4)
G4vG4 <- G4vG4[order(G4vG4[,3],-G4vG4[,2]),]
G4vG4[,4] <- 1:length(G4vG4[,1])
G4vG4[,5] <- 1:length(G4vG4[,1])
for (x in 1:length(G4vG4[,4])) {G4vG4[x,5] = tax_table(ps_final)[G4vG4[x,"taxa"],"genus"]}
for (x in 1:length(G4vG4[,4])) {if (is.na(G4vG4[x,5])) {G4vG4[x,5] = tax_table(ps_final)[G4vG4[x,"taxa"],"family"]}}
anyNA(G4vG4[,5])
for (x in 1:length(G4vG4[,4])) {if (is.na(G4vG4[x,5])) {G4vG4[x,5] = tax_table(ps_final)[G4vG4[x,"taxa"],"order"]}}
anyNA(G4vG4[,5])
for (x in 1:length(G4vG4[,4])) {if (is.na(G4vG4[x,5])) {G4vG4[x,5] = tax_table(ps_final)[G4vG4[x,"taxa"],"class"]}}
anyNA(G4vG4[,5])
for (x in 1:length(G4vG4[,4])) {if (is.na(G4vG4[x,5])) {G4vG4[x,5] = tax_table(ps_final)[G4vG4[x,"taxa"],"phylum"]}}
anyNA(G4vG4[,5])
for (x in 1:length(G4vG4[,4])) {if (is.na(G4vG4[x,5])) {G4vG4[x,5] = tax_table(ps_final)[G4vG4[x,"taxa"],"domain"]}}
anyNA(G4vG4[,5])

G4vG4[,6] <- c(se_G4DOXLvG4H2OL,se_G4ATCLvG4H2OL,se_G4DOXLvG4ATCL,se_G4DOXAvG4H2OA,se_G4ATCAvG4H2OA,se_G4DOXAvG4ATCA)
for (x in 1:length(G4vG4[,1])) {if (G4vG4[x,2] < 0) {G4vG4[x,6] = G4vG4[x,6] * -1}}
for (x in 1:length(G4vG4[,1])) {if (G4vG4[x,2] > 0) {G4vG4[x,6] = abs(G4vG4[x,6])}}

colnames(G4vG4)[4:6] <- c("order","ID","se")
orders <- ordered(G4vG4$comp, levels = c("G4DOXLvG4H2OL", "G4ATCLvG4H2OL", "G4DOXLvG4ATCL","G4DOXAvG4H2OA","G4ATCAvG4H2OA","G4DOXAvG4ATCA"))
G4vG4 <- G4vG4[order(orders),]

G4vG4small <- subset(G4vG4, G4vG4$comp %in% c("G4DOXLvG4H2OL","G4ATCLvG4H2OL","G4DOXLvG4ATCL","G4ATCAvG4H2OA"))
orders <- ordered(G4vG4small$comp, levels = c("G4DOXLvG4H2OL","G4ATCLvG4H2OL","G4DOXLvG4ATCL","G4ATCAvG4H2OA"))
G4vG4small <- G4vG4small[order(orders),]
G4vG4small[c(1,3,5,6,14),5] <- c("Myroides (1)","Myroides (2)","Holophagae Subgroup 7","Xanthomonadaceae (1)","Xanthomonadaceae (2)")
G4vG4big <- subset(G4vG4, G4vG4$comp %in% c("G4DOXAvG4H2OA","G4DOXAvG4ATCA"))
orders <- ordered(G4vG4big$comp, levels = c("G4DOXAvG4H2OA","G4DOXAvG4ATCA"))
G4vG4big <- G4vG4big[order(orders),]

G4vG3 <- data_frame(taxa=c(taxa_G4DOXLvG3DOXL,taxa_G4ATCLvG3ATCL,taxa_G4H2OLvG3H2OL,taxa_G4DOXAvG3DOXA,taxa_G4ATCAvG3ATCA,taxa_G4H2OAvG3H2OA),lfc=c(lfc_G4DOXLvG3DOXL,lfc_G4ATCLvG3ATCL,lfc_G4H2OLvG3H2OL,lfc_G4DOXAvG3DOXA,lfc_G4ATCAvG3ATCA,lfc_G4H2OAvG3H2OA),comp=c(rep("G4DOXLvG3DOXL",length(taxa_G4DOXLvG3DOXL)),rep("G4ATCLvG3ATCL",length(taxa_G4ATCLvG3ATCL)),rep("G4H2OLvG3H2OL",length(taxa_G4H2OLvG3H2OL)),rep("G4DOXAvG3DOXA",length(taxa_G4DOXAvG3DOXA)),rep("G4ATCAvG3ATCA",length(taxa_G4ATCAvG3ATCA)),rep("G4H2OAvG3H2OA",length(taxa_G4H2OAvG3H2OA))))
G4vG3 <- as.data.frame(G4vG3)
G4vG3 <- G4vG3[order(G4vG3[,3],-G4vG3[,2]),]
G4vG3[,4] <- 1:length(G4vG3[,1])
G4vG3[,5] <- 1:length(G4vG3[,1])
for (x in 1:length(G4vG3[,4])) {G4vG3[x,5] = tax_table(ps_final)[G4vG3[x,"taxa"],"genus"]}
anyNA(G4vG3[,5])
for (x in 1:length(G4vG3[,4])) {if (is.na(G4vG3[x,5])) {G4vG3[x,5] = tax_table(ps_final)[G4vG3[x,"taxa"],"family"]}}
anyNA(G4vG3[,5])
for (x in 1:length(G4vG3[,4])) {if (is.na(G4vG3[x,5])) {G4vG3[x,5] = tax_table(ps_final)[G4vG3[x,"taxa"],"order"]}}
anyNA(G4vG3[,5])
for (x in 1:length(G4vG3[,4])) {if (is.na(G4vG3[x,5])) {G4vG3[x,5] = tax_table(ps_final)[G4vG3[x,"taxa"],"class"]}}
anyNA(G4vG3[,5])
for (x in 1:length(G4vG3[,4])) {if (is.na(G4vG3[x,5])) {G4vG3[x,5] = tax_table(ps_final)[G4vG3[x,"taxa"],"phylum"]}}
anyNA(G4vG3[,5])
for (x in 1:length(G4vG3[,4])) {if (is.na(G4vG3[x,5])) {G4vG3[x,5] = tax_table(ps_final)[G4vG3[x,"taxa"],"domain"]}}
anyNA(G4vG3[,5])

G4vG3[,6] <- c(se_G4DOXLvG3DOXL,se_G4ATCLvG3ATCL,se_G4H2OLvG3H2OL,se_G4DOXAvG3DOXA,se_G4ATCAvG3ATCA,se_G4H2OAvG3H2OA)
for (x in 1:length(G4vG3[,1])) {if (G4vG3[x,2] < 0) {G4vG3[x,6] = G4vG3[x,6] * -1}}
for (x in 1:length(G4vG3[,1])) {if (G4vG3[x,2] > 0) {G4vG3[x,6] = abs(G4vG3[x,6])}}

colnames(G4vG3)[4:6] <- c("order","ID","se")
orders <- ordered(G4vG3$comp, levels = c("G4H2OLvG3H2OL", "G4DOXLvG3DOXL", "G4ATCLvG3ATCL","G4H2OAvG3H2OA","G4DOXAvG3DOXA","G4ATCAvG3ATCA"))
G4vG3 <- G4vG3[order(orders),]

G4vG3small <- subset(G4vG3, G4vG3$comp %in% c("G4DOXLvG3DOXL","G4ATCLvG3ATCL","G4H2OLvG3H2OL","G4H2OAvG3H2OA"))
orders <- ordered(G4vG3small$comp, levels = c("G4H2OLvG3H2OL", "G4DOXLvG3DOXL", "G4ATCLvG3ATCL","G4H2OAvG3H2OA"))
G4vG3small <- G4vG3small[order(orders),]
G4vG3small[1:10,5] <- c("Lactobacillales (2)","Acinetobacter (1)", "Acinetobacter (2)", "Vagococcus (2)","Myroides (1)","Vagococcus (3)","Acinetobacter (1)","Myroides (2)","Acinetobacter (2)", "Pseudomonas (3)")
G4vG3big <- subset(G4vG3, G4vG3$comp %in% c("G4DOXAvG3DOXA","G4ATCAvG3ATCA"))
orders <- ordered(G4vG3big$comp, levels = c("G4DOXAvG3DOXA","G4ATCAvG3ATCA"))
G4vG3big <- G4vG3big[order(orders),]

G4vG4big[c(4,187),5] <- c("Lactobacillales (1)","Vagococcus (1)")
top12_G4DOXAvG4H2OA <- G4vG4big[c(1:6,186:191),]
G4vG4big[c(196,298,300,301,302),5] <- c("Enterobacterales (1)","Vagococcus (1)", "Myroides (2)", "Pseudomonas (1)","Pseudomonas (2)")
top12_G4DOXAvG4ATCA <- G4vG4big[c(192:197,297:302),]
G4vG3big[c(2,3,267,268,269,270),5] <- c("Enterobacterales (1)", "Enterobacterales (2)","Myroides (2)", "Vagococcus (1)","Acinetobacter (3)","Stenotrophomonas (2)")
top12_G4DOXAvG3DOXA <- G4vG3big[c(1:6,265:270),]
G4vG3big[c(271,274,388,390),5] <- c("Pseudomonas (1)","Pseudomonas (3)","Xanthomonadaceae (2)","Acinetobacter (3)")
top12_G4ATCAvG3ATCA <- G4vG3big[c(271:276,385:390),]

taxa59 <- sort(unique(c(G3vG3[,5],G4vG4small[,5],G4vG3small[,5],G4vG4big[c(1:6,186:191),5],G4vG4big[c(192:197,297:302),5],G4vG3big[c(1:6,265:270),5],G4vG3big[c(271:276,385:390),5])))
taxa59key <- as.data.frame(data_frame(ID=taxa59, colors=c("white","#FFECFA","#DCD6FE","aquamarine2","aquamarine2","aquamarine2","#9FD47B","#FFECFA","#FFECFA","#FFECFA","#9FD47B","#DCD6FE","aquamarine2","#FFECFA","#9FD47B","aquamarine2","aquamarine2","aquamarine2","aquamarine2","#F2DE4D","white","#FFC18F","#DCD6FE","#F2DE4D","#FFECFA","#FFECFA","#FFECFA","aquamarine2","#F2DE4D","aquamarine2","#FFECFA","#9FD47B","#9FD47B","#F2DE4D","#FFECFA","aquamarine2","#9FD47B","#9FD47B","aquamarine2","aquamarine2","aquamarine2","#F2DE4D","#FFC18F","aquamarine2","#FFECFA","#FFECFA","aquamarine2","aquamarine2","aquamarine2","#FFECFA","#FFECFA","white","#FFECFA","#FFECFA","#FFECFA","#FFECFA","white","aquamarine2","aquamarine2")))
taxa59key$colors <- sub("aquamarine2", "#F2DE2D", taxa59key$colors)
taxa59key$colors <- sub("#F2DE4D", "#FFAEE8", taxa59key$colors)
taxa59key$colors <- sub("#FFECFA", "#8FdfFB", taxa59key$colors)
taxa59key$colors <- sub("#DCD6FE", "#DCC1FD", taxa59key$colors)
taxa59key$colors <- sub("#9FD47B", "#DAFF8F", taxa59key$colors)

#colors <- c(rep("orchid3",7),rep("#df4a53",6),rep("orange",3),rep("khaki1",3),rep("greenyellow",2),rep("aquamarine2",1),rep("darkcyan",1),rep("skyblue2",1),rep("navy",3),rep("purple3",6),rep("#FCDCFF",1),rep("#D9FFD2",7),rep("#FFDCCA",1),rep("#B7C6FF",6),rep("#BA9C06",2),rep("#6E3636",2),rep("gray90",4),rep("gray50",1),rep("gray20",2))
#patternsvec <- c("none","stripe","circle","wave","crosshatch","stripe","wave")
#patterns <- c(patternsvec[1:7],patternsvec[1:6],patternsvec[1:3],patternsvec[1:3],patternsvec[1:2],patternsvec[1],patternsvec[1],patternsvec[1],patternsvec[1:3],patternsvec[1:6],patternsvec[1],patternsvec[1:7],patternsvec[1],patternsvec[1:6],patternsvec[1:2],patternsvec[1:2],patternsvec[1:4],patternsvec[1],patternsvec[1:2])
#patterns_fillvec <- c("orchid3","black","black","black","black","black","black","#df4a53","black","black","black","black","black","orange","black","black","khaki1","black","black","greenyellow","black","aquamarine2","darkcyan","skyblue2","navy","white","white","purple3","white","white","white","white","white","#FCDCFF","#D9FFD2","black","black","black","black","black","black","#FFDCCA","#B7C6FF","black","black","black","black","black","#BA9C06","black","#6E3636","white","gray90","black","black","black","gray50","gray20","white")
#anglesvec <- c(30,30,30,30,30,120,30)
#angles <- c(anglesvec[1:7],anglesvec[1:6],anglesvec[1:3],anglesvec[1:3],anglesvec[1:2],anglesvec[1],anglesvec[1],anglesvec[1],anglesvec[1:3],anglesvec[1:6],anglesvec[1],anglesvec[1:7],anglesvec[1],anglesvec[1:6],anglesvec[1:2],anglesvec[1:2],anglesvec[1:4],anglesvec[1],anglesvec[1:2])
#density_vec <- c(1, 0.3, 0.3, 0.15, 0.1, 0.3, 0.2)
#densities <- c(density_vec[1:7],density_vec[1:6],density_vec[1:3],density_vec[1:3],density_vec[1:2],density_vec[1],density_vec[1],density_vec[1],density_vec[1:3],density_vec[1:6],density_vec[1],density_vec[1:7],density_vec[1],density_vec[1:6],density_vec[1:2],density_vec[1:2],density_vec[1:4],density_vec[1],density_vec[1:2])
#type_vec <- c(rep(NA,6),"sine")
#types <- c(type_vec[1:7],type_vec[1:6],type_vec[1:3],type_vec[1:3],type_vec[1:2],type_vec[1],type_vec[1],type_vec[1],type_vec[1:3],type_vec[1:6],type_vec[1],type_vec[1:7],type_vec[1],type_vec[1:6],type_vec[1:2],type_vec[1:2],type_vec[1:4],type_vec[1],type_vec[1:2])
#taxa59 <- sort(unique(c(G3vG3[,5],G4vG4small[,5],G4vG3small[,5],G4vG4big[c(1:6,186:191),5],G4vG4big[c(192:197,297:302),5],G4vG3big[c(1:6,265:270),5],G4vG3big[c(271:276,385:390),5])))
#taxa59key <- as.data.frame(data_frame(ID=taxa59, colors=colors, patterns=patterns,patternsfill=patterns_fillvec,angles=angles, densities=densities, types=types))



G3vG3merged <- merge(taxa59key, G3vG3, by="ID")
G4vG4smallmerged  <- merge(taxa59key, G4vG4small, by="ID")
G4vG3smallmerged  <- merge(taxa59key, G4vG3small, by="ID")
top12_G4DOXAvG4H2OAmerged  <- merge(taxa59key, top12_G4DOXAvG4H2OA, by="ID")
top12_G4DOXAvG4ATCAmerged  <- merge(taxa59key, top12_G4DOXAvG4ATCA, by="ID")
top12_G4DOXAvG3DOXAmerged  <- merge(taxa59key, top12_G4DOXAvG3DOXA, by="ID")
top12_G4ATCAvG3ATCAmerged  <- merge(taxa59key, top12_G4ATCAvG3ATCA, by="ID")

combo1 <- rbind(G3vG3merged,G4vG4smallmerged,top12_G4DOXAvG4H2OAmerged,top12_G4DOXAvG4ATCAmerged)
combo1$order[5:44] <- combo1$order[5:44]+4
combo1 <- arrange(combo1, factor(comp, levels=c("G3DOXLvG3H2OL","G3ATCLvG3H2OL","G3DOXAvG3ATCA","G4DOXLvG4H2OL","G4ATCLvG4H2OL","G4DOXLvG4ATCL","G4ATCAvG4H2OA","G4DOXAvG4H2OA","G4DOXAvG4ATCA")))

comboplot1 <- ggplot(combo1,aes(factor(order),lfc)) + 
  geom_errorbar(aes(ymin=lfc, ymax=lfc+se), width=c(0.46,0.81,0.81,0.46,0.46,rep(0.8,39)),color="gray20") +
  theme(plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill="gray99", colour = "black"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.y = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.y=element_line(), axis.ticks.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="transparent")) +
  geom_bar(stat="identity",fill="gray20",color="gray20",width=c(0.46,0.81,0.81,0.46,0.46,rep(0.8,39))) +
  scale_y_continuous(limits = c(-7.8, 7.8),breaks=seq(-7,7,by=1),expand=c(0,0)) +
  facet_grid(~factor(comp,levels=c("G3DOXLvG3H2OL","G3ATCLvG3H2OL","G3DOXAvG3ATCA","G4DOXLvG4H2OL","G4ATCLvG4H2OL","G4DOXLvG4ATCL","G4ATCAvG4H2OA","G4DOXAvG4H2OA","G4DOXAvG4ATCA")),scales="free_x",space="free_x") +
  labs(y="log fold change") +
  geom_text(aes(label=ID,color=factor(order)),color="white",size=c(rep(2,8),1.6,1.4,1.6,2,1.6,2,1.6,1.6,1.6,2,1.6,1.6,rep(2,24)),position=position_stack(vjust=0.5),vjust=0.4, angle=90,fontface="bold") +
  geom_hline(yintercept = 0,colour="gray60") +
  force_panelsizes(rows=rep(1,8),cols=c(0.21,0.21,0.21,0.21,0.21,1.3,1.2,1.2),respect=FALSE)
comboplot1

#comboplot1 <- ggplot(combo1,aes(factor(order),lfc,color=ID)) + 
#theme(plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill="white", colour = "black"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.y = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.y=element_line(), axis.ticks.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="transparent")) +
#geom_bar_pattern(aes(fill=factor(order),pattern=factor(order),pattern_fill=factor(order),pattern_color=factor(order),color=factor(order),pattern_angle=factor(order),pattern_type=factor(order),pattern_density=factor(order),pattern_spacing=factor(order),pattern_yoffset=factor(order)),color="black",fill=combo1$colors, pattern=combo1$patterns, pattern_fill=combo1$patternsfill, pattern_color=combo1$patternsfill,pattern_angle=as.numeric(combo1$angles), pattern_type=combo1$types, pattern_density=as.numeric(combo1$densities), pattern_spacing=as.numeric(combo1$spacing), pattern_yoffset=as.numeric(combo1$y.offset), stat="identity",width=c(0.45,0.8,0.8,0.45,0.45,rep(0.8,39))) +
#scale_y_continuous(limits = c(-7.8, 7.8),breaks=seq(-7,7,by=1),expand=c(0,0)) +
#facet_grid(~factor(comp,levels=c("G3DOXLvG3H2OL","G3ATCLvG3H2OL","G3DOXAvG3ATCA","G4DOXLvG4H2OL","G4ATCLvG4H2OL","G4DOXLvG4ATCL","G4ATCAvG4H2OA","G4DOXAvG4H2OA","G4DOXAvG4ATCA")),scales="free_x",space="free_x") +
#labs(y="log fold change") +
#geom_text(aes(label=ID),color="black",position=position_stack(vjust=0.5), angle=90) +
#geom_hline(yintercept = 0,colour="gray60") +
# geom_errorbar(aes(ymin=lfc, ymax=lfc+se), width=c(0.45,0.8,0.8,0.45,0.45,rep(0.8,39)),color="black") +
#force_panelsizes(rows=rep(1,8),cols=c(0.2,0.2,0.2,0.2,0.2,1.3,1.2,1.2),respect=FALSE)

comboplot1

combo2 <- rbind(G4vG3smallmerged,top12_G4DOXAvG3DOXAmerged,top12_G4ATCAvG3ATCAmerged)
combo2 <- arrange(combo2, factor(comp, levels=c("G4H2OLvG3H2OL","G4DOXLvG3DOXL","G4ATCLvG3ATCL","G4H2OAvG3H2OA","G4DOXAvG3DOXA","G4ATCAvG3ATCA")))

comboplot2 <- ggplot(combo2,aes(factor(order),lfc)) + 
  theme(plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill="white", colour = "black"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.y = element_text(size=12), axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.y=element_line(), axis.ticks.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="transparent")) +
  geom_errorbar(aes(ymin=lfc, ymax=lfc+se), width=0.8,color="gray20") +
  geom_bar(stat="identity",fill="gray20",color="gray20",width=0.8) +
  scale_y_continuous(limits = c(-7.8, 7.8),breaks=seq(-7,7,by=1),expand=c(0,0)) +
  facet_grid(~factor(comp,levels=c("G4H2OLvG3H2OL","G4DOXLvG3DOXL","G4ATCLvG3ATCL","G4H2OAvG3H2OA","G4DOXAvG3DOXA","G4ATCAvG3ATCA")),scales="free_x",space="free_x") +
  labs(y="log fold change") +
  geom_text(aes(label=ID,color=factor(order)),color="white",size=rep(2,38),position=position_stack(vjust=0.5),vjust=0.4, angle=90,fontface="bold") +
  geom_hline(yintercept = 0,colour="gray60") +
  force_panelsizes(rows=rep(1,6),cols=c(0.3,0.2,0.4,0.5,1.2,1.2),respect=FALSE) +
  theme(legend.position="none")

comboplot2

G4vG4bigplot <- ggplot(G4vG4big,aes(as.factor(order),lfc)) + 
  geom_bar(stat="identity",position="stack",width=1,fill=c(rep("gray20",6),rep("gray80",179),rep("gray20",6),rep("gray20",6),rep("gray80",99),rep("gray20",6))) + 
  theme(panel.ontop = TRUE,plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill=NA, colour = "black"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="transparent")) +
  scale_y_continuous(limits = c(-7.8, 7.8),breaks=seq(-7,7,by=1),expand=c(0,0)) +
  facet_grid(~factor(comp,levels=c("G4DOXAvG4H2OA","G4DOXAvG4ATCA")),scales="free_x",space="free") +
  geom_hline(yintercept = 0,colour="gray60")

G4vG3bigplot <- ggplot(G4vG3big,aes(as.factor(order),lfc)) + 
  geom_bar(stat="identity",position="dodge",width=1,fill=c(rep("gray20",6),rep("gray80",258),rep("gray20",6),rep("gray20",6),rep("gray80",108),rep("gray20",6))) + 
  theme(panel.ontop = TRUE,plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill=NA, colour = "black"), panel.grid.minor=element_blank(), panel.grid.major = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="transparent")) +
  scale_y_continuous(limits = c(-7.8, 7.8),breaks=seq(-7,7,by=1),expand=c(0,0)) +
  facet_grid(~factor(comp,levels=c("G4DOXAvG3DOXA","G4ATCAvG3ATCA")),scales="free_x",space="free") +
  scale_x_discrete(labels=G4vG3big[,5]) +
  geom_hline(yintercept = 0,colour="gray60")

top_row <- plot_grid(comboplot1, G4vG4bigplot, ncol=2,axis="bt",align="h",rel_widths=c(1,0.57))
bottom_row <- plot_grid(comboplot2, G4vG3bigplot, ncol=2,axis="bt",align="h",rel_widths=c(0.89,0.8))
lfc_plots <- plot_grid(top_row,bottom_row,ncol=1)
lfc_plots

############## qPCR data ############## 

# Import qPCR data file #
qPCRstats <- read.csv("qPCRdata2.csv")
qPCRstats_L <- qPCRstats[qPCRstats$stage == "Larvae",]
qPCRstats_A <- qPCRstats[qPCRstats$stage == "Adults",]

### Statistical analysis ###
# we assume that qPCR data are normally distributed
# we can check with Shapiro-Wilk normality test on our control datasets (the H2O samples)
shapiro.test(qPCRstats_A$ddCt16S[c(1:6,19:24)])  # 16S adults, p-value = 0.9999
shapiro.test(qPCRstats_A$ddCtCOXI[c(1:6,19:24)])  # COXI adults, p-value = 0.5298
shapiro.test(qPCRstats_L$ddCt16S[c(1:3,10:12)]) # 16S larvae, p-value = 0.8666
shapiro.test(qPCRstats_L$ddCtCOXI[c(1:3,10:12)]) # COXI larvae, p-value = 0.2282

## 16S gene: ANOVA (global result) and pairwise-t tests for group comparisons
## BH correction applied to account for multiple comparisons

# First, run code that allows us to display degrees of freedom and t-value (test statistic) using a custom function

pairwise.t.test.with.t.and.df <- function (x, g, p.adjust.method = p.adjust.methods, pool.sd = !paired, 
                                           paired = FALSE, alternative = c("two.sided", "less", "greater"), 
                                           ...) 
{
  if (paired & pool.sd) 
    stop("pooling of SD is incompatible with paired tests")
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)
  p.adjust.method <- match.arg(p.adjust.method)
  alternative <- match.arg(alternative)
  if (pool.sd) {
    METHOD <- "t tests with pooled SD"
    xbar <- tapply(x, g, mean, na.rm = TRUE)
    s <- tapply(x, g, sd, na.rm = TRUE)
    n <- tapply(!is.na(x), g, sum)
    degf <- n - 1
    total.degf <- sum(degf)
    pooled.sd <- sqrt(sum(s^2 * degf)/total.degf)
    compare.levels <- function(i, j) {
      dif <- xbar[i] - xbar[j]
      se.dif <- pooled.sd * sqrt(1/n[i] + 1/n[j])
      t.val <- dif/se.dif
      if (alternative == "two.sided") 
        2 * pt(-abs(t.val), total.degf)
      else pt(t.val, total.degf, lower.tail = (alternative == 
                                                 "less"))
    }
    compare.levels.t <- function(i, j) {
      dif <- xbar[i] - xbar[j]
      se.dif <- pooled.sd * sqrt(1/n[i] + 1/n[j])
      t.val = dif/se.dif 
      t.val
    }       
  }
  else {
    METHOD <- if (paired) 
      "paired t tests"
    else "t tests with non-pooled SD"
    compare.levels <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$p.value
    }
    compare.levels.t <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$statistic
    }
    compare.levels.df <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$parameter
    }
  }
  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  TVAL <- pairwise.table.t(compare.levels.t, levels(g), p.adjust.method)
  if (pool.sd) 
    DF <- total.degf
  else
    DF <- pairwise.table.t(compare.levels.df, levels(g), p.adjust.method)           
  ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, 
              p.adjust.method = p.adjust.method, t.value = TVAL, dfs = DF)
  class(ans) <- "pairwise.htest"
  ans
}
pairwise.table.t <- function (compare.levels.t, level.names, p.adjust.method) 
{
  ix <- setNames(seq_along(level.names), level.names)
  pp <- outer(ix[-1L], ix[-length(ix)], function(ivec, jvec) sapply(seq_along(ivec), 
                                                                    function(k) {
                                                                      i <- ivec[k]
                                                                      j <- jvec[k]
                                                                      if (i > j)
                                                                        compare.levels.t(i, j)               
                                                                      else NA
                                                                    }))
  pp[lower.tri(pp, TRUE)] <- pp[lower.tri(pp, TRUE)]
  pp
}

res.aov.16S <- aov(ddCt16S ~ group, data = qPCRstats)
summary(res.aov.16S)
comps_16S <- pairwise.t.test(qPCRstats$ddCt16S, qPCRstats$group,p.adjust.method="BH")
comps_16S$p.value

# Adults, 16S gene
res.aov.16S_A <- aov(ddCt16S ~ group, data = qPCRstats_A)
summary(res.aov.16S_A)
comps_16S_A <- pairwise.t.test(qPCRstats_A$ddCt16S, qPCRstats_A$group,p.adjust.method="BH")
comps_16S_A$p.value

# Larvae, 16S gene
res.aov.16S_L <- aov(ddCt16S ~ group, data = qPCRstats_L)
summary(res.aov.16S_L)
comps_16S_L <- pairwise.t.test(qPCRstats_L$ddCt16S, qPCRstats_L$group,p.adjust.method="BH")
comps_16S_L$p.value

# Adults, COXI gene
res.aov.COXI_A <- aov(ddCtCOXI ~ group, data = qPCRstats_A)
summary(res.aov.COXI_A)
comps_COXI_A <- pairwise.t.test(qPCRstats_A$ddCtCOXI, qPCRstats_A$group,p.adjust.method="BH")
comps_COXI_A$p.value

# Larvae, COXI gene
res.aov.COXI_L <- aov(ddCtCOXI ~ group, data = qPCRstats_L)
summary(res.aov.COXI_L)
comps_COXI_L <- pairwise.t.test(qPCRstats_L$ddCtCOXI, qPCRstats_L$group,p.adjust.method="BH")
comps_COXI_L$p.value

result_16S_A <- pairwise.t.test.with.t.and.df(qPCRstats_A$ddCt16S, qPCRstats_A$group, p.adjust.method="BH")
result_16S_A[[5]] # t-values 
result_16S_A[[6]] # df

result_16S_L <- pairwise.t.test.with.t.and.df(qPCRstats_L$ddCt16S, qPCRstats_L$group, p.adjust.method="BH")
result_16S_L[[5]] # t-values 
result_16S_L[[6]] # df

result_COXI_A <- pairwise.t.test.with.t.and.df(qPCRstats_A$ddCtCOXI, qPCRstats_A$group, p.adjust.method="BH")
result_COXI_A[[5]] # t-values 
result_COXI_A[[6]] # df

result_COXI_L <- pairwise.t.test.with.t.and.df(qPCRstats_L$ddCtCOXI, qPCRstats_L$group, p.adjust.method="BH")
result_COXI_L[[5]] # t-values 
result_COXI_L[[6]] # df

## Export output from the pairwise t-tests

# p-values (adjusted)
write.csv(comps_16S_A$p.value,"qPCR_16S_adults_pairwise_comps.csv") # 16S adults
write.csv(comps_16S_L$p.value,"qPCR_16S_larvae_pairwise_comps.csv") # 16S larvae
write.csv(comps_COXI_A$p.value,"qPCR_COXI_adults_pairwise_comps.csv") # COXI adults
write.csv(comps_COXI_L$p.value,"qPCR_COXI_larvae_pairwise_comps.csv") # COXI larvae

# t-values
write.csv(result_16S_A[[5]],"qPCR_16SA_tvals.csv")
write.csv(result_16S_L[[5]],"qPCR_16SL_tvals.csv")
write.csv(result_COXI_A[[5]],"qPCR_COXIA_tvals.csv")
write.csv(result_COXI_L[[5]],"qPCR_COXIL_tvals.csv")

# prepare p-value data frame for plotting significance groups
pvals_16S_A_df <- as.data.frame(comps_16S_A$p.value)
pvals_16S_L_df <- as.data.frame(comps_16S_L$p.value)
pvals_COXI_A_df <- as.data.frame(comps_COXI_A$p.value)
pvals_COXI_L_df <- as.data.frame(comps_COXI_L$p.value)

pvals_16S_A_df <- add_row(pvals_16S_A_df,.before = 1)
pvals_16S_L_df <- add_row(pvals_16S_L_df,.before = 1)
pvals_COXI_A_df <- add_row(pvals_COXI_A_df,.before = 1)
pvals_COXI_L_df <- add_row(pvals_COXI_L_df,.before = 1)

rownames(pvals_16S_A_df)[1] <- "G3 Adults ATC"
rownames(pvals_COXI_A_df)[1] <- "G3 Adults ATC"
rownames(pvals_16S_L_df)[1] <- "G3 Larvae ATC"
rownames(pvals_COXI_L_df)[1] <- "G3 Larvae ATC"

pvals_16S_A_df <- cbind(pvals_16S_A_df,empty_column=NA)
pvals_COXI_A_df <- cbind(pvals_COXI_A_df,empty_column=NA)
pvals_16S_L_df <- cbind(pvals_16S_L_df,empty_column=NA)
pvals_COXI_L_df <- cbind(pvals_COXI_L_df,empty_column=NA)

colnames(pvals_16S_A_df)[6] <- "G4 Adults H2O"
colnames(pvals_COXI_A_df)[6] <- "G4 Adults H2O"
colnames(pvals_16S_L_df)[6] <- "G4 Larvae H2O"
colnames(pvals_COXI_L_df)[6] <- "G4 Larvae H2O"

# make the matrices of p-values symmetrical
pvals_16S_A_df.sym <- Matrix::forceSymmetric(as.matrix(pvals_16S_A_df),uplo="L")
pvals_16S_L_df.sym <- Matrix::forceSymmetric(as.matrix(pvals_16S_L_df),uplo="L")
pvals_COXI_A_df.sym <- Matrix::forceSymmetric(as.matrix(pvals_COXI_A_df),uplo="L")
pvals_COXI_L_df.sym <- Matrix::forceSymmetric(as.matrix(pvals_COXI_L_df),uplo="L")

# add new columns to the qPCRstats data
qPCRstats[,10:11] <- c(NA,NA)

# get letters denoting significantly different groups (p<0.05)

letters1 <- multcompLetters(pvals_16S_A_df.sym)$Letters
letters2 <- multcompLetters(pvals_16S_L_df.sym)$Letters
for (x in 1:54) {qPCRstats[x, 10] <- c(letters1, letters2)[[qPCRstats$group[x]]]}

letters3 <- multcompLetters(pvals_COXI_A_df.sym)$Letters
letters4 <- multcompLetters(pvals_COXI_L_df.sym)$Letters
for (x in 1:54) {qPCRstats[x, 11] <- c(letters3, letters4)[[qPCRstats$group[x]]]}

# make new columns in qPCRstats with max values -- this will help us plot the letters above data points on the graph
for (x in 1:54) {qPCRstats[x, 12] <- max(qPCRstats[qPCRstats$group == qPCRstats$group[x],8])}
for (x in 1:54) {qPCRstats[x, 13] <- max(qPCRstats[qPCRstats$group == qPCRstats$group[x],9])}

# rename columns for easier plotting
colnames(qPCRstats)[10:13] <- c("letters_16S","letters_COXI","maxes_16S","maxes_COXI") 

### qPCR plots ###

group_colors <- c("G3 Larvae DOX" = "#CD96CD","G4 Larvae DOX" = "#8B668B","G3 Adults DOX" = "#CD96CD","G4 Adults DOX" = "#8B668B","G3 Adults ATC" = "#35C47D","G4 Adults ATC" = "#199151","G3 Larvae ATC" = "#35C47D","G4 Larvae ATC" = "#199151","G3 Larvae H2O" = "#5CACEE","G4 Larvae H2O" = "#1874CD","G3 Adults H2O" = "#5CACEE","G4 Adults H2O" = "#1874CD")
# Plot the COXI data
final_plot_COXI <- ggplot(qPCRstats, aes(forcats::fct_relevel(G3diet, "DOX", "ATC", "H2O"), -1*ddCtCOXI, group=generation, color=group)) + 
  scale_color_manual(values = group_colors) + 
  geom_point(size=2.5,shape=16,alpha=1,position = position_jitterdodge(jitter.width = 0, dodge.width = .7)) + theme(strip.background = element_rect(colour = "black"),panel.border = element_rect(linetype = "solid", fill = NA), panel.background = element_rect(fill = "white", colour = "grey30"), panel.grid.minor=element_blank(), axis.title.x = element_text(size=15),  axis.title.y = element_text(size=15), axis.text.y = element_text(size=18), axis.text.x = element_text(size=12, angle = 0, vjust = 0, hjust=0.5), strip.text.x= element_text(size=20), legend.text = element_text(size=10)) + 
  scale_y_continuous(limits = c(-2, 2)) + 
  stat_summary(fun.data=mean_se, color=c(rep(c("#CD96CD","#35C47D","#5CACEE","#8B668B","#199151","#1874CD"),2)),position=position_dodge(width=0.7), size=0.6, width=0.4, geom="errorbar") + stat_summary(fun.data=mean_se, color=c(rep(c("#CD96CD","#35C47D","#5CACEE","#8B668B","#199151","#1874CD"),2)),geom="point",position=position_dodge(width=0.7), size=2.5, shape=15) + 
  xlab("G1-G3 Diet") + ylab("-\u0394\u0394Ct: COXI") + 
  geom_hline(yintercept=0,alpha=0.3,linetype="dashed") + 
  geom_text(aes(label = str_trim(letters_COXI), y = rep(1.8,54)), vjust = -0.5, position=position_dodge(width=0.7),size=5,color="black") + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  facet_wrap2(~forcats::fct_relevel(stage, "Larvae", "Adults"), strip = strip)

final_plot_COXI

# Plot the 16S data
final_plot_16S <- ggplot(qPCRstats, aes(forcats::fct_relevel(G3diet, "DOX", "ATC", "H2O"), -1*ddCt16S, group=generation, color=group)) + 
  scale_color_manual(values = group_colors) +
  geom_point(size=2.5,shape=16,alpha=1,position = position_jitterdodge(jitter.width = 0, dodge.width = .7)) + 
  theme(strip.background = element_blank(),panel.border = element_rect(linetype = "solid", fill = NA), panel.background = element_rect(fill = "white", colour = "grey30"), panel.grid.minor=element_blank(), axis.title.x = element_text(size=15),  axis.title.y = element_text(size=15), axis.text.y = element_text(size=18), axis.text.x = element_text(size=12, angle = 0, vjust = 0, hjust=0.5), strip.text.x= element_text(size=0), legend.text = element_text(size=10)) + 
  scale_y_continuous(limits = c(-6, 6)) + 
  stat_summary(fun.data=mean_se, color=c(rep(c("#CD96CD","#35C47D","#5CACEE","#8B668B","#199151","#1874CD"),2)),position=position_dodge(width=0.7), size=0.6, width=0.4, geom="errorbar") + 
  stat_summary(fun.data=mean_se, color=c(rep(c("#CD96CD","#35C47D","#5CACEE","#8B668B","#199151","#1874CD"),2)),geom="point",position=position_dodge(width=0.7), size=2.5, shape=15) + 
  xlab("G1-G3 Diet") + 
  ylab("-\u0394\u0394Ct: 16S") + 
  geom_hline(yintercept=0,alpha=0.3,linetype="dashed") +
  theme(legend.position = "none", axis.title.x = element_blank(),axis.text.x = element_blank()) + 
  facet_wrap2(~forcats::fct_relevel(stage, "Larvae", "Adults")) + 
  geom_text(aes(label = str_trim(letters_16S), y = rep(5.6,54)), vjust = -0.5, position=position_dodge(width=0.7),size=5,color="black")
final_plot_16S 

mean(qPCRstats$ddCt16S[c(49:54)]) # G4 DOX A
mean(qPCRstats$ddCt16S[c(22:27)]) # G3 DOX A
mean(qPCRstats$ddCt16S[c(43:48)]) # G4 ATC A
mean(qPCRstats$ddCt16S[c(16:21)]) # G3 ATC A

# Plot both graphs in vertical stack
qPCR_plots_vertical <- plot_grid(final_plot_COXI, final_plot_16S, ncol=1, align="v")
qPCR_plots_vertical
# export portrait, 11x8"
