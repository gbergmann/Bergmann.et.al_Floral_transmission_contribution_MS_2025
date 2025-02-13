# Bioinformatics workflow to process field survey (FS) data
  # This script has three components: denoising with DADA2, compiling data into a phyloseq object, and filtering out contaminants and low-abundance ASVs.

#### loading packages ####
library(tidyverse)
library(magrittr)
library(dada2)
library(ShortRead)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(stringr)
library(readxl)
library(tibble)
library(foreach)
library(DECIPHER)
library(vegan)
library(ggthemes)
library(decontam)

#### Denoising with DADA2 ####
#path for output
out.path <- "data_output_seq/Field_survey/dada"

# Identify paths to trimmed fwd and rev reads
in.path <- file.path("data_raw/Field.survey_NovaSeq/")
path.fwd <- in.path %>% list.files(pattern="_R1.fastq.gz",full.names = T) %>% sort 
path.rev <- in.path %>% list.files(pattern="_R2.fastq.gz",full.names = T) %>% sort 

# make file names / paths for trimmed and filtered files
filt.path <- file.path(out.path, "filt")
filt.fwd <- path.fwd %>% gsub(in.path,filt.path,.)
filt.rev <- path.rev %>% gsub(in.path,filt.path,.)

# preview read quality of 12 randomly selected samples before trimming
qual.samps <- sample(1:length(path.fwd),6)
plotQualityProfile(path.fwd[qual.samps]) 
plotQualityProfile(path.rev[qual.samps])

# Filter and trim, dereplicate identical sequences
trim <- filterAndTrim(path.fwd, filt.fwd,
                      path.rev, filt.rev,
                      maxN=0, maxEE=c(2,2), truncQ=2,
                      multithread=TRUE)

# Update list of trimmed file paths to exclude samples with no reads passing filters
filt.fwd <- list.files(filt.path, pattern = "_R1.fastq.gz",full.names = T)
filt.rev <- list.files(filt.path, pattern = "_R2.fastq.gz",full.names = T)

# Check quality of trimmed and filtered reads
qual.samps <- sample(1:length(filt.fwd),12)
plotQualityProfile(filt.fwd[qual.samps])
plotQualityProfile(filt.rev[qual.samps]) 

# Dereplicate identical sequences
derep.fwd <- derepFastq(filt.fwd, verbose=F)
derep.rev <- derepFastq(filt.rev, verbose=F) 

# Trim names of derep objects to just sample names
names(derep.fwd) %<>% gsub("_R1.fastq.gz","",.)
names(derep.rev) %<>% gsub("_R2.fastq.gz","",.)

# Learn errors, denoise, merge sequences into table
err.fwd <- learnErrors(filt.fwd, multithread=TRUE, nbases = 3e08) 
plotErrors(err.fwd, nominalQ=TRUE)
err.rev <- learnErrors(filt.rev, multithread=TRUE, nbases = 3e08)
plotErrors(err.rev, nominalQ=TRUE)

# Denoise
dada.fwd <- dada(derep.fwd, err=err.fwd, multithread=TRUE)
dada.rev <- dada(derep.rev, err=err.rev, multithread=TRUE)

# Merge reads
merged <- mergePairs(dada.fwd, derep.fwd, dada.rev, derep.rev, trimOverhang = T)

# Make sequence table
seqtab <- merged %>% makeSequenceTable

# Remove chimeras, make summary report
seqtab.nonChimeras <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

# Make summary report of what happened
getN <- function(x) sum(getUniques(x))
trim.summary <- trim %>% data.frame %>% rownames_to_column("Sample") 
trim.summary$Sample %<>% str_sub(., 1, 12)
track <- cbind(sapply(dada.fwd, getN), 
               sapply(dada.rev, getN), 
               sapply(merged, getN),
               rowSums(seqtab.nonChimeras)) %>% 
  data.frame %>% rownames_to_column("Sample") 
track %<>% left_join(trim.summary,.)
colnames(track) <- c("sample", "input", "filtered", "denoised.fwd", "denoised.rev", "merged", "ChimeraFiltered")
file.path(out.path,"dadaSummary.csv") %>% write.csv(track,.,row.names = F)

# Save output 
file.path(out.path,"seqTab.rds") %>% saveRDS(seqtab.nonChimeras,.)

#### Assigning taxonomy, compiling into phyloseq object ####
# read in sample metadata (row.names refers to column one, where the sample IDs are)
meta <- read.csv("data_raw/Field.survey_NovaSeq/clean.bact.rel.meta.csv",as.is=T,row.names=1)

# read in denoised sequence table
seqTab <- readRDS("data_output_seq/Field_survey/dada/seqTab.rds") # alternatively, use seqTab <- seqtab

# extract sequences from OTU table 
seqs <- getSequences(seqTab) %>% DNAStringSet

# rename sequences
ASV.names <- paste0("ASV.",1:ncol(seqTab))
colnames(seqTab) <- ASV.names
names(seqs) <- ASV.names

# Create data frame for summary of processing
compile.summary <- data.frame(SampID=rownames(seqTab),
                              denoised=rowSums(seqTab))

# create scratch directory
dir.create("data_output_seq/Field_survey/scratch")
# write temp file
writeXStringSet(seqs,file="data_output_seq/Field_survey/scratch/tmp.fasta",width=1000)

# load a version of the SILVA database 
load("data_raw/Field.survey_NovaSeq/SILVA_SSU_r138_2019.RData") # this is not provided, download the latest SILVA release from DECIPHER here: https://www2.decipher.codes/Downloads.html

# remove sequences that are shorter than 50 bp
seqs.16S <- readDNAStringSet("data_output_seq/Field_surveyscratch/tmp.fasta") %>% .[.@ranges@width > 50]
seqs.16S %>% writeXStringSet(., file="data_output_seq/Field_survey/scratch/16S.tmp.fasta", width=1000)

# remove ASVs from seqTab that are too short
seqTab %<>% .[,names(seqs.16S)]

# collapse any identical sequences after isolating ITS2 
colnames(seqTab) <- seqs.16S %>% as.character %>% unname
seqTab %<>% collapseNoMismatch

# assign taxonomy with IdTaxa from the DECIPHER package (recommended by DADA2 developer, faster than assignTaxonomy) - A Murali et al. (2018) "IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences." Microbiome, doi:10.1186/s40168-018-0521-5.
ids <- IdTaxa(seqs.16S, trainingSet, strand = "top", processors = NULL, verbose = F)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks 
rownames(taxid) <- getSequences(seqTab)

# set final ASV names
final.names <- paste0("ASV.",1:ncol(seqTab))

final.seqs <- getSequences(seqTab) %>% DNAStringSet()
names(final.seqs) <- final.names
colnames(seqTab) <- final.names
rownames(taxid) <- final.names

# ake phyloseq object
phy <- phyloseq(otu_table(seqTab,taxa_are_rows = F), 
                sample_data(meta),
                tax_table(taxid),
                refseq(final.seqs))

# remove non-bacterial taxa from phyloseq object
min.reads <- 0
keep.samples <- sample_sums(phy) > min.reads
phy %<>% subset_taxa(., domain=="Bacteria") %>% prune_samples(keep.samples,.) %>% filter_taxa(function(x) {sum(x) > 0}, TRUE)

#Save phyloseq object
dir.create("data_output_seq/Field_survey/compiled")
saveRDS(phy,"data_output_seq/Field_survey/compiled/phy.rds")

#Save components for possible manual inspection
otu_table(phy) %>% write.csv("data_output_seq/Field_survey/compiled/ASV.table.csv")
tax_table(phy) %>% write.csv("data_output_seq/Field_survey/compiled/taxonomy.table.csv")

#### Filtering out putative contaminants, low-depth samples and singleton ASVs ####
# removing host DNA (identified as chloroplast or mitochondria)
min.reads <- 0
keep.samples <- sample_sums(phy) > min.reads
bact.phy <- subset_taxa(phy, order!="Chloroplast") %>% 
  subset_taxa(., family!="Mitochondria") %>%
  prune_samples(keep.samples,.) %>% 
  filter_taxa(function(x) {sum(x) > 0}, TRUE)

phy.taxa <- tax_table(bact.phy) %>% as.data.frame

# prepare sample and control data for decontam
df <- as.data.frame(sample_data(bact.phy))
df$LibrarySize <- sample_sums(bact.phy)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=control)) + geom_point()

# identify negative controls
sample_data(bact.phy)$is.neg <- sample_data(bact.phy)$control == "yes"

# run decontam with threshold 0.2
contamdf.prev2 <- isContaminant(bact.phy, method="prevalence", neg="is.neg", threshold=0.2) 
table(contamdf.prev2$contaminant) 
  # 0.2: FALSE  TRUE 
  # 0.2:  4012   9   
head(which(contamdf.prev2$contaminant))

# show contaminant ASVs with associated taxonomy
contamdf.prev.list2 <- contamdf.prev2 %>% 
  filter(contaminant == TRUE)
contamdf.prev.list.taxo2 <-  merge(contamdf.prev.list2, 
                                   as.data.frame(tax_table(bact.phy)), 
                                   by="row.names")

contamdf.prev.list.taxo2 %>% write.csv("output/compiled/decontam_results_threshold02.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(bact.phy, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$control == "yes", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$control == "no", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa2 <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                     contaminant=contamdf.prev2$contaminant)
ggplot(df.pa2, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

contamdf.prev2[contamdf.prev2$contaminant == TRUE, ]  #9 true contaminants
# Prune TRUE contaminants in phyloseq objects
phy_v2 <- prune_taxa(!contamdf.prev2$contaminant, bact.phy)
phy_v2 #4012 taxa and 82 samples

phy.asv <- otu_table(phy_v2) %>% as.data.frame()
phy.asv %>% write.csv("data_output_seq/Field_survey/compiled/post.decontam.asv.csv")

# Look at sequencing depth
plotSeqDepth <- function(phy,cutoff){
  pctPass <- phy_v2 %>% sample_data %>% data.frame %>% 
    mutate(SeqDepth=sample_sums(phy_v2),
           pass=SeqDepth>cutoff) %>% 
    summarise(nSamps=n(),pctPass=round(100*sum(pass)/n(),1))
  phy_v2 %>% sample_data %>% data.frame %>%  
    mutate(SeqDepth=sample_sums(phy_v2)) %>% 
    .[order(.$SeqDepth),] %>%
    mutate(Index=factor(seq(nrow(.)))) %>%
    ggplot(aes(x=Index, y=SeqDepth,color=SeqDepth<cutoff)) + 
    geom_point(position=position_jitter(width=50),alpha=0.5,shape=19) +
    geom_text(aes(x=-Inf,y=Inf,label=paste0(nSamps," samples")),
              data=pctPass,inherit.aes = F, hjust=-0.1, vjust=1.5) +
    geom_text(aes(x=-Inf,y=Inf,label=paste0(pctPass, "% > ",cutoff, " seqs")),
              data=pctPass,inherit.aes = F, hjust=-0.1, vjust=3) +
    scale_color_calc(name=paste0("Sequencing depth > ",cutoff,"\n")) +
    theme_classic() + 
    theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank()) 
}
plotSeqDepth(phy_v2,500)
minDepth <- 500 # this minDepth is based off the cutoffs tested with the plot

# Remove samples below sequencing depth cutoff
(phy_v2 %<>% prune_samples(sample_sums(.)>minDepth,.))

# Remove low abundance ASVs by prevalence - Only keep ASVs present in at least 2.5% of samples
phy_v3 <- filter_taxa(phy_v2, function (x) {sum(x > 0) > 1}, prune=TRUE)
bact.rich.rel.phy <- transform_sample_counts(phy_v3, function(x){(x/sum(x))*100})
bact.rel.phy <- transform_sample_counts(phy_v2, function(x){(x/sum(x))*100})

# save cleaned phyloseq objects (optional here)
saveRDS(phy_v2, "data_output_seq/Field_survey/clean.bact.raw.phy.rds")
saveRDS(bact.rel.phy, "data_output_seq/Field_survey/clean.bact.rel.phy.rds")
saveRDS(phy_v3, "data_output_seq/Field_survey/clean.bact.rich.raw.phy.rds")
saveRDS(bact.rich.rel.phy, "data_output_seq/Field_survey/clean.bact.rich.rel.phy.rds")
