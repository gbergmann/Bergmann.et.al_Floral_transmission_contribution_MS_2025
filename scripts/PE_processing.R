# Bioinformatics workflow to process pollinator experiment (PE) dataset

#### loading required packages ####
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
library(ggthemes)
library(decontam)

#### Part 1: cutadapt (done in command line) ####

  #for i in `cat group`; do cutadapt --discard-untrimmed -o $i.gyrB.R1.fq -p $i.gyrB.R2.fq -g   #MGNCCNGSNATGTAYATHGG -G ACNCCRTGNARDCCDCCNGA -e 0.1  -O 20 $i*L001_R1_001.fastq.gz $i*_L001_R2_001.fastq.gz #; done
  #
  #for i in `cat group`; do cutadapt --discard-untrimmed -o $i.ITS1.R1.fq -p $i.ITS1.R2.fq -g   #CTTGGTCATTTAGAGGAAGTAA  -G GCTGCGTTCTTCATCGATGC -e 0 -O 20 $i*"L001_R1_001.fastq.gz" $i*"L001_R2_001.fastq  #.gz"; done

#### Part 2: denoising ####
path <- "data_raw/Pollinator.Exp_MiSeq/" 
list.files(path)
fnFs <- sort(list.files(path, pattern="gyrB.R1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="gyrB.R2.fq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:8]) #quality drop around 130 (dimers ?) 
plotQualityProfile(fnRs[1:8]) #quality drop around 100 (dimers ?)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(160,130),
                     maxN=0, maxEE=c(2,2),  rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # Need to change this for increasing the number of sequences
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) #without pooling or pseudo-pooling (no need to detect rare ASV)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
#Since gyrB is a protein-coding genes only triplets should be conserved (244-247-250-253-256-259-262-265-268)
seqtab244 <- seqtab[,nchar(colnames(seqtab)) %in% 244]
seqtab247 <- seqtab[,nchar(colnames(seqtab)) %in% 247]
seqtab250 <- seqtab[,nchar(colnames(seqtab)) %in% 250]
seqtab253 <- seqtab[,nchar(colnames(seqtab)) %in% 253]
seqtab256 <- seqtab[,nchar(colnames(seqtab)) %in% 256]
seqtab259 <- seqtab[,nchar(colnames(seqtab)) %in% 259]
seqtab262 <- seqtab[,nchar(colnames(seqtab)) %in% 262]
seqtab265 <- seqtab[,nchar(colnames(seqtab)) %in% 265]
seqtab268 <- seqtab[,nchar(colnames(seqtab)) %in% 268]
#Merge all files
seq.final <- cbind(seqtab244, seqtab247, seqtab250, seqtab253, seqtab256, seqtab259, seqtab262, seqtab265, seqtab268) 
dim(seq.final)
sum(seq.final)/sum(seqtab)
#Detect/Remove chimera
seqtab.nochim <- removeBimeraDenovo(seq.final, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seq.final)
#Summary
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
saveRDS(seqtab.nochim, "Gillian_ASV_gyrB.rds")

allgyrB<-seqtab.nochim[, order(colSums(-seqtab.nochim))]
taxa.g <- assignTaxonomy(allgyrB, "C:/Documents and Settings/mbarret/Mes documents/DB/train_set_gyrB_v5.fa.gz", multithread=TRUE)
taxa.print.g <- taxa.g # Removing sequence rownames for display only
rownames(taxa.print.g) <- NULL
head(taxa.print.g)
saveRDS(taxa.g, "Gillian_gyrB_taxo.rds")

#### Part 3: Compile data into phyloseq object #### 
asvgyrB <- readRDS("Gillian_ASV_gyrB.rds")
taxgyrB <- readRDS("Gillian_gyrB_taxo.rds")
design <- read.csv("design.csv", sep = ";",  check.names=FALSE, row.names=1)
psgyrB_0 <- phyloseq(tax_table(taxgyrB), sample_data(design),
                     otu_table(asvgyrB, taxa_are_rows = FALSE)) 
psgyrB_0# 3359 taxa and 132 samples
psgyrB_1 <- subset_taxa(psgyrB_0, !is.na(Kingdom) & !Kingdom %in% c("parE") & !is.na(Phylum)) psgyrB_1 # 2117 taxa and 132 samples 

# Renaming ASVs
dna.gyrB <- Biostrings::DNAStringSet(taxa_names(psgyrB_1))
names(dna.gyrB) <- taxa_names(psgyrB_1)
psgyrB_2 <- merge_phyloseq(psgyrB_1, dna.gyrB)
taxa_names(psgyrB_2) <- paste0("ASV", seq(ntaxa(psgyrB_2)))
gyrB.phy <- psgyrB_2

#### Part 4: filtering out contaminants and low-depth samples ####
df <- as.data.frame(sample_data(gyrB.phy))
df$LibrarySize <- sample_sums(gyrB.phy)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

#identify negative controls
sample_data(gyrB.phy)$is.neg <- sample_data(gyrB.phy)$Sample_or_Control == "Control"

#run decontam with threshold 0.5
contamdf.prev <- isContaminant(gyrB.phy, method="prevalence", neg="is.neg", threshold=0.2) 
table(contamdf.prev$contaminant) 
# show contaminant ASVs with associated taxonomy
contamdf.prev.list <- contamdf.prev %>% 
  filter(contaminant == TRUE)
contamdf.prev.list.taxo <-  merge(contamdf.prev.list, 
                                   as.data.frame(tax_table(gyrB.phy)), 
                                   by="row.names")
ps.pa <- transform_sample_counts(gyrB.phy, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$control == "yes", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$control == "no", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

contamdf.prev[contamdf.prev$contaminant == TRUE, ]  #9 true contaminants
# Prune TRUE contaminants in phyloseq objects
gyrB.phy.v2 <- prune_taxa(!contamdf.prev$contaminant, gyrB.phy)
gyrB.phy.v2 
gyrB.phy.asv <- otu_table(gyrB.phy.v2) %>% as.data.frame()
gyrB.phy.asv %>% write.csv("data_output_seq/Pollinator_exp/post.decontam.asv.csv")

plotSeqDepth <- function(phy,cutoff){
  pctPass <- gyrB.phy.v2 %>% sample_data %>% data.frame %>% 
    mutate(SeqDepth=sample_sums(phy),
           pass=SeqDepth>cutoff) %>% 
    #group_by(DataSet) %>%
    summarise(nSamps=n(),pctPass=round(100*sum(pass)/n(),1))
  gyrB.phy.v2 %>% sample_data %>% data.frame %>%  
    mutate(SeqDepth=sample_sums(gyrB.phy.v2)) %>% 
    .[order(.$SeqDepth),] %>%
    mutate(Index=factor(seq(nrow(.)))) %>%
    ggplot(aes(x=Index, y=SeqDepth,color=SeqDepth<cutoff)) + 
    geom_point(position=position_jitter(width=50),alpha=0.5,shape=19) +
    geom_text(aes(x=-Inf,y=Inf,label=paste0(nSamps," samples")),data=pctPass,inherit.aes = F, hjust=-0.1, vjust=1.5) +
    geom_text(aes(x=-Inf,y=Inf,label=paste0(pctPass, "% > ",cutoff, " seqs")),data=pctPass,inherit.aes = F, hjust=-0.1, vjust=3) +
    scale_color_calc(name=paste0("Sequencing depth > ",cutoff,"\n")) +
    #facet_wrap(~DataSet,scales="free") +
    theme_classic() + 
    theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank()) 
}
plotSeqDepth(gyrB.phy.v2,2000)
minDepth <- 2000 # this minDepth is based off the cutoffs tested with the plot

# Remove samples below sequencing depth cutoff
(gyrB.phy.v2 %<>% prune_samples(sample_sums(.)>minDepth,.))
gyrB.rel.phy <- transform_sample_counts(gyrB.phy, function(x){(x/sum(x))*100})

saveRDS(gyrB.phy, "data_output_seq/Pollinator_exp/PE.gyrB.raw.phy.rds")
saveRDS(gyrB.rel.phy, "data_output_seq/Pollinator_exp/PE.gyrB.rel.phy.rds")
