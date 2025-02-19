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
