# Analysis of field survey (FS) dataset

#### loading packages ####
library(tidyverse)
library(phyloseq)
library(vegan)
library(betapart)
library(viridis)
library(plotly)
library(pheatmap)
library(microbiomeutilities) 
library(pairwiseAdonis)
library(fantaxtic)
library(indicspecies)
library(reshape2)
library(RColorBrewer)
library(microbiome)
library(picante)
library(ggpubr) 
library(rcartocolor)
library(hrbrthemes)
library(cowplot)
library(khroma)
library(ANCOMBC)
library(betapart)
library(phangorn) 
library(DECIPHER)

#### loading data ####
bact.rich.rel.phy <- readRDS("data_output_seq/Field_survey/clean.bact.rich.rel.phy.rds")
bact.rel.phy <- readRDS("data_output_seq/Field_survey/clean.bact.rel.phy.rds")

#### alpha diversity (observed richness and Faith's PD) ####

#### NMDS and PerMANOVA ####

#### overlap between seeds and stigmas across all samples ####

#### visualizing most abundant shared + unique genera ####

#### Differential abundance analysis with ANCOM-BC2 ####