# Analysis of pollinator experiment (PE) dataset

#### loading packages and data ####
library(tidyverse)
library(phyloseq)
library(vegan)
library(viridis)
library(microbiomeutilities) 
library(pairwiseAdonis)
library(fantaxtic)
library(reshape2)
library(cowplot)
library(microbiome)
library(rstatix)
library(ggpubr)
library(picante)
library(khroma)
library(ANCOMBC)
library(UpSetR)
library(ComplexUpset)
library(paletteer)
library(cowplot)
library(DECIPHER)
library(phangorn)

# and data
gyrB.rel.phy <- readRDS("data_output_seq/Pollinator_exp/PE.gyrB.rel.phy.rds") 

#### UpSet plots bee vs. seed ASVs in each treatment (Fig 4A) ####
hand.bee.phy <- subset_samples(gyrB.rel.phy, treatment == "hand"|tissue_type == "bee_pooled")
hand.bee.phy <- subset_samples(hand.bee.phy, tissue_type != "stigma")
hand.bee.phy <- prune_taxa(taxa_sums(hand.bee.phy)>0, hand.bee.phy)

hand.var1 = as.character(get_variable(hand.bee.phy, "tissue_type"))
sample_data(hand.bee.phy)$tissue <- mapply(paste0, hand.var1, collapse = "_")
merge_samples(hand.bee.phy, "tissue")
merged_hand.bee_all <- merge_samples(hand.bee.phy, "tissue")
#Fix continuous variable to discrete value
sample_data(hand.bee.phy)$tissue <- factor(sample_names(hand.bee.phy),
                                           levels = c("bee_pooled", "seed_pooled"))
#Extract count datasets from phyloseq
count_hand.bee_all <- t(otu_table(merged_hand.bee_all))
#Transform in binary matrix
count_hand.bee_all[count_hand.bee_all> 0] <- 1
#Convert to dataframe
df_hand.bee_all <- as.data.frame(count_hand.bee_all)
#Read abundance per ASV
hand.bee.rel.abu <- data.frame(taxa_sums(hand.bee.phy))
hand.bee.rel.abu$log <- log10(hand.bee.rel.abu$taxa_sums.hand.bee.phy.)
#Merge ASV dataframe and ASV abundance
df_hand.bee_all_upset <-merge.data.frame(df_hand.bee_all, hand.bee.rel.abu, by="row.names",
                                         all.x=TRUE)
#Append taxonomy to upset dataframe
df_hand.bee_taxa <- tax_table(hand.bee.phy) %>% data.frame() %>% rownames_to_column()
names(df_hand.bee_taxa)[names(df_hand.bee_taxa)=="rowname"] <- "Row.names"
df_hand.bee_all_upset2 <- full_join(df_hand.bee_all_upset, df_hand.bee_taxa, by="Row.names")

# create the plot
df_hand.bee_all_upset2 <- df_hand.bee_all_upset2[, c(1,2,3,4,5,6,7,8,9,10,11,12)]    
tissues <- colnames(df_hand.bee_all_upset2)[2:3]
df_hand.bee_all_upset2[tissues] = df_hand.bee_all_upset2[tissues] == 1
t(head(df_hand.bee_all_upset2[tissues], 3))
df_hand.bee_all_upset2[is.na(df_hand.bee_all_upset2)] <- "unclassified"

Upset_hand_bee <- ComplexUpset::upset(
  df_hand.bee_all_upset2,
  tissues,
  name = 'Tissue type',
  width_ratio=0.1
)
Upset_hand_bee2 <- Upset_hand_bee + theme(axis.title = element_text(size=20), axis.text = element_text(size=15))

# doing the same thing for bees and seeds in the insect-pollinated treatment
ins.bee.phy <- subset_samples(gyrB.rel.phy, treatment == "insect"|tissue_type == "bee_pooled")
ins.bee.phy <- subset_samples(ins.bee.phy, tissue_type != "stigma")
ins.bee.phy <- prune_taxa(taxa_sums(ins.bee.phy)>0, ins.bee.phy)

ins.var1 = as.character(get_variable(ins.bee.phy, "tissue_type"))
sample_data(ins.bee.phy)$tissue <- mapply(paste0, ins.var1, collapse = "_")
merge_samples(ins.bee.phy, "tissue")
merged_ins.bee_all <- merge_samples(ins.bee.phy, "tissue")
#Fix continuous variable to discrete value
sample_data(ins.bee.phy)$tissue <- factor(sample_names(ins.bee.phy),
                                          levels = c("bee_pooled", "seed_pooled"))
#Extract count datasets from phyloseq
count_ins.bee_all <- t(otu_table(merged_ins.bee_all))
#Transform in binary matrix
count_ins.bee_all[count_ins.bee_all> 0] <- 1
#Convert to dataframe
df_ins.bee_all <- as.data.frame(count_ins.bee_all)
#Read abundance per ASV
ins.bee.rel.abu <- data.frame(taxa_sums(ins.bee.phy))
ins.bee.rel.abu$log <- log10(ins.bee.rel.abu$taxa_sums.ins.bee.phy.)
#Merge ASV dataframe and ASV abundance
df_ins.bee_all_upset <-merge.data.frame(df_ins.bee_all, ins.bee.rel.abu, by="row.names",
                                        all.x=TRUE)
#Append taxonomy to upset dataframe
df_ins.bee_taxa <- tax_table(ins.bee.phy) %>% data.frame() %>% rownames_to_column()
names(df_ins.bee_taxa)[names(df_ins.bee_taxa)=="rowname"] <- "Row.names"
df_ins.bee_all_upset2 <- full_join(df_ins.bee_all_upset, df_ins.bee_taxa, by="Row.names")

# create the plot
df_ins.bee_all_upset2 <- df_ins.bee_all_upset2[, c(1,2,3,4,5,6,7,8,9,10,11,12)]    
tissues <- colnames(df_ins.bee_all_upset2)[2:3]
df_ins.bee_all_upset2[tissues] = df_ins.bee_all_upset2[tissues] == 1
t(head(df_ins.bee_all_upset2[tissues], 3))
df_ins.bee_all_upset2[is.na(df_ins.bee_all_upset2)] <- "unclassified"

Upset_ins_bee <- ComplexUpset::upset(
  df_ins.bee_all_upset2,
  tissues,
  name = 'Tissue type',
  width_ratio=0.1)
Upset_ins_bee2 <- Upset_ins_bee + theme(axis.title = element_text(size=20), axis.text = element_text(size=15)) 

Fig4A <- plot_grid(Upset_hand_bee2, Upset_ins_bee2, ncol=1, scale = 0.9)
Fig4A

#### compare Faith's phylogenetic diversity between treatments (Fig. 4B) ####
seqs.gyrB <- refseq(gyrB.rel.phy)
alignment.gyrB <- AlignTranslation(seqs.gyrB, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.gyrB <- phyDat(as(alignment.gyrB, "matrix"), type="DNA")
dm.gyrB <- dist.ml(phang.align.gyrB)
treeNJ.gyrB <- NJ(dm.gyrB) 
is.rooted(treeNJ.gyrB)
gyrB.tree <- root(treeNJ.gyrB, outgroup = 1, resolve.root = T)
is.rooted(gyrB.tree)
dgyrB.pt <- phyloseq(tax_table(gyrB.rel.phy), sample_data(gyrB.rel.phy),
                     otu_table(gyrB.rel.phy, taxa_are_rows = FALSE), 
                     refseq(gyrB.rel.phy), phy_tree(gyrB.tree))

seed.dgyrB.pt <- subset_samples(dgyrB.pt, tissue_type == "seed_pooled")
seed.dgyrB.pt <- prune_taxa(taxa_sums(seed.dgyrB.pt)>0, seed.dgyrB.pt)

seed.pt.asv.df <- as.data.frame(seed.dgyrB.pt@otu_table)
seed.gyrB.tree <- phy_tree(seed.dgyrB.pt)
seed.gyrB.pd <- pd(seed.pt.asv.df,seed.gyrB.tree, include.root = T)

seed.pt.meta <- as.data.frame(seed.dgyrB.pt@sam_data)
seed.pt.meta %<>% data.frame()
seed.pt.meta$Faiths_PD <- seed.gyrB.pd$PD
seed.pt.meta$treatment <- factor(seed.pt.meta$treatment, levels = c("hand", "insect"))

seed.pd.plot <- seed.pt.meta %>%
  ggplot(aes(x=treatment, y=Faiths_PD)) +
  geom_boxplot() +
  geom_jitter(size=4, width = 0.25) +
  theme_light()+
  theme(axis.title = element_text(size=30), 
        axis.text = element_text(size = 25),) +
  labs(x="Pollination treatment", y="Faith's PD") +
  scale_x_discrete(labels=c("hand" = "Hand", "insect" = "Hand + bee"))
seed.pd.plot 

Fig4B <- seed.pd.plot + stat_compare_means(method = "t.test")
Fig4B

t.test(Faiths_PD~treatment, data = seed.pt.meta)

#### compare community composition between treatments (Fig. 4C) ####
seed.gyrB.uni <- seed.dgyrB.pt %>% phyloseq::distance("wunifrac") %>% sqrt
seed.gyrB.nmds <- metaMDS(seed.gyrB.uni, trymax=200, parallel=10, k=2) 
seed.gyrB.nmds.dat <- scores(seed.gyrB.nmds, display = "site") %>% data.frame %>% 
  rownames_to_column("sampID") %>%
  full_join(sample_data(seed.dgyrB.pt) %>% data.frame %>% rownames_to_column("sampID"))

Fig4C <- ggplot(seed.gyrB.nmds.dat,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(color=treatment),size=5)+
  scale_color_discrete(type = c("#32648EFF","#74D055FF"),
                       labels=c("hand" = "Hand", "insect" = "Hand + bee"))+
  stat_ellipse(aes(color=treatment), lwd=2, linetype = 2) +
  ggthemes::theme_few() +
  labs(fill="Pollination\ntreatment", color = "Pollination\ntreatment") +
  theme(axis.title = element_text(size=30), axis.text = element_text(size=25), 
        legend.title = element_text(size=30), legend.text = element_text(size=25),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))
Fig4C

# PerMANOVA
seed.gyrB.dist <- seed.dgyrB.pt %>% phyloseq::distance("wunifrac")  
seed.perm <- adonis(seed.gyrB.dist~sample_data(seed.dgyrB.pt)$treatment)
print(seed.perm$aov.tab)

#### visualize taxonomic composition in each treatment (Fig. S5A) ####
seed.merged.phy <- dgyrB.pt %>%
  aggregate_top_taxa2(., 20, "Genus")

brewerPlus <- distinct_palette()

FigS5A <- plot_bar(seed.merged.phy, x="Sample",fill = "Genus") +
  geom_bar(aes(fill=Genus), stat="identity", position="stack", color="gray33") +
  facet_grid(~treatment, scales = "free_x", space = "free_x") +
  theme(legend.position = "right",legend.text = element_text(size=15), 
        axis.text.x = element_blank(), axis.text.y = element_text(size=25),
        axis.title = element_text(size=30), 
        strip.text = element_text(size=30))+  
  guides(fill=guide_legend(ncol = 1)) +
  labs(y="Relative abundance (%)") +
  scale_fill_manual(values = brewerPlus)
FigS5A

#### visualize the differences in beta-dispersion between treatments (Fig. S5B) ####
treatments <- factor(c("hand", "hand", "hand", "insect", "hand", "hand", "hand","insect","insect",
                   "hand", "hand", "hand","insect","insect","hand","hand","insect","insect"))

gyrB.dist <- vegdist(otu_table(seed.dgyrB.pt))
gyrB.disper <- betadisper(gyrB.dist, treatments)

betadisper_distances <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(distances = unlist(model$distances))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# getting eigenvalues out of betadisper() object
betadisper_eigenvalue <- function(model){
  temp <- data.frame(eig = unlist(model$eig))
  temp$PCoA <- row.names(temp)
  row.names(temp) <- NULL
  return(temp)
}

# getting the eigenvectors out of a betadisper() object
betadisper_eigenvector <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(unlist(model$vectors))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# get centroids
betadisper_centroids <- function(model){
  temp <- data.frame(unlist(model$centroids))
  temp$group <- row.names(temp)
  temp <- dplyr::select(temp, group, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# betadisper data
get_betadisper_data <- function(model){
  temp <- list(distances = betadisper_distances(model),
               eigenvalue = betadisper_eigenvalue(model),
               eigenvector = betadisper_eigenvector(model),
               centroids = betadisper_centroids(model))
  return(temp)
}

# get betadisper data ####
gyrb.disper.dat <- get_betadisper_data(gyrB.disper)

# do some transformations on the data
gyrb.disper.dat$eigenvalue <- mutate(gyrb.disper.dat$eigenvalue, percent = eig/sum(eig))

# make the nicer boxplot
gyrB.disper.dat.df <- gyrb.disper.dat$distances # need to join with metadata
gyrB.disper.dat.df

Figs5B<- gyrB.disper.dat.df %>%
  ggplot(aes(group, distances)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, 
               position = position_dodge(width = 0.55)) +
  stat_compare_means()+
  scale_fill_viridis(option = "C", discrete = T) +
  scale_color_viridis(option = "C", discrete = T) +
  ylab('Distance to centroid') +
  theme(legend.position = 'right', axis.title = element_text(size=20), 
        axis.text = element_text(size = 15)) +
  xlab('Sampling tissue') 
Figs5B
