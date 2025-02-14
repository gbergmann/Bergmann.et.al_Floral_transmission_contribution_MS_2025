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
library(UpSetR)
library(ComplexUpset)

#### loading data ####
bact.rich.rel.phy <- readRDS("data_output_seq/Field_survey/clean.bact.rich.rel.phy.rds")
bact.rel.phy <- readRDS("data_output_seq/Field_survey/clean.bact.rel.phy.rds")

#### alpha diversity (observed richness and Faith's PD) ####
# observed richness (Figure 1A)
bact.rich.rel.meta <- sample_data(bact.rich.rel.phy) %>% data.frame() %>% rownames_to_column(var = "SampID")
bact.rich.rel.asv <- otu_table(bact.rich.rel.phy) %>% data.frame() 
bact.rel.asv.binary <- bact.rel.asv %>%
  mutate_if(is.numeric, ~1 * (. != 0))
bact.rel.asv.binary <- bact.rel.asv.binary %>%
  dplyr::mutate(sum = rowSums(across(where(is.numeric))))
bact.rel.asv.binary <- rownames_to_column(bact.rel.asv.binary, var = "SampID")
obs.rich <- bact.rel.asv.binary %>%
  select(SampID, sum) %>%
  merge(x=., y=bact.rel.meta, by="SampID")
obs.rich$tissue_type <- factor(obs.rich$tissue_type, levels = c("stigma", "seed_pooled"))
obs.rich %>%
  group_by(tissue_type) %>%
  dplyr::summarize(mean_rich = mean(sum))

richness.aov <- aov(sum~tissue_type+Updated_Field_Name, data = obs.rich)
summary(richness.aov)
richness.ttest <- t.test(sum ~ tissue_type, data = obs.rich)
richness.ttest 

# creating the boxplot
Fig1A <- obs.rich %>%
  ggplot(aes(x=tissue_type, y=sum)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  facet_wrap(~Updated_Field_Name, nrow = 1) +
  ylim(0,350)+
  stat_compare_means(label = "p.format", 
                     label.x.npc = 0.7, label.y.npc = 0.95) +
  labs(x="Tissue type", y="Observed richness") +
  scale_x_discrete(labels=c("stigma" = "Stigma", "seed_pooled" = "Pooled\nseeds")) +
  theme_light() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15))
Fig1A

# phylogenetic alpha diversity
# generate phylogenetic tree
seqs.bact <- refseq(bact.rich.rel.phy)
alignment.bact <- AlignTranslation(seqs.bact, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.bact <- phyDat(as(alignment.bact, "matrix"), type="DNA")
dm.bact <- dist.ml(phang.align.bact)
treeNJ.bact <- NJ(dm.bact)# Note, tip order != sequence order 
is.rooted(treeNJ.bact)
bact.tree <- root(treeNJ.bact, outgroup = 1, resolve.root = T)
is.rooted(bact.tree)
dbact.pt <- phyloseq(tax_table(bact.rich.rel.phy), sample_data(bact.rich.rel.phy),
                     otu_table(bact.rich.rel.phy, taxa_are_rows = FALSE), refseq(bact.rich.rel.phy), phy_tree(bact.tree))
saveRDS(dbact.pt, "output/analysis/clean.bact.rich.rel.phylo.tree.rds")

# calculate Faith's PD in picante, merge with metadata
bact.pd <- pd(dbact.pt@otu_table, dbact.pt@phy_tree, include.root = T)
dbact.meta <- as.data.frame(dbact.pt@sam_data)
dbact.meta$Faiths_PD <- bact.pd$PD
dbact.meta$tissue_type <- factor(dbact.meta$tissue_type, levels = c("stigma", "seed_pooled"))

pd.plot <- dbact.meta %>%
  ggplot(aes(x=tissue_type, y=Faiths_PD)) +
  geom_boxplot() +
  geom_jitter(width=0.25) +
  facet_wrap(~Updated_Field_Name, nrow=1) +
  labs(x="Tissue type", y="Faith's PD") +
  scale_x_discrete(labels=c("stigma" = "Stigma", "seed_pooled" = "Pooled\nseeds")) +
  theme_light() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15))

FigS1 <- pd.plot + stat_compare_means(label = "p.format", 
                                         label.x.npc = 0.7, label.y.npc = 0.9)

#### NMDS and PerMANOVA ####
# visualizing differences between tissue types (Fig. 1B)
bact.bray <- dbact.pt %>% phyloseq::distance("wunifrac") %>% sqrt
bact.nmds <- metaMDS(bact.bray, trymax=200, parallel=10, k=2) 
bact.nmds.dat <- scores(bact.nmds, display = "site") %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(dbact.pt) %>% data.frame %>% rownames_to_column("sampID"))

nmds.plot <- ggplot(bact.nmds.dat,aes(x=NMDS1,y=NMDS2, label = sampID))+ 
  geom_point(aes(shape=Updated_Field_Name, color=tissue_type),size=3)+
  geom_text() +
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few() +
  labs(fill="Tissue type", color = "Tissue type", shape="Field")
nmds.plot #outliers: SJ1_1_8, Co2_1_8, Co1_1_3
ggplotly(nmds.plot)

# removing outliers, redoing the visualization
dbact.pt2 <- subset_samples(dbact.pt, sample_names(dbact.pt) !="Co1_1_3" & 
                              sample_names(dbact.pt) !="SJ1_1_8" &
                              sample_names(dbact.pt) != "Co2_1_8")
dbact.pt2 <- prune_taxa(taxa_sums(dbact.pt2)>0, dbact.pt2)
sample_data(dbact.pt2)$tissue_type <- factor(sample_data(dbact.pt2)$tissue_type,
                                             levels = c("stigma", "seed_pooled"))

bact.bray2 <- dbact.pt2 %>% phyloseq::distance("wunifrac") %>% sqrt
bact.nmds2 <- metaMDS(bact.bray2, trymax=200, parallel=10, k=2) 
bact.nmds.dat2 <- scores(bact.nmds2, display = "site") %>% 
  data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(dbact.pt2) %>% data.frame %>% rownames_to_column("sampID"))

Fig1B <- ggplot(bact.nmds.dat2,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(shape=Updated_Field_Name, color=tissue_type),size=4)+
  scale_color_discrete(labels = c("Stigmas", "Pooled seeds")) +
  stat_ellipse(aes(color=tissue_type),lwd = 2)+
  scale_shape_manual(values = c(15,16,17,8)) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  ggthemes::theme_few() +
  labs(fill="Tissue type", color = "Tissue type", shape="Field") +
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size=16),
        axis.title.y = element_blank(),
        legend.position = "bottom", 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16)) +
  guides(fill=guide_legend(nrow = 2)) 
Fig1B

# Looking at potential field effects within each tissue type (Fig. 1C)
stigma.pt.phy <- subset_samples(dbact.pt2, tissue_type=="stigma")
stigma.pt.phy <- prune_taxa(taxa_sums(stigma.pt.phy)>0, stigma.pt.phy)
stigma.uni <- stigma.pt.phy %>% phyloseq::distance("wunifrac") %>% sqrt
stigma.uni.nmds <- metaMDS(stigma.uni, trymax=200, parallel=10, k=2) 
stigma.uni.nmds.dat <- scores(stigma.uni.nmds, display = "site") %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(stigma.pt.phy) %>% data.frame %>% rownames_to_column("sampID"))

stigma.uni.nmds.plot <- ggplot(stigma.uni.nmds.dat,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(color=Updated_Field_Name),size=4)+
  stat_ellipse(aes(color=Updated_Field_Name),lwd = 2)+
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few() +
  labs(color = "Field") +
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size=16),
        axis.title.y = element_blank(),
        legend.position = "right", 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))
stigma.uni.nmds.plot

seed.pt.phy <- subset_samples(dbact.pt2, tissue_type=="seed_pooled")
seed.pt.phy <- prune_taxa(taxa_sums(seed.pt.phy)>0, seed.pt.phy)
seed.uni <- seed.pt.phy %>% phyloseq::distance("wunifrac") %>% sqrt
seed.uni.nmds <- metaMDS(seed.uni, trymax=200, parallel=10, k=2) 
seed.uni.nmds.dat <- scores(seed.uni.nmds, display = "site") %>% data.frame %>% rownames_to_column("sampID") %>%
  full_join(sample_data(seed.pt.phy) %>% data.frame %>% rownames_to_column("sampID"))

seed.uni.nmds.plot <- ggplot(seed.uni.nmds.dat,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(color=Updated_Field_Name),size=4)+
  stat_ellipse(aes(color=Updated_Field_Name),lwd = 2)+
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few() +
  labs(color = "Field") +
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size=16),
        axis.title.y = element_blank(),
        legend.position = "right", 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))
seed.uni.nmds.plot

Fig1C <- ggpubr::ggarrange(stigma.uni.nmds.plot,seed.uni.nmds.plot, 
                                            nrow=1, common.legend = T, legend = "bottom")
Fig1C

# PerMANOVAs between tissue types and by field within each tissue type
  # by tissue type
bact.dist <- dbact.pt2 %>% phyloseq::distance("wunifrac") 
perm2 <- adonis2(bact.dist~sample_data(dbact.pt2)$tissue_type)
perm2

pairwise.adonis(bact.dist, sample_data(dbact.pt2)$tissue_type)

  # by field for stigmas
stigma.uni.dist <- stigma.pt.phy %>% phyloseq::distance("wunifrac")
stigma.perm <- adonis2(stigma.uni.dist~sample_data(stigma.pt.phy)$Updated_Field_Name)
stigma.perm

  # by field for seeds
seed.uni.dist <- seed.pt.phy %>% phyloseq::distance("wunifrac")
seed.perm <- adonis2(seed.uni.dist~sample_data(seed.pt.phy)$Updated_Field_Name)
seed.perm

#### overlap between seeds and stigmas across all samples (Fig. 2A) ####
sample_data(bact.rel.phy)$tissue_type <- factor(sample_data(bact.rel.phy)$tissue_type,
                                                levels = c("stigma", "seed_pooled"))
variable1 = as.character(get_variable(bact.rel.phy, "tissue_type"))
sample_data(bact.rel.phy)$tissue <- mapply(paste0, variable1, collapse = "_")

merged.16S_all <- merge_samples(bact.rel.phy, "tissue")
#Fix continuous variable to discrete value
sample_data(merged.16S_all)$tissue <- factor(sample_names(merged.16S_all),
                                             levels = c("stigma","seed_pooled"))
#Extract count datasets from phyloseq
count_16S_all <- t(otu_table(merged.16S_all))
#Transform in binary matrix
count_16S_all[count_16S_all> 0] <- 1
#Convert to dataframe
df_16S_all <- as.data.frame(count_16S_all)
#Read abundance per ASV
rel.abu <- data.frame(taxa_sums(bact.rel.phy))
rel.abu$log <- log10(rel.abu$taxa_sums.bact.rel.phy.)
#Merge ASV dataframe and ASV abundance
df_16S_all_upset <-merge.data.frame(df_16S_all, rel.abu, by="row.names",
                                    all.x=TRUE)
#Append taxonomy to upset dataframe
df_16S_taxa <- tax_table(bact.rel.phy) %>% data.frame() %>% rownames_to_column()
names(df_16S_taxa)[names(df_16S_taxa)=="rowname"] <- "Row.names"

df_16S_all_upset2 <- full_join(df_16S_all_upset, df_16S_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_16S_all <- UpSetR::upset(df_16S_all_upset2,
                               nsets = 2, #adjust to number of columns to compare, here 2 tissue columns
                               nintersects = NA,
                               keep.order = T,
                               order.by = "freq",
                               mainbar.y.label = "Number of ASVs",
                               sets.x.label="total tissue\nASV richness",
                               decreasing ="TRUE", 
                               text.scale = 2
)
Upset_16S_all 

# preparing something fancier with ComplexUpset
df_16S_all_upset3 <- df_16S_all_upset2[, c(1,3,2,4,5,6,7,8,9,10,11)]    
tissues <- colnames(df_16S_all_upset3)[2:3]
df_16S_all_upset3[tissues] = df_16S_all_upset3[tissues] == 1
t(head(df_16S_all_upset3[tissues], 2))
df_16S_all_upset3 <- df_16S_all_upset3 %>%
  filter(!is.na(stigma), !is.na(seed_pooled))

Upset_16S_all2 <- ComplexUpset::upset(
  df_16S_all_upset3,
  tissues,
  base_annotations=list(
    'Intersection size'=intersection_size(
      text=list(size=6))),
  min_size=10,
  width_ratio=0.1
)

Fig2A <- Upset_16S_all2 + theme(axis.title = element_text(size=20), axis.text = element_text(size=15))
Fig2A 

#### Pairwise comparisons of taxa overlaps (Fig. 2B, Fig. S2) ####
bact.rel.meta <- sample_data(bact.rel.phy) %>% data.frame()
bact.rel.meta2 <- bact.rel.meta
bact.rel.meta2$Plant_rep <- c('1','10','3','4','5','6','8',
                              '1','10','3','4','5','6','7','8',
                              '1','2','3','4','5','6','8','9',
                              '1','10','2','4','6','7','8','9',
                              '1','2','3','4','5','6','7','8',
                              '1','10','2','3','4','5','6','7','8','9',
                              '10','2','3','4','5','6','7',
                              '1','2','3','5','6','9')
sample_data(bact.rel.phy) <- bact.rel.meta2

# SJ1.1
SJ1.1.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ1"&Plant_rep == "1") 
SJ1.1.phy <- prune_taxa(taxa_sums(SJ1.1.phy)>0, SJ1.1.phy)

#Extract count datasets from phyloseq
count_SJ1.1_all <- t(otu_table(SJ1.1.phy))
#Transform in binary matrix
count_SJ1.1_all[count_SJ1.1_all> 0] <- 1
#Convert to dataframe
df_SJ1.1_all <- as.data.frame(count_SJ1.1_all)
#Read abundance per ASV
SJ1.1_rel.abu <- data.frame(taxa_sums(SJ1.1.phy))
SJ1.1_rel.abu$log <- log10(SJ1.1_rel.abu$taxa_sums.SJ1.1.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ1.1_all_upset <-merge.data.frame(df_SJ1.1_all, SJ1.1_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ1.1_taxa <- tax_table(SJ1.1.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ1.1_taxa)[names(df_SJ1.1_taxa)=="rowname"] <- "Row.names"

df_SJ1.1_all_upset2 <- full_join(df_SJ1.1_all_upset, df_SJ1.1_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ1.1_all <- UpSetR::upset(df_SJ1.1_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE")
Upset_SJ1.1_all # 17 ASVs out of 215 stigma, 72 seed

# SJ1.2
SJ1.2.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ1"&Plant_rep == "2") 
SJ1.2.phy <- prune_taxa(taxa_sums(SJ1.2.phy)>0, SJ1.2.phy)

#Extract count datasets from phyloseq
count_SJ1.2_all <- t(otu_table(SJ1.2.phy))
#Transform in binary matrix
count_SJ1.2_all[count_SJ1.2_all> 0] <- 1
#Convert to dataframe
df_SJ1.2_all <- as.data.frame(count_SJ1.2_all)
#Read abundance per ASV
SJ1.2_rel.abu <- data.frame(taxa_sums(SJ1.2.phy))
SJ1.2_rel.abu$log <- log10(SJ1.2_rel.abu$taxa_sums.SJ1.2.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ1.2_all_upset <-merge.data.frame(df_SJ1.2_all, SJ1.2_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ1.2_taxa <- tax_table(SJ1.2.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ1.2_taxa)[names(df_SJ1.2_taxa)=="rowname"] <- "Row.names"

df_SJ1.2_all_upset2 <- full_join(df_SJ1.2_all_upset, df_SJ1.2_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ1.2_all <- UpSetR::upset(df_SJ1.2_all_upset2,
                                 nsets = 2,
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE")
Upset_SJ1.2_all # 11 ASVs out of 176 stigma, 87 seed

# SJ1.3 
SJ1.3.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ1"&Plant_rep == "3") 
SJ1.3.phy <- prune_taxa(taxa_sums(SJ1.3.phy)>0, SJ1.3.phy)

#Extract count datasets from phyloseq
count_SJ1.3_all <- t(otu_table(SJ1.3.phy))
#Transform in binary matrix
count_SJ1.3_all[count_SJ1.3_all> 0] <- 1
#Convert to dataframe
df_SJ1.3_all <- as.data.frame(count_SJ1.3_all)
#Read abundance per ASV
SJ1.3_rel.abu <- data.frame(taxa_sums(SJ1.3.phy))
SJ1.3_rel.abu$log <- log10(SJ1.3_rel.abu$taxa_sums.SJ1.3.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ1.3_all_upset <-merge.data.frame(df_SJ1.3_all, SJ1.3_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ1.3_taxa <- tax_table(SJ1.3.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ1.3_taxa)[names(df_SJ1.3_taxa)=="rowname"] <- "Row.names"

df_SJ1.3_all_upset2 <- full_join(df_SJ1.3_all_upset, df_SJ1.3_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ1.3_all <- UpSetR::upset(df_SJ1.3_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE")
Upset_SJ1.3_all # 41 ASVs out of 225 stigma, 222 seed

# SJ1.4
SJ1.4.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ1"&Plant_rep == "4") 
SJ1.4.phy <- prune_taxa(taxa_sums(SJ1.4.phy)>0, SJ1.4.phy)

#Extract count datasets from phyloseq
count_SJ1.4_all <- t(otu_table(SJ1.4.phy))
#Transform in binary matrix
count_SJ1.4_all[count_SJ1.4_all> 0] <- 1
#Convert to dataframe
df_SJ1.4_all <- as.data.frame(count_SJ1.4_all)
#Read abundance per ASV
SJ1.4_rel.abu <- data.frame(taxa_sums(SJ1.4.phy))
SJ1.4_rel.abu$log <- log10(SJ1.4_rel.abu$taxa_sums.SJ1.4.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ1.4_all_upset <-merge.data.frame(df_SJ1.4_all, SJ1.4_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ1.4_taxa <- tax_table(SJ1.4.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ1.4_taxa)[names(df_SJ1.4_taxa)=="rowname"] <- "Row.names"

df_SJ1.4_all_upset2 <- full_join(df_SJ1.4_all_upset, df_SJ1.4_taxa, by="Row.names")

# Upset graph will all comparisons
Upset_SJ1.4_all <- UpSetR::upset(df_SJ1.4_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE")
Upset_SJ1.4_all # 13 ASVs out of 99 stigma, 45 seed

# SJ1.5 
SJ1.5.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ1"&Plant_rep == "5") 
SJ1.5.phy <- prune_taxa(taxa_sums(SJ1.5.phy)>0, SJ1.5.phy)

#Extract count datasets from phyloseq
count_SJ1.5_all <- t(otu_table(SJ1.5.phy))
#Transform in binary matrix
count_SJ1.5_all[count_SJ1.5_all> 0] <- 1
#Convert to dataframe
df_SJ1.5_all <- as.data.frame(count_SJ1.5_all)
#Read abundance per ASV
SJ1.5_rel.abu <- data.frame(taxa_sums(SJ1.5.phy))
SJ1.5_rel.abu$log <- log10(SJ1.5_rel.abu$taxa_sums.SJ1.5.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ1.5_all_upset <-merge.data.frame(df_SJ1.5_all, SJ1.5_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ1.5_taxa <- tax_table(SJ1.5.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ1.5_taxa)[names(df_SJ1.5_taxa)=="rowname"] <- "Row.names"

df_SJ1.5_all_upset2 <- full_join(df_SJ1.5_all_upset, df_SJ1.5_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ1.5_all <- UpSetR::upset(df_SJ1.5_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE")
Upset_SJ1.5_all # 7 ASVs out of 37 stigma, 78 seed

# SJ1.6
SJ1.6.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ1"&Plant_rep == "6") 
SJ1.6.phy <- prune_taxa(taxa_sums(SJ1.6.phy)>0, SJ1.6.phy)

#Extract count datasets from phyloseq
count_SJ1.6_all <- t(otu_table(SJ1.6.phy))
#Transform in binary matrix
count_SJ1.6_all[count_SJ1.6_all> 0] <- 1
#Convert to dataframe
df_SJ1.6_all <- as.data.frame(count_SJ1.6_all)
#Read abundance per ASV
SJ1.6_rel.abu <- data.frame(taxa_sums(SJ1.6.phy))
SJ1.6_rel.abu$log <- log10(SJ1.6_rel.abu$taxa_sums.SJ1.6.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ1.6_all_upset <-merge.data.frame(df_SJ1.6_all, SJ1.6_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ1.6_taxa <- tax_table(SJ1.6.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ1.6_taxa)[names(df_SJ1.6_taxa)=="rowname"] <- "Row.names"

df_SJ1.6_all_upset2 <- full_join(df_SJ1.6_all_upset, df_SJ1.6_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ1.6_all <- UpSetR::upset(df_SJ1.6_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_SJ1.6_all # 16 ASVs out of 251 stigma, 114 seed

# SJ1.7
SJ1.7.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ1"&Plant_rep == "7") 
SJ1.7.phy <- prune_taxa(taxa_sums(SJ1.7.phy)>0, SJ1.7.phy)

#Extract count datasets from phyloseq
count_SJ1.7_all <- t(otu_table(SJ1.7.phy))
#Transform in binary matrix
count_SJ1.7_all[count_SJ1.7_all> 0] <- 1
#Convert to dataframe
df_SJ1.7_all <- as.data.frame(count_SJ1.7_all)
#Read abundance per ASV
SJ1.7_rel.abu <- data.frame(taxa_sums(SJ1.7.phy))
SJ1.7_rel.abu$log <- log10(SJ1.7_rel.abu$taxa_sums.SJ1.7.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ1.7_all_upset <-merge.data.frame(df_SJ1.7_all, SJ1.7_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ1.7_taxa <- tax_table(SJ1.7.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ1.7_taxa)[names(df_SJ1.7_taxa)=="rowname"] <- "Row.names"

df_SJ1.7_all_upset2 <- full_join(df_SJ1.7_all_upset, df_SJ1.7_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ1.7_all <- UpSetR::upset(df_SJ1.7_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_SJ1.7_all # 9 ASVs out of 47 stigma, 140 seed

# SJ1.8
SJ1.8.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ1"&Plant_rep == "8") 
SJ1.8.phy <- prune_taxa(taxa_sums(SJ1.8.phy)>0, SJ1.8.phy)

#Extract count datasets from phyloseq
count_SJ1.8_all <- t(otu_table(SJ1.8.phy))
#Transform in binary matrix
count_SJ1.8_all[count_SJ1.8_all> 0] <- 1
#Convert to dataframe
df_SJ1.8_all <- as.data.frame(count_SJ1.8_all)
#Read abundance per ASV
SJ1.8_rel.abu <- data.frame(taxa_sums(SJ1.8.phy))
SJ1.8_rel.abu$log <- log10(SJ1.8_rel.abu$taxa_sums.SJ1.8.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ1.8_all_upset <-merge.data.frame(df_SJ1.8_all, SJ1.8_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ1.8_taxa <- tax_table(SJ1.8.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ1.8_taxa)[names(df_SJ1.8_taxa)=="rowname"] <- "Row.names"

df_SJ1.8_all_upset2 <- full_join(df_SJ1.8_all_upset, df_SJ1.8_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ1.8_all <- UpSetR::upset(df_SJ1.8_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_SJ1.8_all # 2 ASVs out of 25 stigma, 137 seed


# SJ2.2 
SJ2.2.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ2"&Plant_rep == "2") 
SJ2.2.phy <- prune_taxa(taxa_sums(SJ2.2.phy)>0, SJ2.2.phy)

#Extract count datasets from phyloseq
count_SJ2.2_all <- t(otu_table(SJ2.2.phy))
#Transform in binary matrix
count_SJ2.2_all[count_SJ2.2_all> 0] <- 1
#Convert to dataframe
df_SJ2.2_all <- as.data.frame(count_SJ2.2_all)
#Read abundance per ASV
SJ2.2_rel.abu <- data.frame(taxa_sums(SJ2.2.phy))
SJ2.2_rel.abu$log <- log10(SJ2.2_rel.abu$taxa_sums.SJ2.2.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ2.2_all_upset <-merge.data.frame(df_SJ2.2_all, SJ2.2_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ2.2_taxa <- tax_table(SJ2.2.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ2.2_taxa)[names(df_SJ2.2_taxa)=="rowname"] <- "Row.names"

df_SJ2.2_all_upset2 <- full_join(df_SJ2.2_all_upset, df_SJ2.2_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ2.2_all <- UpSetR::upset(df_SJ2.2_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_SJ2.2_all # 12 ASVs out of 287 stigma, 75 seed

# SJ2.3
SJ2.3.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ2"&Plant_rep == "3") 
SJ2.3.phy <- prune_taxa(taxa_sums(SJ2.3.phy)>0, SJ2.3.phy)

#Extract count datasets from phyloseq
count_SJ2.3_all <- t(otu_table(SJ2.3.phy))
#Transform in binary matrix
count_SJ2.3_all[count_SJ2.3_all> 0] <- 1
#Convert to dataframe
df_SJ2.3_all <- as.data.frame(count_SJ2.3_all)
#Read abundance per ASV
SJ2.3_rel.abu <- data.frame(taxa_sums(SJ2.3.phy))
SJ2.3_rel.abu$log <- log10(SJ2.3_rel.abu$taxa_sums.SJ2.3.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ2.3_all_upset <-merge.data.frame(df_SJ2.3_all, SJ2.3_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ2.3_taxa <- tax_table(SJ2.3.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ2.3_taxa)[names(df_SJ2.3_taxa)=="rowname"] <- "Row.names"

df_SJ2.3_all_upset2 <- full_join(df_SJ2.3_all_upset, df_SJ2.3_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ2.3_all <- UpSetR::upset(df_SJ2.3_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_SJ2.3_all # 11 ASVs out of 121 stigma, 66 seed

# SJ2.5
SJ2.5.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ2"&Plant_rep == "5") 
SJ2.5.phy <- prune_taxa(taxa_sums(SJ2.5.phy)>0, SJ2.5.phy)

#Extract count datasets from phyloseq
count_SJ2.5_all <- t(otu_table(SJ2.5.phy))
#Transform in binary matrix
count_SJ2.5_all[count_SJ2.5_all> 0] <- 1
#Convert to dataframe
df_SJ2.5_all <- as.data.frame(count_SJ2.5_all)
#Read abundance per ASV
SJ2.5_rel.abu <- data.frame(taxa_sums(SJ2.5.phy))
SJ2.5_rel.abu$log <- log10(SJ2.5_rel.abu$taxa_sums.SJ2.5.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ2.5_all_upset <-merge.data.frame(df_SJ2.5_all, SJ2.5_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ2.5_taxa <- tax_table(SJ2.5.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ2.5_taxa)[names(df_SJ2.5_taxa)=="rowname"] <- "Row.names"

df_SJ2.5_all_upset2 <- full_join(df_SJ2.5_all_upset, df_SJ2.5_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ2.5_all <- UpSetR::upset(df_SJ2.5_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_SJ2.5_all # 18 ASVs out of 106 stigma, 130 seed

# SJ2.6
SJ2.6.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "SJ2"&Plant_rep == "6") 
SJ2.6.phy <- prune_taxa(taxa_sums(SJ2.6.phy)>0, SJ2.6.phy)

#Extract count datasets from phyloseq
count_SJ2.6_all <- t(otu_table(SJ2.6.phy))
#Transform in binary matrix
count_SJ2.6_all[count_SJ2.6_all> 0] <- 1
#Convert to dataframe
df_SJ2.6_all <- as.data.frame(count_SJ2.6_all)
#Read abundance per ASV
SJ2.6_rel.abu <- data.frame(taxa_sums(SJ2.6.phy))
SJ2.6_rel.abu$log <- log10(SJ2.6_rel.abu$taxa_sums.SJ2.6.phy.)
#Merge ASV dataframe and ASV abundance
df_SJ2.6_all_upset <-merge.data.frame(df_SJ2.6_all, SJ2.6_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_SJ2.6_taxa <- tax_table(SJ2.6.phy) %>% data.frame() %>% rownames_to_column()
names(df_SJ2.6_taxa)[names(df_SJ2.6_taxa)=="rowname"] <- "Row.names"

df_SJ2.6_all_upset2 <- full_join(df_SJ2.6_all_upset, df_SJ2.6_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_SJ2.6_all <- UpSetR::upset(df_SJ2.6_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_SJ2.6_all # 18 ASVs out of 210 stigma, 50 seed

# Co1.1
Co1.1.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co1"&Plant_rep == "1") 
Co1.1.phy <- prune_taxa(taxa_sums(Co1.1.phy)>0, Co1.1.phy)

#Extract count datasets from phyloseq
count_Co1.1_all <- t(otu_table(Co1.1.phy))
#Transform in binary matrix
count_Co1.1_all[count_Co1.1_all> 0] <- 1
#Convert to dataframe
df_Co1.1_all <- as.data.frame(count_Co1.1_all)
#Read abundance per ASV
Co1.1_rel.abu <- data.frame(taxa_sums(Co1.1.phy))
Co1.1_rel.abu$log <- log10(Co1.1_rel.abu$taxa_sums.Co1.1.phy.)
#Merge ASV dataframe and ASV abundance
df_Co1.1_all_upset <-merge.data.frame(df_Co1.1_all, Co1.1_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co1.1_taxa <- tax_table(Co1.1.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co1.1_taxa)[names(df_Co1.1_taxa)=="rowname"] <- "Row.names"

df_Co1.1_all_upset2 <- full_join(df_Co1.1_all_upset, df_Co1.1_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co1.1_all <- UpSetR::upset(df_Co1.1_all_upset2,
                                 nsets = 2,
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE",
                                 text.scale = 2)
Upset_Co1.1_all # 7 ASVs out of 128 stigma, 23 seed

# Co1.3
Co1.3.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co1"&Plant_rep == "3") 
Co1.3.phy <- prune_taxa(taxa_sums(Co1.3.phy)>0, Co1.3.phy)

#Extract count datasets from phyloseq
count_Co1.3_all <- t(otu_table(Co1.3.phy))
#Transform in binary matrix
count_Co1.3_all[count_Co1.3_all> 0] <- 1
#Convert to dataframe
df_Co1.3_all <- as.data.frame(count_Co1.3_all)
#Read abundance per ASV
Co1.3_rel.abu <- data.frame(taxa_sums(Co1.3.phy))
Co1.3_rel.abu$log <- log10(Co1.3_rel.abu$taxa_sums.Co1.3.phy.)
#Merge ASV dataframe and ASV abundance
df_Co1.3_all_upset <-merge.data.frame(df_Co1.3_all, Co1.3_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co1.3_taxa <- tax_table(Co1.3.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co1.3_taxa)[names(df_Co1.3_taxa)=="rowname"] <- "Row.names"

df_Co1.3_all_upset2 <- full_join(df_Co1.3_all_upset, df_Co1.3_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co1.3_all <- UpSetR::upset(df_Co1.3_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_Co1.3_all # 4 ASVs out of 26 stigma, 57 seed

# Co1.4
Co1.4.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co1"&Plant_rep == "4") 
Co1.4.phy <- prune_taxa(taxa_sums(Co1.4.phy)>0, Co1.4.phy)

#Extract count datasets from phyloseq
count_Co1.4_all <- t(otu_table(Co1.4.phy))
#Transform in binary matrix
count_Co1.4_all[count_Co1.4_all> 0] <- 1
#Convert to dataframe
df_Co1.4_all <- as.data.frame(count_Co1.4_all)
#Read abundance per ASV
Co1.4_rel.abu <- data.frame(taxa_sums(Co1.4.phy))
Co1.4_rel.abu$log <- log10(Co1.4_rel.abu$taxa_sums.Co1.4.phy.)
#Merge ASV dataframe and ASV abundance
df_Co1.4_all_upset <-merge.data.frame(df_Co1.4_all, Co1.4_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co1.4_taxa <- tax_table(Co1.4.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co1.4_taxa)[names(df_Co1.4_taxa)=="rowname"] <- "Row.names"

df_Co1.4_all_upset2 <- full_join(df_Co1.4_all_upset, df_Co1.4_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co1.4_all <- UpSetR::upset(df_Co1.4_all_upset2,
                                 nsets = 2,
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_Co1.4_all # 8 ASVs out of 129 stigma, 81 seed


# Co1.5
Co1.5.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co1"&Plant_rep == "5") 
Co1.5.phy <- prune_taxa(taxa_sums(Co1.5.phy)>0, Co1.5.phy)

#Extract count datasets from phyloseq
count_Co1.5_all <- t(otu_table(Co1.5.phy))
#Transform in binary matrix
count_Co1.5_all[count_Co1.5_all> 0] <- 1
#Convert to dataframe
df_Co1.5_all <- as.data.frame(count_Co1.5_all)
#Read abundance per ASV
Co1.5_rel.abu <- data.frame(taxa_sums(Co1.5.phy))
Co1.5_rel.abu$log <- log10(Co1.5_rel.abu$taxa_sums.Co1.5.phy.)
#Merge ASV dataframe and ASV abundance
df_Co1.5_all_upset <-merge.data.frame(df_Co1.5_all, Co1.5_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co1.5_taxa <- tax_table(Co1.5.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co1.5_taxa)[names(df_Co1.5_taxa)=="rowname"] <- "Row.names"

df_Co1.5_all_upset2 <- full_join(df_Co1.5_all_upset, df_Co1.5_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co1.5_all <- UpSetR::upset(df_Co1.5_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_Co1.5_all # 9 ASVs out of 65 stigma, 71 seed

# Co1.6
Co1.6.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co1"&Plant_rep == "6") 
Co1.6.phy <- prune_taxa(taxa_sums(Co1.6.phy)>0, Co1.6.phy)

#Extract count datasets from phyloseq
count_Co1.6_all <- t(otu_table(Co1.6.phy))
#Transform in binary matrix
count_Co1.6_all[count_Co1.6_all> 0] <- 1
#Convert to dataframe
df_Co1.6_all <- as.data.frame(count_Co1.6_all)
#Read abundance per ASV
Co1.6_rel.abu <- data.frame(taxa_sums(Co1.6.phy))
Co1.6_rel.abu$log <- log10(Co1.6_rel.abu$taxa_sums.Co1.6.phy.)
#Merge ASV dataframe and ASV abundance
df_Co1.6_all_upset <-merge.data.frame(df_Co1.6_all, Co1.6_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co1.6_taxa <- tax_table(Co1.6.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co1.6_taxa)[names(df_Co1.6_taxa)=="rowname"] <- "Row.names"

df_Co1.6_all_upset2 <- full_join(df_Co1.6_all_upset, df_Co1.6_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co1.6_all <- UpSetR::upset(df_Co1.6_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_Co1.6_all # 7 ASVs out of 68 stigma, 82 seed

# Co1.8
Co1.8.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co1"&Plant_rep == "8") 
Co1.8.phy <- prune_taxa(taxa_sums(Co1.8.phy)>0, Co1.8.phy)

#Extract count datasets from phyloseq
count_Co1.8_all <- t(otu_table(Co1.8.phy))
#Transform in binary matrix
count_Co1.8_all[count_Co1.8_all> 0] <- 1
#Convert to dataframe
df_Co1.8_all <- as.data.frame(count_Co1.8_all)
#Read abundance per ASV
Co1.8_rel.abu <- data.frame(taxa_sums(Co1.8.phy))
Co1.8_rel.abu$log <- log10(Co1.8_rel.abu$taxa_sums.Co1.8.phy.)
#Merge ASV dataframe and ASV abundance
df_Co1.8_all_upset <-merge.data.frame(df_Co1.8_all, Co1.8_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co1.8_taxa <- tax_table(Co1.8.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co1.8_taxa)[names(df_Co1.8_taxa)=="rowname"] <- "Row.names"

df_Co1.8_all_upset2 <- full_join(df_Co1.8_all_upset, df_Co1.8_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co1.8_all <- UpSetR::upset(df_Co1.8_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_Co1.8_all # 3 ASVs out of 197 stigma, 47 seed

# Co1.10
Co1.10.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co1"&Plant_rep == "10") 
Co1.10.phy <- prune_taxa(taxa_sums(Co1.10.phy)>0, Co1.10.phy)

#Extract count datasets from phyloseq
count_Co1.10_all <- t(otu_table(Co1.10.phy))
#Transform in binary matrix
count_Co1.10_all[count_Co1.10_all> 0] <- 1
#Convert to dataframe
df_Co1.10_all <- as.data.frame(count_Co1.10_all)
#Read abundance per ASV
Co1.10_rel.abu <- data.frame(taxa_sums(Co1.10.phy))
Co1.10_rel.abu$log <- log10(Co1.10_rel.abu$taxa_sums.Co1.10.phy.)
#Merge ASV dataframe and ASV abundance
df_Co1.10_all_upset <-merge.data.frame(df_Co1.10_all, Co1.10_rel.abu, by="row.names",
                                       all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co1.10_taxa <- tax_table(Co1.10.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co1.10_taxa)[names(df_Co1.10_taxa)=="rowname"] <- "Row.names"

df_Co1.10_all_upset2 <- full_join(df_Co1.10_all_upset, df_Co1.10_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co1.10_all <- UpSetR::upset(df_Co1.10_all_upset2,
                                  nsets = 2, 
                                  nintersects = NA,
                                  order.by = "freq",
                                  mainbar.y.label = "Number of ASVs",
                                  sets.x.label="total tissue ASV richness",
                                  decreasing ="TRUE", 
                                  text.scale = 2)
Upset_Co1.10_all # 12 ASVs out of 94 stigma, 59 seed

# Co2.1
Co2.1.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co2"&Plant_rep == "1") 
Co2.1.phy <- prune_taxa(taxa_sums(Co2.1.phy)>0, Co2.1.phy)

#Extract count datasets from phyloseq
count_Co2.1_all <- t(otu_table(Co2.1.phy))
#Transform in binary matrix
count_Co2.1_all[count_Co2.1_all> 0] <- 1
#Convert to dataframe
df_Co2.1_all <- as.data.frame(count_Co2.1_all)
#Read abundance per ASV
Co2.1_rel.abu <- data.frame(taxa_sums(Co2.1.phy))
Co2.1_rel.abu$log <- log10(Co2.1_rel.abu$taxa_sums.Co2.1.phy.)
#Merge ASV dataframe and ASV abundance
df_Co2.1_all_upset <-merge.data.frame(df_Co2.1_all, Co2.1_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co2.1_taxa <- tax_table(Co2.1.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co2.1_taxa)[names(df_Co2.1_taxa)=="rowname"] <- "Row.names"

df_Co2.1_all_upset2 <- full_join(df_Co2.1_all_upset, df_Co2.1_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co2.1_all <- UpSetR::upset(df_Co2.1_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_Co2.1_all # 10 ASVs out of 137 stigma, 43 seed

# Co2.2
Co2.2.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co2"&Plant_rep == "2") 
Co2.2.phy <- prune_taxa(taxa_sums(Co2.2.phy)>0, Co2.2.phy)

#Extract count datasets from phyloseq
count_Co2.2_all <- t(otu_table(Co2.2.phy))
#Transform in binary matrix
count_Co2.2_all[count_Co2.2_all> 0] <- 1
#Convert to dataframe
df_Co2.2_all <- as.data.frame(count_Co2.2_all)
#Read abundance per ASV
Co2.2_rel.abu <- data.frame(taxa_sums(Co2.2.phy))
Co2.2_rel.abu$log <- log10(Co2.2_rel.abu$taxa_sums.Co2.2.phy.)
#Merge ASV dataframe and ASV abundance
df_Co2.2_all_upset <-merge.data.frame(df_Co2.2_all, Co2.2_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co2.2_taxa <- tax_table(Co2.2.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co2.2_taxa)[names(df_Co2.2_taxa)=="rowname"] <- "Row.names"

df_Co2.2_all_upset2 <- full_join(df_Co2.2_all_upset, df_Co2.2_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co2.2_all <- UpSetR::upset(df_Co2.2_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_Co2.2_all # 7 ASVs out of 76 stigma, 45 seed

# Co2.4
Co2.4.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co2"&Plant_rep == "4") 
Co2.4.phy <- prune_taxa(taxa_sums(Co2.4.phy)>0, Co2.4.phy)

#Extract count datasets from phyloseq
count_Co2.4_all <- t(otu_table(Co2.4.phy))
#Transform in binary matrix
count_Co2.4_all[count_Co2.4_all> 0] <- 1
#Convert to dataframe
df_Co2.4_all <- as.data.frame(count_Co2.4_all)
#Read abundance per ASV
Co2.4_rel.abu <- data.frame(taxa_sums(Co2.4.phy))
Co2.4_rel.abu$log <- log10(Co2.4_rel.abu$taxa_sums.Co2.4.phy.)
#Merge ASV dataframe and ASV abundance
df_Co2.4_all_upset <-merge.data.frame(df_Co2.4_all, Co2.4_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co2.4_taxa <- tax_table(Co2.4.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co2.4_taxa)[names(df_Co2.4_taxa)=="rowname"] <- "Row.names"

df_Co2.4_all_upset2 <- full_join(df_Co2.4_all_upset, df_Co2.4_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co2.4_all <- UpSetR::upset(df_Co2.4_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE",
                                 text.scale = 2)
Upset_Co2.4_all # 9 ASVs out of 86 stigma, 25 seed

# Co2.6
Co2.6.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co2"&Plant_rep == "6") 
Co2.6.phy <- prune_taxa(taxa_sums(Co2.6.phy)>0, Co2.6.phy)

#Extract count datasets from phyloseq
count_Co2.6_all <- t(otu_table(Co2.6.phy))
#Transform in binary matrix
count_Co2.6_all[count_Co2.6_all> 0] <- 1
#Convert to dataframe
df_Co2.6_all <- as.data.frame(count_Co2.6_all)
#Read abundance per ASV
Co2.6_rel.abu <- data.frame(taxa_sums(Co2.6.phy))
Co2.6_rel.abu$log <- log10(Co2.6_rel.abu$taxa_sums.Co2.6.phy.)
#Merge ASV dataframe and ASV abundance
df_Co2.6_all_upset <-merge.data.frame(df_Co2.6_all, Co2.6_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co2.6_taxa <- tax_table(Co2.6.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co2.6_taxa)[names(df_Co2.6_taxa)=="rowname"] <- "Row.names"

df_Co2.6_all_upset2 <- full_join(df_Co2.6_all_upset, df_Co2.6_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co2.6_all <- UpSetR::upset(df_Co2.6_all_upset2,
                                 nsets = 2, 
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE",
                                 text.scale = 2)
Upset_Co2.6_all # 7 ASVs out of 51 stigma, 33 seed

# Co2.8: 39.08096
Co2.8.phy <- subset_samples(bact.rel.phy, Updated_Field_Name == "Co2"&Plant_rep == "8") 
Co2.8.phy <- prune_taxa(taxa_sums(Co2.8.phy)>0, Co2.8.phy)

#Extract count datasets from phyloseq
count_Co2.8_all <- t(otu_table(Co2.8.phy))
#Transform in binary matrix
count_Co2.8_all[count_Co2.8_all> 0] <- 1
#Convert to dataframe
df_Co2.8_all <- as.data.frame(count_Co2.8_all)
#Read abundance per ASV
Co2.8_rel.abu <- data.frame(taxa_sums(Co2.8.phy))
Co2.8_rel.abu$log <- log10(Co2.8_rel.abu$taxa_sums.Co2.8.phy.)
#Merge ASV dataframe and ASV abundance
df_Co2.8_all_upset <-merge.data.frame(df_Co2.8_all, Co2.8_rel.abu, by="row.names",
                                      all.x=TRUE)
#Append taxonomy to upset dataframe
df_Co2.8_taxa <- tax_table(Co2.8.phy) %>% data.frame() %>% rownames_to_column()
names(df_Co2.8_taxa)[names(df_Co2.8_taxa)=="rowname"] <- "Row.names"

df_Co2.8_all_upset2 <- full_join(df_Co2.8_all_upset, df_Co2.8_taxa, by="Row.names")

#Upset graph will all comparisons
Upset_Co2.8_all <- UpSetR::upset(df_Co2.8_all_upset2,
                                 nsets = 2,
                                 nintersects = NA,
                                 order.by = "freq",
                                 mainbar.y.label = "Number of ASVs",
                                 sets.x.label="total tissue ASV richness",
                                 decreasing ="TRUE", 
                                 text.scale = 2)
Upset_Co2.8_all # 5 ASVs out of 32 stigma, 54 seed

# Summarizing these comparisons into a plot - data tabulated in Excel based on ASV counts commented above and from viewing relative abundances of shared ASVs for each pair in ASV tables
  ## This will be updated to include code for data tabulation directly in R
pairwise.stats <- read.csv("data_output_seq/Field_survey/clean_pairwise_comparison_stats.csv")


seed.pair.rich <- pairwise.stats %>%
  ggplot(aes(x=Field, y=shared_seed_asv_percent)) +
  geom_boxplot() +
  geom_jitter(width = 0.25, size=2.5) +
  labs(x="Field",y="Proportion of richness\nfrom shared taxa(%)",
       color="Number of\nshared ASVs", fill = "Tissue type") +
  scale_fill_manual(labels = c("Stigmas", "Pooled seeds"), values = c("#F0E442", "#CC79A7")) +
  theme_light() +
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=14), 
        legend.position = "right", 
        legend.text = element_text(size=18))
seed.pair.rich

seed.pair.abund <- pairwise.stats %>%
  ggplot(aes(x=Field, y=shared_seed_asv_abund)) +
  geom_boxplot() +
  geom_jitter(width = 0.25, size=2.5) +
  labs(x="Field",y="Proportion of reads\nfrom shared taxa (%)",
       color="Number of\nshared ASVs", fill="Tissue type") +
  scale_fill_manual(labels = c("Stigmas", "Pooled seeds"), values = c("#F0E442", "#CC79A7")) +
  theme_light() +
  theme(axis.title = element_text(size=16), 
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=15),
        legend.position = "right")
seed.pair.abund

Fig2B <- ggarrange(seed.pair.rich, seed.pair.abund, common.legend = T)
Fig2B 

# making the same type of plots for stigmas
FigS2A <- pairwise.stats %>%
  ggplot(aes(x=Field, y=shared_stigma_asv_percent)) +
  geom_boxplot() +
  geom_jitter(width = 0.25, size=2.5) +
  labs(x="Field",y="Proportion of richness\nfrom shared taxa(%)") +
  theme_light() +
  theme(axis.title = element_text(size=18), 
        axis.text = element_text(size=15),
        legend.position = "bottom")
FigS2A

FigS2B <- pairwise.stats %>%
  ggplot(aes(x=Field, y=shared_stigma_asv_abund)) +
  geom_boxplot() +
  geom_jitter(width = 0.25, size=2.5) +
  labs(x="Field",y="Proportion of reads\nfrom shared taxa (%)",
       color="Number of\nshared ASVs") +
  theme_light() +
  theme(axis.title = element_text(size=18), 
        axis.title.y = element_text(size=16),
        axis.text = element_text(size=15), 
        legend.title = element_text(size=18), 
        legend.text = element_text(size=15),
        legend.position = "bottom")
FigS2B

#### Partitioning beta-diversity for individual sample pairs (Fig. S3) ####
#Extract count datasets from phyloseq
SJ1.1.count <- t(otu_table(SJ1.1.phy))
SJ1.2.count <- t(otu_table(SJ1.2.phy))
SJ1.3.count <- t(otu_table(SJ1.3.phy))
SJ1.4.count <- t(otu_table(SJ1.4.phy))
SJ1.5.count <- t(otu_table(SJ1.5.phy))
SJ1.6.count <- t(otu_table(SJ1.6.phy))
SJ1.7.count <- t(otu_table(SJ1.7.phy))
SJ1.8.count <- t(otu_table(SJ1.8.phy))
SJ2.2.count <- t(otu_table(SJ2.2.phy))
SJ2.3.count <- t(otu_table(SJ2.3.phy))
SJ2.5.count <- t(otu_table(SJ2.5.phy))
SJ2.6.count <- t(otu_table(SJ2.6.phy))
Co1.1.count <- t(otu_table(Co1.1.phy))
Co1.3.count <- t(otu_table(Co1.3.phy))
Co1.4.count <- t(otu_table(Co1.4.phy))
Co1.5.count <- t(otu_table(Co1.5.phy))
Co1.6.count <- t(otu_table(Co1.6.phy))
Co1.8.count <- t(otu_table(Co1.8.phy))
Co1.10.count <- t(otu_table(Co1.10.phy))
Co2.1.count <- t(otu_table(Co2.1.phy))
Co2.2.count <- t(otu_table(Co2.2.phy))
Co2.4.count <- t(otu_table(Co2.4.phy))
Co2.6.count <- t(otu_table(Co2.6.phy))
Co2.8.count <- t(otu_table(Co2.8.phy))
Co2.9.count <- t(otu_table(Co2.9.phy))

#Transform in binary matrix - skip for abundance-weighted metrics
SJ1.1.count[SJ1.1.count>0] <- 1
SJ1.2.count[SJ1.2.count>0] <- 1
SJ1.3.count[SJ1.3.count>0] <- 1
SJ1.4.count[SJ1.4.count>0] <- 1
SJ1.5.count[SJ1.5.count>0] <- 1
SJ1.6.count[SJ1.6.count>0] <- 1
SJ1.7.count[SJ1.7.count>0] <- 1
SJ1.8.count[SJ1.8.count>0] <- 1
SJ2.2.count[SJ2.2.count>0] <- 1
SJ2.3.count[SJ2.3.count>0] <- 1
SJ2.5.count[SJ2.5.count>0] <- 1
SJ2.6.count[SJ2.6.count>0] <- 1
Co1.1.count[Co1.1.count>0] <- 1
Co1.3.count[Co1.3.count>0] <- 1
Co1.4.count[Co1.4.count>0] <- 1
Co1.5.count[Co1.5.count>0] <- 1
Co1.6.count[Co1.6.count>0] <- 1
Co1.8.count[Co1.8.count>0] <- 1
Co1.10.count[Co1.10.count>0] <- 1
Co2.1.count[Co2.1.count>0] <- 1
Co2.2.count[Co2.2.count>0] <- 1
Co2.4.count[Co2.4.count>0] <- 1
Co2.6.count[Co2.6.count>0] <- 1
Co2.8.count[Co2.8.count>0] <- 1
Co2.9.count[Co2.9.count>0] <- 1

#Convert to dataframe
df_SJ1.1_count <- as.data.frame(t(SJ1.1.count))
df_SJ1.2_count <- as.data.frame(t(SJ1.2.count))
df_SJ1.3_count <- as.data.frame(t(SJ1.3.count))
df_SJ1.4_count <- as.data.frame(t(SJ1.4.count))
df_SJ1.5_count <- as.data.frame(t(SJ1.5.count))
df_SJ1.6_count <- as.data.frame(t(SJ1.6.count))
df_SJ1.7_count <- as.data.frame(t(SJ1.7.count))
df_SJ1.8_count <- as.data.frame(t(SJ1.8.count))
df_SJ2.2_count <- as.data.frame(t(SJ2.2.count))
df_SJ2.3_count <- as.data.frame(t(SJ2.3.count))
df_SJ2.5_count <- as.data.frame(t(SJ2.5.count))
df_SJ2.6_count <- as.data.frame(t(SJ2.6.count))
df_Co1.1_count <- as.data.frame(t(Co1.1.count))
df_Co1.3_count <- as.data.frame(t(Co1.3.count))
df_Co1.4_count <- as.data.frame(t(Co1.4.count))
df_Co1.5_count <- as.data.frame(t(Co1.5.count))
df_Co1.6_count <- as.data.frame(t(Co1.6.count))
df_Co1.8_count <- as.data.frame(t(Co1.8.count))
df_Co1.10_count <- as.data.frame(t(Co1.10.count))
df_Co2.1_count <- as.data.frame(t(Co2.1.count))
df_Co2.2_count <- as.data.frame(t(Co2.2.count))
df_Co2.4_count <- as.data.frame(t(Co2.4.count))
df_Co2.6_count <- as.data.frame(t(Co2.6.count))
df_Co2.8_count <- as.data.frame(t(Co2.8.count))
df_Co2.9_count <- as.data.frame(t(Co2.9.count))

# Betapart - generate phylogenetic trees from the phyloseq objects 
seqs.SJ1.1 <- refseq(SJ1.1.phy)
alignment.SJ1.1 <- AlignTranslation(seqs.SJ1.1, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ1.1 <- phyDat(as(alignment.SJ1.1, "matrix"), type="DNA")
dm.SJ1.1 <- dist.ml(phang.align.SJ1.1)
treeNJ.SJ1.1 <- NJ(dm.SJ1.1)
phylo.SJ1.1 <- root(treeNJ.SJ1.1, outgroup = 1, resolve.root = T)

seqs.SJ1.2 <- refseq(SJ1.2.phy)
alignment.SJ1.2 <- AlignTranslation(seqs.SJ1.2, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ1.2 <- phyDat(as(alignment.SJ1.2, "matrix"), type="DNA")
dm.SJ1.2 <- dist.ml(phang.align.SJ1.2)
treeNJ.SJ1.2 <- NJ(dm.SJ1.2)
phylo.SJ1.2 <- root(treeNJ.SJ1.2, outgroup = 1, resolve.root = T)

seqs.SJ1.3 <- refseq(SJ1.3.phy)
alignment.SJ1.3 <- AlignTranslation(seqs.SJ1.3, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ1.3 <- phyDat(as(alignment.SJ1.3, "matrix"), type="DNA")
dm.SJ1.3 <- dist.ml(phang.align.SJ1.3)
treeNJ.SJ1.3 <- NJ(dm.SJ1.3)
phylo.SJ1.3 <- root(treeNJ.SJ1.3, outgroup = 1, resolve.root = T)

seqs.SJ1.4 <- refseq(SJ1.4.phy)
alignment.SJ1.4 <- AlignTranslation(seqs.SJ1.4, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ1.4 <- phyDat(as(alignment.SJ1.4, "matrix"), type="DNA")
dm.SJ1.4 <- dist.ml(phang.align.SJ1.4)
treeNJ.SJ1.4 <- NJ(dm.SJ1.4)
phylo.SJ1.4 <- root(treeNJ.SJ1.4, outgroup = 1, resolve.root = T)

seqs.SJ1.5 <- refseq(SJ1.5.phy)
alignment.SJ1.5 <- AlignTranslation(seqs.SJ1.5, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ1.5 <- phyDat(as(alignment.SJ1.5, "matrix"), type="DNA")
dm.SJ1.5 <- dist.ml(phang.align.SJ1.5)
treeNJ.SJ1.5 <- NJ(dm.SJ1.5)
phylo.SJ1.5 <- root(treeNJ.SJ1.5, outgroup = 1, resolve.root = T)

seqs.SJ1.6 <- refseq(SJ1.6.phy)
alignment.SJ1.6 <- AlignTranslation(seqs.SJ1.6, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ1.6 <- phyDat(as(alignment.SJ1.6, "matrix"), type="DNA")
dm.SJ1.6 <- dist.ml(phang.align.SJ1.6)
treeNJ.SJ1.6 <- NJ(dm.SJ1.6)
phylo.SJ1.6 <- root(treeNJ.SJ1.6, outgroup = 1, resolve.root = T)

seqs.SJ1.7 <- refseq(SJ1.7.phy)
alignment.SJ1.7 <- AlignTranslation(seqs.SJ1.7, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ1.7 <- phyDat(as(alignment.SJ1.7, "matrix"), type="DNA")
dm.SJ1.7 <- dist.ml(phang.align.SJ1.7)
treeNJ.SJ1.7 <- NJ(dm.SJ1.7)
phylo.SJ1.7 <- root(treeNJ.SJ1.7, outgroup = 1, resolve.root = T)

seqs.SJ1.8 <- refseq(SJ1.8.phy)
alignment.SJ1.8 <- AlignTranslation(seqs.SJ1.8, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ1.8 <- phyDat(as(alignment.SJ1.8, "matrix"), type="DNA")
dm.SJ1.8 <- dist.ml(phang.align.SJ1.8)
treeNJ.SJ1.8 <- NJ(dm.SJ1.8)
phylo.SJ1.8 <- root(treeNJ.SJ1.8, outgroup = 1, resolve.root = T)

seqs.SJ2.2 <- refseq(SJ2.2.phy)
alignment.SJ2.2 <- AlignTranslation(seqs.SJ2.2, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ2.2 <- phyDat(as(alignment.SJ2.2, "matrix"), type="DNA")
dm.SJ2.2 <- dist.ml(phang.align.SJ2.2)
treeNJ.SJ2.2 <- NJ(dm.SJ2.2)
phylo.SJ2.2 <- root(treeNJ.SJ2.2, outgroup = 1, resolve.root = T)

seqs.SJ2.3 <- refseq(SJ2.3.phy)
alignment.SJ2.3 <- AlignTranslation(seqs.SJ2.3, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ2.3 <- phyDat(as(alignment.SJ2.3, "matrix"), type="DNA")
dm.SJ2.3 <- dist.ml(phang.align.SJ2.3)
treeNJ.SJ2.3 <- NJ(dm.SJ2.3)
phylo.SJ2.3 <- root(treeNJ.SJ2.3, outgroup = 1, resolve.root = T)

seqs.SJ2.5 <- refseq(SJ2.5.phy)
alignment.SJ2.5 <- AlignTranslation(seqs.SJ2.5, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ2.5 <- phyDat(as(alignment.SJ2.5, "matrix"), type="DNA")
dm.SJ2.5 <- dist.ml(phang.align.SJ2.5)
treeNJ.SJ2.5 <- NJ(dm.SJ2.5)
phylo.SJ2.5 <- root(treeNJ.SJ2.5, outgroup = 1, resolve.root = T)

seqs.SJ2.6 <- refseq(SJ2.6.phy)
alignment.SJ2.6 <- AlignTranslation(seqs.SJ2.6, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.SJ2.6 <- phyDat(as(alignment.SJ2.6, "matrix"), type="DNA")
dm.SJ2.6 <- dist.ml(phang.align.SJ2.6)
treeNJ.SJ2.6 <- NJ(dm.SJ2.6)
phylo.SJ2.6 <- root(treeNJ.SJ2.6, outgroup = 1, resolve.root = T)

seqs.Co1.1 <- refseq(Co1.1.phy)
alignment.Co1.1 <- AlignTranslation(seqs.Co1.1, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co1.1 <- phyDat(as(alignment.Co1.1, "matrix"), type="DNA")
dm.Co1.1 <- dist.ml(phang.align.Co1.1)
treeNJ.Co1.1 <- NJ(dm.Co1.1)
phylo.Co1.1 <- root(treeNJ.Co1.1, outgroup = 1, resolve.root = T)

seqs.Co1.3 <- refseq(Co1.3.phy) 
alignment.Co1.3 <- AlignTranslation(seqs.Co1.3, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co1.3 <- phyDat(as(alignment.Co1.3, "matrix"), type="DNA")
dm.Co1.3 <- dist.ml(phang.align.Co1.3)
treeNJ.Co1.3 <- NJ(dm.Co1.3)
phylo.Co1.3 <- root(treeNJ.Co1.3, outgroup = 1, resolve.root = T)

seqs.Co1.4 <- refseq(Co1.4.phy)
alignment.Co1.4 <- AlignTranslation(seqs.Co1.4, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co1.4 <- phyDat(as(alignment.Co1.4, "matrix"), type="DNA")
dm.Co1.4 <- dist.ml(phang.align.Co1.4)
treeNJ.Co1.4 <- NJ(dm.Co1.4)
phylo.Co1.4 <- root(treeNJ.Co1.4, outgroup = 1, resolve.root = T)

seqs.Co1.5 <- refseq(Co1.5.phy)
alignment.Co1.5 <- AlignTranslation(seqs.Co1.5, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co1.5 <- phyDat(as(alignment.Co1.5, "matrix"), type="DNA")
dm.Co1.5 <- dist.ml(phang.align.Co1.5)
treeNJ.Co1.5 <- NJ(dm.Co1.5)
phylo.Co1.5 <- root(treeNJ.Co1.5, outgroup = 1, resolve.root = T)

seqs.Co1.6 <- refseq(Co1.6.phy)
alignment.Co1.6 <- AlignTranslation(seqs.Co1.6, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co1.6 <- phyDat(as(alignment.Co1.6, "matrix"), type="DNA")
dm.Co1.6 <- dist.ml(phang.align.Co1.6)
treeNJ.Co1.6 <- NJ(dm.Co1.6)
phylo.Co1.6 <- root(treeNJ.Co1.6, outgroup = 1, resolve.root = T)

seqs.Co1.8 <- refseq(Co1.8.phy)
alignment.Co1.8 <- AlignTranslation(seqs.Co1.8, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co1.8 <- phyDat(as(alignment.Co1.8, "matrix"), type="DNA")
dm.Co1.8 <- dist.ml(phang.align.Co1.8)
treeNJ.Co1.8 <- NJ(dm.Co1.8)
phylo.Co1.8 <- root(treeNJ.Co1.8, outgroup = 1, resolve.root = T)

seqs.Co1.10 <- refseq(Co1.10.phy)
alignment.Co1.10 <- AlignTranslation(seqs.Co1.10, sense = "+", readingFrame = 2, 
                                     type ="DNAStringSet") 
phang.align.Co1.10 <- phyDat(as(alignment.Co1.10, "matrix"), type="DNA")
dm.Co1.10 <- dist.ml(phang.align.Co1.10)
treeNJ.Co1.10 <- NJ(dm.Co1.10)
phylo.Co1.10 <- root(treeNJ.Co1.10, outgroup = 1, resolve.root = T)

seqs.Co2.1 <- refseq(Co2.1.phy)
alignment.Co2.1 <- AlignTranslation(seqs.Co2.1, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co2.1 <- phyDat(as(alignment.Co2.1, "matrix"), type="DNA")
dm.Co2.1 <- dist.ml(phang.align.Co2.1)
treeNJ.Co2.1 <- NJ(dm.Co2.1)
phylo.Co2.1 <- root(treeNJ.Co2.1, outgroup = 1, resolve.root = T)

seqs.Co2.2 <- refseq(Co2.2.phy)
alignment.Co2.2 <- AlignTranslation(seqs.Co2.2, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co2.2 <- phyDat(as(alignment.Co2.2, "matrix"), type="DNA")
dm.Co2.2 <- dist.ml(phang.align.Co2.2)
treeNJ.Co2.2 <- NJ(dm.Co2.2)
phylo.Co2.2 <- root(treeNJ.Co2.2, outgroup = 1, resolve.root = T)

seqs.Co2.4 <- refseq(Co2.4.phy)
alignment.Co2.4 <- AlignTranslation(seqs.Co2.4, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co2.4 <- phyDat(as(alignment.Co2.4, "matrix"), type="DNA")
dm.Co2.4 <- dist.ml(phang.align.Co2.4)
treeNJ.Co2.4 <- NJ(dm.Co2.4)
phylo.Co2.4 <- root(treeNJ.Co2.4, outgroup = 1, resolve.root = T)

seqs.Co2.6 <- refseq(Co2.6.phy)
alignment.Co2.6 <- AlignTranslation(seqs.Co2.6, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co2.6 <- phyDat(as(alignment.Co2.6, "matrix"), type="DNA")
dm.Co2.6 <- dist.ml(phang.align.Co2.6)
treeNJ.Co2.6 <- NJ(dm.Co2.6)
phylo.Co2.6 <- root(treeNJ.Co2.6, outgroup = 1, resolve.root = T)

seqs.Co2.8 <- refseq(Co2.8.phy)
alignment.Co2.8 <- AlignTranslation(seqs.Co2.8, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co2.8 <- phyDat(as(alignment.Co2.8, "matrix"), type="DNA")
dm.Co2.8 <- dist.ml(phang.align.Co2.8)
treeNJ.Co2.8 <- NJ(dm.Co2.8)
phylo.Co2.8 <- root(treeNJ.Co2.8, outgroup = 1, resolve.root = T)

seqs.Co2.9 <- refseq(Co2.9.phy)
alignment.Co2.9 <- AlignTranslation(seqs.Co2.9, sense = "+", readingFrame = 2, type ="DNAStringSet") 
phang.align.Co2.9 <- phyDat(as(alignment.Co2.9, "matrix"), type="DNA")
dm.Co2.9 <- dist.ml(phang.align.Co2.9)
treeNJ.Co2.9 <- NJ(dm.Co2.9)
phylo.Co2.9 <- root(treeNJ.Co2.9, outgroup = 1, resolve.root = T)

# Betapart - run the analysis 
SJ1.1.betapart <- phylo.beta.pair(df_SJ1.1_count, phylo.SJ1.1, index.family = "jaccard")
SJ1.2.betapart <- phylo.beta.pair(df_SJ1.2_count, phylo.SJ1.2, index.family = "jaccard")
SJ1.3.betapart <- phylo.beta.pair(df_SJ1.3_count, phylo.SJ1.3, index.family = "jaccard")
SJ1.4.betapart <- phylo.beta.pair(df_SJ1.4_count, phylo.SJ1.4, index.family = "jaccard")
SJ1.5.betapart <- phylo.beta.pair(df_SJ1.5_count, phylo.SJ1.5, index.family = "jaccard")
SJ1.6.betapart <- phylo.beta.pair(df_SJ1.6_count, phylo.SJ1.6, index.family = "jaccard")
SJ1.7.betapart <- phylo.beta.pair(df_SJ1.7_count, phylo.SJ1.7, index.family = "jaccard")
SJ1.8.betapart <- phylo.beta.pair(df_SJ1.8_count, phylo.SJ1.8, index.family = "jaccard")
SJ2.2.betapart <- phylo.beta.pair(df_SJ2.2_count, phylo.SJ2.2, index.family = "jaccard")
SJ2.3.betapart <- phylo.beta.pair(df_SJ2.3_count, phylo.SJ2.3, index.family = "jaccard")
SJ2.5.betapart <- phylo.beta.pair(df_SJ2.5_count, phylo.SJ2.5, index.family = "jaccard")
SJ2.6.betapart <- phylo.beta.pair(df_SJ2.6_count, phylo.SJ2.6, index.family = "jaccard")
Co1.1.betapart <- phylo.beta.pair(df_Co1.1_count, phylo.Co1.1, index.family = "jaccard")
Co1.3.betapart <- phylo.beta.pair(df_Co1.3_count, phylo.Co1.3, index.family = "jaccard")
Co1.4.betapart <- phylo.beta.pair(df_Co1.4_count, phylo.Co1.4, index.family = "jaccard")
Co1.5.betapart <- phylo.beta.pair(df_Co1.5_count, phylo.Co1.5, index.family = "jaccard")
Co1.6.betapart <- phylo.beta.pair(df_Co1.6_count, phylo.Co1.6, index.family = "jaccard")
Co1.8.betapart <- phylo.beta.pair(df_Co1.8_count, phylo.Co1.8, index.family = "jaccard")
Co1.10.betapart <- phylo.beta.pair(df_Co1.10_count, phylo.Co1.10, index.family = "jaccard")
Co2.1.betapart <- phylo.beta.pair(df_Co2.1_count, phylo.Co2.1, index.family = "jaccard")
Co2.2.betapart <- phylo.beta.pair(df_Co2.2_count, phylo.Co2.2, index.family = "jaccard")
Co2.4.betapart <- phylo.beta.pair(df_Co2.4_count, phylo.Co2.4, index.family = "jaccard")
Co2.6.betapart <- phylo.beta.pair(df_Co2.6_count, phylo.Co2.6, index.family = "jaccard")
Co2.8.betapart <- phylo.beta.pair(df_Co2.8_count, phylo.Co2.8, index.family = "jaccard")
Co2.9.betapart <- phylo.beta.pair(df_Co2.9_count, phylo.Co2.9, index.family = "jaccard")

# visualizing the decomposed dispersion
SJ1.1.betapart.df <- as.data.frame(SJ1.1.betapart,row.names = "SJ1.1") %>% rownames_to_column(var = "SampID")
SJ1.2.betapart.df <- as.data.frame(SJ1.2.betapart, row.names = "SJ1.2") %>% 
  rownames_to_column(var = "SampID")
SJ1.3.betapart.df <- as.data.frame(SJ1.3.betapart, row.names = "SJ1.3") %>% 
  rownames_to_column(var = "SampID")
SJ1.4.betapart.df <- as.data.frame(SJ1.4.betapart, row.names = "SJ1.4") %>% 
  rownames_to_column(var = "SampID")
SJ1.5.betapart.df <- as.data.frame(SJ1.5.betapart, row.names = "SJ1.5") %>% 
  rownames_to_column(var = "SampID")
SJ1.6.betapart.df <- as.data.frame(SJ1.6.betapart, row.names = "SJ1.6") %>% 
  rownames_to_column(var = "SampID")
SJ1.7.betapart.df <- as.data.frame(SJ1.7.betapart, row.names = "SJ1.7") %>% 
  rownames_to_column(var = "SampID")
SJ1.8.betapart.df <- as.data.frame(SJ1.8.betapart, row.names = "SJ1.8") %>% 
  rownames_to_column(var = "SampID")
SJ2.2.betapart.df <- as.data.frame(SJ2.2.betapart, row.names = "SJ2.2") %>% 
  rownames_to_column(var = "SampID")
SJ2.3.betapart.df <- as.data.frame(SJ2.3.betapart, row.names = "SJ2.3") %>% 
  rownames_to_column(var = "SampID")
SJ2.5.betapart.df <- as.data.frame(SJ2.5.betapart, row.names = "SJ2.5") %>% 
  rownames_to_column(var = "SampID")
SJ2.6.betapart.df <- as.data.frame(SJ2.6.betapart, row.names = "SJ2.6") %>% 
  rownames_to_column(var = "SampID")
Co1.1.betapart.df <- as.data.frame(Co1.1.betapart, row.names = "Co1.1") %>% 
  rownames_to_column(var = "SampID")
Co1.3.betapart.df <- as.data.frame(Co1.3.betapart, row.names = "Co1.3") %>% 
  rownames_to_column(var = "SampID")
Co1.4.betapart.df <- as.data.frame(Co1.4.betapart, row.names = "Co1.4") %>% 
  rownames_to_column(var = "SampID")
Co1.5.betapart.df <- as.data.frame(Co1.5.betapart, row.names = "Co1.5") %>% 
  rownames_to_column(var = "SampID")
Co1.6.betapart.df <- as.data.frame(Co1.6.betapart, row.names = "Co1.6") %>% 
  rownames_to_column(var = "SampID")
Co1.8.betapart.df <- as.data.frame(Co1.8.betapart, row.names = "Co1.8") %>% 
  rownames_to_column(var = "SampID")
Co1.10.betapart.df <- as.data.frame(Co1.10.betapart, row.names = "Co1.10") %>%
  rownames_to_column(var = "SampID")
Co2.1.betapart.df <- as.data.frame(Co2.1.betapart, row.names = "Co2.1") %>%
  rownames_to_column(var = "SampID")
Co2.2.betapart.df <- as.data.frame(Co2.2.betapart, row.names = "Co2.2") %>%
  rownames_to_column(var = "SampID")
Co2.4.betapart.df <- as.data.frame(Co2.4.betapart, row.names = "Co2.4") %>%
  rownames_to_column(var = "SampID")
Co2.6.betapart.df <- as.data.frame(Co2.6.betapart, row.names = "Co2.6") %>%
  rownames_to_column(var = "SampID")
Co2.8.betapart.df <- as.data.frame(Co2.8.betapart, row.names = "Co2.8") %>%
  rownames_to_column(var = "SampID")
Co2.9.betapart.df <- as.data.frame(Co2.9.betapart, row.names = "Co2.9") %>%
  rownames_to_column(var = "SampID")

all.pairs.beta <- SJ1.1.betapart.df %>%
  bind_rows(SJ1.2.betapart.df) %>%
  bind_rows(SJ1.3.betapart.df) %>%
  bind_rows(SJ1.4.betapart.df) %>%
  bind_rows(SJ1.5.betapart.df) %>%
  bind_rows(SJ1.6.betapart.df) %>%
  bind_rows(SJ1.7.betapart.df) %>%
  bind_rows(SJ1.8.betapart.df) %>%
  bind_rows(SJ2.2.betapart.df) %>%
  bind_rows(SJ2.3.betapart.df) %>%
  bind_rows(SJ2.5.betapart.df) %>%
  bind_rows(SJ2.6.betapart.df) %>%
  bind_rows(Co1.1.betapart.df) %>%
  bind_rows(Co1.10.betapart.df) %>%
  bind_rows(Co1.3.betapart.df) %>%
  bind_rows(Co1.4.betapart.df) %>%
  bind_rows(Co1.5.betapart.df) %>%
  bind_rows(Co1.6.betapart.df) %>%
  bind_rows(Co1.8.betapart.df) %>%
  bind_rows(Co2.1.betapart.df) %>%
  bind_rows(Co2.2.betapart.df) %>%
  bind_rows(Co2.4.betapart.df) %>%
  bind_rows(Co2.6.betapart.df) %>%
  bind_rows(Co2.8.betapart.df)

all.pairs.beta <- all.pairs.beta %>% 
  mutate(jtu.percent = (phylo.beta.jtu/phylo.beta.jac)*100) %>%
  mutate(jne.percent = (phylo.beta.jne/phylo.beta.jac)*100)

all.pairs.beta$Field <- c("SJ1", "SJ1", "SJ1", "SJ1","SJ1","SJ1","SJ1","SJ1","SJ2","SJ2","SJ2","SJ2",
                          "Co1","Co1","Co1","Co1","Co1","Co1","Co1","Co2","Co2","Co2","Co2","Co2")

pairs.beta.long <- all.pairs.beta %>%
  select(SampID, Field, jtu.percent, jne.percent) %>%
  pivot_longer(cols = c(jtu.percent, jne.percent), names_to = "beta_type", values_to = "percent")

FigS3A <- pairs.beta.long %>%
  ggplot(aes(x=SampID, y=percent, fill=beta_type)) +
  geom_col(position = "stack") +
  labs(x="Pair",y="Percent of beta-diversity (Jaccard)", fill="Diversity component")+
  facet_wrap(~Field, scales = "free_x") + 
  scale_fill_discrete(name = "Diversity\ncomponent", labels = c("Nestedness", "Turnover")) +
  theme_light() +
  theme(axis.title = element_text(size=18),
        axis.text.y = element_text(size=16), 
        axis.text.x = element_text(size=14),
        legend.title = element_text(size=18), 
        legend.text = element_text(size=16))
FigS3A 

#### visualizing the abundant taxa in pairwise comparisons (Fig. S3B) ####


#### visualizing most abundant shared + unique genera (Fig. 3A, Fig. S4) ####
# shared taxa
bact.rel.taxa <- tax_table(bact.rel.phy) %>% data.frame() %>% rownames_to_column()
names(bact.rel.taxa)[names(bact.rel.taxa) == "rowname"] <- "Row.names"
bact.rel.taxa2 <- merge(bact.rel.taxa, df_16S_all_upset3, by="Row.names")
bact.rel.taxa2 <- column_to_rownames(bact.rel.taxa2, var = "Row.names")
tax_table(bact.rel.phy) <- as.matrix(bact.rel.taxa2)
shared.rel.phy <- subset_taxa(bact.rel.phy, stigma=="TRUE"&seed_pooled=="TRUE")
shared.rel.phy <- prune_samples(sample_sums(shared.rel.phy)>0, shared.rel.phy)
sample_data(shared.rel.phy)$tissue_type <- factor(sample_data(shared.rel.phy)$tissue_type,
                                                  levels = c("stigma", "seed_pooled"))
merged.share.phy <- shared.rel.phy %>%
  aggregate_top_taxa2(., 20, "genus.x")

scales::show_col(brewerPlus)

Fig3A <- plot_composition(merged.share.phy, 
                          verbose = F, otu.sort = "abundance", group_by = "tissue_type") +
  scale_fill_manual(values = brewerPlus) +
  labs(x="Sample", y="Relative abundance (%)", fill="Genus") +
  #facet_grid(~tissue_type, scales = "free_x", space = "free_x") +
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=16), 
        axis.title = element_text(size=18),
        legend.text = element_text(face = "italic", size=16), 
        legend.title = element_text(size=18),
        legend.position = "right", 
        strip.text = element_text(size=20)) +
  guides(fill=guide_legend(ncol = 1)) 
Fig3A 

# taxa unique to stigmas (Fig. S4A)
stigma.only.phy <- subset_taxa(bact.rel.phy, stigma=="TRUE"&seed_pooled=="FALSE")
stigma.only.phy <- prune_samples(sample_sums(stigma.only.phy)>0, stigma.only.phy)

merged.stigma.phy <- stigma.only.phy %>%
  aggregate_top_taxa2(., 20, "genus.x")
merged.stigma.taxa <- tax_table(merged.stigma.phy) %>% as.data.frame()

# need to make a manual color-scheme where colors are shared for overlapping genera, but others are unique
kelly <- distinct_palette(pal = "kelly")
scales::show_col(kelly)

FigS4A <- plot_composition(merged.stigma.phy, 
                           verbose = F, otu.sort = "abundance", group_by = "tissue_type") +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#f3c300", "#875692", "#f38400","#be0032",
                               "#c2b280","#848482","#008856","#e68fac", "#1F78B4","#f99379",
                               "#604e97","#f6a600","#b3446c","#666666","#dcd300","#882d17",
                               "#8db600","#654522","#e25822","#2b3d26")) + # Corynebacterium is 16
  labs(x="Sample", y="Relative abundance (%)", fill="Genus") +
  ylim(0,100)+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=16), 
        axis.title = element_text(size=18),
        legend.text = element_text(face = "italic", size=16), 
        legend.title = element_text(size=18),
        legend.position = "right", 
        strip.text = element_text(size=20)) +
  guides(fill=guide_legend(ncol = 1)) 
FigS4A

# taxa that are unique to seeds (Fig. S4B)
greenArmytage <- distinct_palette(pal = "greenArmytage")
greenArmytage

seed.only.phy <- subset_taxa(bact.rel.phy, stigma=="FALSE"&seed_pooled=="TRUE")
seed.only.phy <- prune_samples(sample_sums(seed.only.phy)>0, seed.only.phy)

merged.seed.phy <- seed.only.phy %>%
  aggregate_top_taxa2(., 20, "genus.x")
merged.seed.taxa <- tax_table(merged.seed.phy) %>% as.data.frame()
# shared with seed-only (colors): unknown ("#FB9A99"), other ("#A6CEE3"), prevotella ("#FFFF99"), corynebacterium ("#666666") 

FigS4B <- plot_composition(merged.seed.phy, 
                           verbose = F, otu.sort = "abundance", group_by = "tissue_type") +
  scale_fill_manual(values = c("#FB9A99", "#0075DC","#A6CEE3","#6A3D9A","#FFFF99","#993F00",
                               "#4C005C","#FFCC99","#5EF1F2","#666666","#94FFB5","#8F7C00",
                               "#9DCC00","#C20088","#003380","#FFA8BB","#FF0010","#00998F",
                               "#E0FF66","#100AFF", "#990000")) + 
  labs(x="Sample", y="Relative abundance (%)", fill="Genus") +
  ylim(0,100)+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=16), 
        axis.title = element_text(size=18),
        legend.text = element_text(face = "italic", size=14), 
        legend.title = element_text(size=16),
        legend.position = "right", 
        strip.text = element_text(size=20)) +
  guides(fill=guide_legend(ncol = 1, override.aes = list(size=4))) 
FigS4B

#### Differential abundance analysis with ANCOM-BC2 (Fig. 3B) ####
output <- ancombc2(data = shared.rel.phy, assay_name = "counts", tax_level = "genus.x",
                   fix_formula = "tissue_type", rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = TRUE,
                   prv_cut = 0.01, lib_cut = 0, s0_perc = 0.05,
                   group = "tissue_type", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                   iter_control = list(tol = 1e-2, max_iter = 20, 
                                       verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl(),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                   trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(-1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE)),
                                        node = list(2, 2, 1),
                                        solver = "ECOS",
                                        B = 100)) # could add field to the fix_formula argument

# primary analysis (differential abundance in seeds)
tab_zero = output$zero_ind
print(head(tab_zero))
res_prim = output$res # differences for stigmas are the (intercept), since it's the reference group

df_seed = res_prim %>%
  dplyr::select(taxon, ends_with("seed_pooled")) 
df_fig_seed = df_seed %>%
  dplyr::filter(diff_tissue_typeseed_pooled == 1) %>% 
  dplyr::arrange(desc(lfc_tissue_typeseed_pooled)) %>%
  dplyr::mutate(direct = ifelse(lfc_tissue_typeseed_pooled > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(passed_ss_tissue_typeseed_pooled == 1, "aquamarine3", "black"))

df_fig_seed$taxon = factor(df_fig_seed$taxon, levels = df_fig_seed$taxon)
df_fig_seed$direct = factor(df_fig_seed$direct, 
                            levels = c("Positive LFC", "Negative LFC"))
# filter out anything that failed the pseudo-count sensitivity test
df_fig_seed2 <- df_fig_seed %>%
  filter(color == "aquamarine3") 
#df_fig_seed2$genus <- c("Shigella", "Kosakonia", "Acidovorax", "Bradyrhizobium")
df_fig_seed2$genus <- factor(df_fig_seed2$genus, 
                             levels = c("Shigella", "Kosakonia", "Acidovorax",
                                        "Bradyrhizobium"))

Fig3B <- df_fig_seed2 %>%
  ggplot(aes(x = taxon, y = lfc_tissue_typeseed_pooled)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_tissue_typeseed_pooled - se_tissue_typeseed_pooled, 
                    ymax = lfc_tissue_typeseed_pooled + se_tissue_typeseed_pooled), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Genus", y = "Log fold change", 
       title = "Differential abundance in seeds") + 
  theme_light() + 
  theme(plot.title = element_text(size=20,hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=15,angle = 60, hjust = 1))
Fig3B
