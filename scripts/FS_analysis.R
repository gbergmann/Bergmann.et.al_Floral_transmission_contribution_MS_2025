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

#### Partitioning beta-diversity for individual sample pairs (Fig. S3)  ####

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
