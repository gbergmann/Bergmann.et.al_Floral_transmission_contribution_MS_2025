# Analysis of seedling viability experiment 

#### loading packages and data #### 
library(tidyverse)
library(cowplot)
library(viridis)
library(FSA)
library(ggsignif)
library(chisq.posthoc.test)
library(ggpubr)
library(rstatix)
library(RColorBrewer)

seed.data <- read.csv("data_raw/Seed_viability/GEB_seed_viability_assay_contd_F2023_bact.only_Field.MS.csv")

#### visualization and statistics on germination ####
seed.data$trt_genus <- factor(seed.data$trt_genus,
                              levels = c("control", "Acidovorax_citrulli","Erwinia","Paraburkholderia",
                                         "Rosenbergiella"))

seed.data <- seed.data %>%
  group_by(trt_strain) %>%
  mutate(germ_avg=mean(germ_yes_percent), pheno_avg=mean(normal_pheno_percent))

germ.plot <- seed.data %>%
  ggplot(aes(x=reorder(trt_strain, germ_yes_percent), y=germ_yes_percent, color=trt_genus)) +
  geom_point(size=3)+
  geom_errorbar(aes(ymin=germ_avg, ymax=germ_avg), linewidth=2) +
  scale_color_brewer(palette = "Paired")+
  coord_flip() +
  geom_hline(yintercept = 98.66667, linetype="dashed")+
  # annotate("text", x=20, y=0.9, label="Average rate\nof sterile control")+
  labs(x="Inoculated strain",y="Seed germination rate (%)", color="Genus") +
  theme_light() +
  theme(axis.title = element_text(size=15), axis.text = element_text(size=13),
        legend.title = element_text(size=15), legend.text = element_text(size=13))
germ.plot 

# calculating a contingency table and running chi-squared test of germination data
germ.tab <- seed.data %>%
  group_by(trt_strain) %>%
  summarise(germ_avg = mean(germ_yes_percent)) %>%
  mutate(no_germ_avg = 100-germ_avg) %>%
  column_to_rownames(var = "trt_strain")

chisq.test(germ.tab)
chisq.posthoc.test(germ.tab, method = "bonferroni")

#### visualization and statistics on seedling phenotypes ####
pheno.plot <- seed.data %>%
  ggplot(aes(x=reorder(trt_strain, normal_pheno_percent), 
             y=normal_pheno_percent, color=trt_genus)) +
  geom_point(size=3)+
  geom_errorbar(aes(ymin=pheno_avg, ymax=pheno_avg),linewidth=2) +
  scale_color_brewer(palette = "Paired")+
  coord_flip() +
  geom_hline(yintercept = 64.76449, linetype="dashed")+
  labs(x="Inoculated strain",y="Normal seedling rate (%)", color="Genus") +
  theme_light() +
  theme(axis.title = element_text(size=15), axis.text = element_text(size=13),
        legend.title = element_text(size=15), legend.text = element_text(size=13))
pheno.plot 

# calculating a contingency table and running chi-squared test of phenotype data
pheno.tab <- seed.data %>%
  group_by(trt_strain) %>%
  summarise(pheno_avg = mean(normal_pheno_percent)) %>%
  mutate(no_germ_avg = 100-pheno_avg) %>%
  column_to_rownames(var = "trt_strain")

chisq.test(pheno.tab)
chisq.posthoc.test(pheno.tab, method = "bonferroni")

#### combining plots for Fig 5B ####

Fig5B <- ggarrange(germ.plot, pheno.plot, nrow = 1, common.legend = T)
Fig5B

# I had difficulty labeling strains in bold when they were used in both the stigma inoculation experiment and this experiment. As such, I manually re-labeled the strains with bold/normal text in PowerPoint when formatting the figures