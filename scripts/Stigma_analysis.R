# Analysis of stigma inoculation experiment

#### loading packages and data ####
library(tidyverse)
library(cowplot)
library(pheatmap)
library(reshape2)
library(plotly)
library(ggpubr)
library(chisq.posthoc.test)
library(viridis)

# data files
slice.turbidity <- read.csv("data_raw/Stigma_inoc/single_strain_reisolation_FieldMS.csv")
fruit.abortion <- read.csv("data_raw/Stigma_inoc/single.strain_fruit.abortion_fieldMS.csv")

#### Calculating average re-isolation frequency per slice, making top part of Fig. 5A ####
slice.turbidity$slice <- factor(slice.turbidity$slice)

# calculating average values for the heatmap
slice.avg <- slice.turbidity %>%
  group_by(Microbe, slice) %>%
  summarise(freq = mean(turbid_percent))

slice.avg$Microbe <- factor(slice.avg$Microbe,
                            levels = c("PBS_control", "AcAAC001R", "AcM6R", "Bacillus", "Erwinia",
                                       "Microbacterium","Paraburkholderia","Rosenbergiella"))

# making the heatmap
slice.heat <- slice.avg %>%
  ggplot(aes(x=Microbe, y=slice, fill=freq)) +
  geom_tile(stat = "identity", color = "black", lwd = 0.5, linetype = 1) +
  geom_text(aes(label=freq), size=5) +
  geom_hline(yintercept = 3.5, color="black", linetype="dashed", linewidth=2) +
  scale_fill_gradient(high = "gray45", low = "white", limits= c(0,50)) +
  labs(x="Inoculum strain", y="Slice number",
       fill="Re-isolation\nfrequency (%)")+
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=15), 
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        axis.title.x = element_blank(),
        legend.position = "top", 
        text = element_text(size=15)) 
slice.heat

#### Calculating fruit abortion frequency per treatment, making bottom part of Fig. 5A ####
strain.abortion2 <- FS.strain.abortion %>%
  mutate(Bacterial_inoc. = recode(Bacterial_inoc., "AcAAC001R" = "A.citrulli_AAC001R",
                                  "AcM6R" = "A.citrulli_M6R",
                                  "Bacillus" = 'Bacillus_T.m10.b1',
                                  "Erwinia" = 'Erwinia_P.s8.b2', 
                                  "Microbacterium"=  'Microbacterium_T.s9.b4',
                                  "Paraburkholderia" = "Paraburkholderia_T.s8.b1",
                                  "Rosenbergiella" = "Rosenbergiella_P.s10.b1"))

strain.abortion2$Bacterial_inoc. <- factor(strain.abortion2$Bacterial_inoc.,
                                           levels = c("PBS_control","A.citrulli_AAC001R",
                                                      "A.citrulli_M6R","Bacillus_T.m10.b1",
                                                      "Erwinia_P.s8.b2","Microbacterium_T.s9.b4",
                                                      "Paraburkholderia_T.s8.b1",
                                                      "Rosenbergiella_P.s10.b1"))

fruit.bar <- fruit.abortion %>%
  ggplot(aes(x=Bacterial_inoc., y=success_rate, fill = Bacterial_inoc.)) +
  geom_col() +
  labs(x="Inoculum strain", y="Proportion of\ndeveloped fruits (%)")+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, 
                                   size=13, face = "italic"), 
        axis.title = element_text(size=18), 
        axis.text.y = element_text(size=15),
        legend.position = "none") +
  scale_x_discrete(labels=c("PBS_control" = "PBS\ncontrol", 
                            "A.citrulli_AAC001R" = "A. citrulli\nAAC001R",
                            "A.citrulli_M6R" = "A. citrulli\nM6R",
                            "Bacillus_T.m10.b1" = "Bacillus\nT.m10.b1R",
                            "Erwinia_P.s8.b2" = "Erwinia\nP.s8.b2R",
                            "Microbacterium_T.s9.b4" = "Microbacterium\nT.s9.b4R",
                            "Paraburkholderia_T.s8.b1" = "Paraburkholderia\nT.s8.b1R",
                            "Rosenbergiella_P.s10.b1" = "Rosenbergiella\nP.s10.b1R")) +
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#1F78B4","gray",
                               "#B2DF8A","gray","#33A02C","#FB9A99"))
fruit.bar

# chi-squared tests of these proportions
abort.chisq <- chisq.test(strain.abortion3)
abort.chisq # X-squared = 234.58, df = 7, p-value < 2.2e-16

# comparing each strain to the controls - Erwinia
control.erw.abort <- strain.abortion %>%
  filter(Bacterial_inoc. == "PBS_control"|Bacterial_inoc. == "Erwinia") %>%
  select(Bacterial_inoc., success_rate, failure_rate) %>%
  column_to_rownames(var = "Bacterial_inoc.")

erw.abort.chisq <- chisq.test(control.erw.abort)
erw.abort.chisq # X-squared = 0, df = 1, p-value = 1

# comparing each strain to the controls - A. citrulli AAC001R
control.AAC.abort <- strain.abortion %>%
  filter(Bacterial_inoc. == "PBS_control"|Bacterial_inoc. == "AcAAC001R") %>%
  select(Bacterial_inoc., success_rate, failure_rate) %>%
  column_to_rownames(var = "Bacterial_inoc.")

AAC.abort.chisq <- chisq.test(control.AAC.abort)
AAC.abort.chisq # X-squared = 5.0491, df = 1, p-value = 0.02464 *

# comparing each strain to the controls - A. citrulli M6R
control.M6R.abort <- strain.abortion %>%
  filter(Bacterial_inoc. == "PBS_control"|Bacterial_inoc. == "AcM6R") %>%
  select(Bacterial_inoc., success_rate, failure_rate) %>%
  column_to_rownames(var = "Bacterial_inoc.")

M6R.abort.chisq <- chisq.test(control.M6R.abort)
M6R.abort.chisq # X-squared = 10.501, df = 1, p-value = 0.001193 **

# comparing each strain to the controls - Bacillus
control.bac.abort <- strain.abortion %>%
  filter(Bacterial_inoc. == "PBS_control"|Bacterial_inoc. == "Bacillus") %>%
  select(Bacterial_inoc., success_rate, failure_rate) %>%
  column_to_rownames(var = "Bacterial_inoc.")

bac.abort.chisq <- chisq.test(control.bac.abort)
bac.abort.chisq # X-squared = 0.2192, df = 1, p-value = 0.6396

# comparing each strain to the controls - Paraburkholderia
control.para.abort <- strain.abortion %>%
  filter(Bacterial_inoc. == "PBS_control"|Bacterial_inoc. == "Paraburkholderia") %>%
  select(Bacterial_inoc., success_rate, failure_rate) %>%
  column_to_rownames(var = "Bacterial_inoc.")

para.abort.chisq <- chisq.test(control.para.abort)
para.abort.chisq # X-squared = 0, df = 1, p-value = 1

# comparing each strain to the controls - Rosenbergiella
control.rose.abort <- strain.abortion %>%
  filter(Bacterial_inoc. == "PBS_control"|Bacterial_inoc. == "Rosenbergiella") %>%
  select(Bacterial_inoc., success_rate, failure_rate) %>%
  column_to_rownames(var = "Bacterial_inoc.")

rose.abort.chisq <- chisq.test(control.rose.abort)
rose.abort.chisq # X-squared = 97.023, df = 1, p-value < 2.2e-16 *** 

# comparing each strain to the controls - Microbacterium
control.micro.abort <- strain.abortion %>%
  filter(Bacterial_inoc. == "PBS_control"|Bacterial_inoc. == "Microbacterium") %>%
  select(Bacterial_inoc., success_rate, failure_rate) %>%
  column_to_rownames(var = "Bacterial_inoc.")

micro.abort.chisq <- chisq.test(control.micro.abort)
micro.abort.chisq # X-squared = 97.023, df = 1, p-value < 2.2e-16 ***

#### Combining both graphs into Fig. 5A ####
Fig5A <- plot_grid(slice.heat, fruit.bar, ncol=1, align = "v")
Fig5A

