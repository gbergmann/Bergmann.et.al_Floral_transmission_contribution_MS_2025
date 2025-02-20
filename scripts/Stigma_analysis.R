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
slice.turbidity <- read.csv("data_raw/Stigma_inoc/")
fruit.abortion <- read.csv("data_raw/Stigma_inoc/")