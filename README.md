# Bergmann et al. Contribution of floral transmission to the assembly of bacterial communities in watermelon seeds 
## Authors: Gillian E. Bergmann, Kacie Lui, Carina Lopez, Alex Velasco, Coralie Marais, Matthieu Barret, Marie Simonin, Rachel L. Vannette, Johan H.J. Leveau

This repository represents a combination of datasets and analyses used to produce the results for a manuscript on the floral transmission of bacterial communities to developing watermelon seeds. This manuscript will be deposited in biorxiv and submitted to New Phytologist for review. 
_______________________
This repository depends on the following R packages:
- tidyverse
- dada2
- phyloseq
- vegan
- ggpubr
- ANCOMBC
- microbiome
- microbiomeutilities
- UpSetR
- ComplexUpset

This repository has the following file structure:

- /scripts: folder with scripts to process and analyse each of the datasets used in the manuscript
- /figures: folder containing PDFs of each figure panel used in the manuscript
- /data_raw: folder with subfolders of raw datasets 
- /data_raw/Field.survey_NovaSeq: folder containing the raw sequence files for the field survey dataset (also available on NCBI Sequence Read Archive)
- /data_raw/Pollinator.Exp_MiSeq: folder containing the raw sequence files for the pollinator experiment dataset (also available on NCBI Sequence Read Archive)
- /data_raw/Seed_viability: folder of raw datasheets generated for the seed viability experiment
- /data_raw/Stigma_inoc: folder of raw datasheets generated for the stigma inoculation experiment
- /data_output_seq: folder containing subfolders of processed sequencing datasets
- /data_output_seq/Field_survey: folder containing phyloseq objects and other processed datasheets generated for the field survey dataset
- /data_output_seq/Pollinator_exp: folder containing phyloseq objects generated for the pollinator experiment dataset

_____________________
This project contains the following files: