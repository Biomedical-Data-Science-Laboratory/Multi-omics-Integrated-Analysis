library(dplyr)
library(tidyr)
library(cluster)
library(factoextra)
library(umap)
library(matrixStats)
library(survival)
library(ggplot2)
library(mixOmics)
library(enrichplot)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(msigdbr)

# Cox Regression Analysis

rna_mirna_protein_surv_obj <- Surv(time = merge_rna_mirna_protein_filter_BRCA_survival$OS.time, event = merge_rna_mirna_protein_filter_BRCA_survival$OS)

rna_mirna_protein_cox_model <- coxph(rna_mirna_protein_surv_obj ~
                     Factor1 + Factor2 + Factor3 + Factor4 + Factor5 +
                     Factor6 + Factor7 + Factor8 + Factor9 + Factor10 +
                     Factor11 + Factor12 + Factor13 + Factor14 + Factor15,
                   data = merge_rna_mirna_protein_filter_BRCA_survival)

summary(rna_mirna_protein_cox_model)
cox.zph(rna_mirna_protein_cox_model)

filter_merge_rna_mirna_protein_filter_BRCA_survival <- merge_rna_mirna_protein_filter_BRCA_survival %>%
  select(patient, Factor3, Factor6, Factor13)

rownames(filter_merge_rna_mirna_protein_filter_BRCA_survival) <- filter_merge_rna_mirna_protein_filter_BRCA_survival$patient
filter_merge_rna_mirna_protein_filter_BRCA_survival$patient <- NULL

# Cluster Analysis

## Optimal k Selection
